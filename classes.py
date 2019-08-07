import re
import pandas as pd
from collections import namedtuple
from collections.abc import Mapping

class GenomicRange(object):
    """ Construct GenomicRange object. This GenomicRange is for chromosomal 
    position data. Coordinates are counted from 1 (e.g. the first nucleotide 
    at a given chromosome is counted as 1).
    """
    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self._range = range(start, end+1)

    @staticmethod
    def from_str(gencoord_str):
        """ Returns GenomicRange object for input string.
        Parameter   
        ---------
        gencoord_str: str
            gencoord_str has to be in the format of 
            "{chr}:{start}..{end}:{version}"
        
        Return
        ------
        GenomicRange object
        """
        chrname = gencoord_str.split(':')[0]
        start = int(gencoord_str.split(':')[1].split('..')[0])
        end = int(gencoord_str.split(':')[1].split('..')[1])

        return GenomicRange(chrname, start, end)

    def get_index(self, coordinate):
        """ Returns an index of a given coordinate within GenomicRange. For 
        example, an index of the first coordinate is 0.
        """
        if coordinate < 0:
            raise IndexError('Invalid coordinate. Negative value is not '\
                'available in GenomicRange.')
        elif coordinate > self.end:
            raise IndexError(f'list index out of range: {coordinate} > {self.end}.')
        elif coordinate < self.start:
            raise IndexError(f'list index out of range: {coordinate} < {self.start}.')

        return coordinate - self.start

    def __len__(self):
        return len(self._range)
    
    def __getitem__(self, index):
        return self._range[index]

    def __eq__(self, compare):
        if isinstance(compare, GenomicRange):
            if self.chromosome != compare.chromosome:
                return False
            elif self.start != compare.start:
                return False
            elif self.end != compare.end:
                return False
            # If compare object passed all comparisons, return True
            return True

        return False

    def __iter__(self):
        return self._range.__iter__()

    def __repr__(self):
        return '{name}(chromosome={chromosome}, start={start}, end={end})'\
            .format(
                name=type(self).__name__, chromosome=self.chromosome, 
                start=self.start, end=self.end
            )

class PairedGenomicRanges(Mapping):
    def __init__(self, keys, ranges, is_inversion, name=None):
        self._keys = tuple(keys)
        self._ranges = tuple(ranges)
        self._pair = dict(zip(keys, ranges))
        self.is_inversion = is_inversion
        self.name = name

    @property
    def keys(self):
        return self._keys

    @property
    def ranges(self):
        return self._ranges

    @staticmethod
    def from_dict(pairedRangeDict, is_inversion):
        keys = pairedRangeDict.keys()
        ranges = pairedRangeDict.values()
        
        return PairedGenomicRanges(keys, ranges, is_inversion)

    def get_another_key(self, key1):
        """ Assumes that only two GenomicRange objects are coupled in 
        this object. """
        ix1 = self.keys.index(key1)
        ix2 = 0 if ix1 == 1 else 1

        return self.keys[ix2]

    def convert_range(self, query_key, query_genomicRange):
        # Get position of the query coordinates
        # Assumes continuous coordinates in both versions
        start_ix = self[query_key].get_index(query_genomicRange.start)
        end_ix = self[query_key].get_index(query_genomicRange.end)

        match_key = self.get_another_key(query_key)

        if self.is_inversion:
            start_ix, end_ix = -(end_ix+1), -(start_ix+1)

        try:
            s2 = self[match_key]._range[start_ix]
        except IndexError:
            raise Exception('IndexError found: {} '\
                'while length is {}'.format(start_ix, len(self[match_key])))
        try:
            e2 = self[match_key]._range[end_ix]
        except IndexError:
            raise Exception('IndexError found: {} '\
                'while length is {}'.format(end_ix, len(self[match_key])))
        
        return PairedGenomicRanges(
            keys=[query_key, match_key],
            ranges=[
                query_genomicRange, 
                GenomicRange(self[match_key].chromosome, s2, e2)],
            is_inversion=self.is_inversion, 
            name=self.name
        )

    def to_Series(self):
        indices = []
        data = []

        for key, genrange in self:
            data += [genrange.chromosome, genrange.start, genrange.end]
            indices = [f'v{key}_chromosome', f'v{key}_start', f'v{key}_end']

        data += [self.is_inversion]
        indices += ['strand']

        return pd.Series(data, index=indices, name=self.name)
    
    def __len__(self):
        return len(self._pair)

    def __iter__(self):
        for item in self._pair.items():
            yield item

    def __getitem__(self, key):
        return self._pair[key]

    def __eq__(self, compare):
        if isinstance(compare, PairedGenomicRanges):
            if self.keys != compare.keys:
                return False
            if self.ranges != compare.ranges:
                return False
            if self.is_inversion != compare.is_inversion:
                return False
            # If compare object passed all comparisons, return True
            return True

        return False

    def __repr__(self):
        return '{pair}; is_inversion={is_inversion}'.format(
            pair=self._pair.__repr__(), is_inversion=self.is_inversion
        )

class Database(Mapping):
    """ This class inherits Mapping class. __iter__, __getitem__ and __len__ 
    functions are overwritten. This is a base class of SFS class. """
    def __init__(self, df, description=''):
        self.df = df
        self.description = description

    @property
    def columns(self):
        return self.df.columns

    @property
    def shape(self):
        return self.df.shape
        
    def filter(self, sort_by='', ascending=True, **kwargs):
        '''
        Search rows which have specifies items from a given dataframe.
        Please pass key words for searching to **kwargs.
        For example, if you want to get items that is greater than equal (>=)
        100 in column "A", please specify **kwargs as "A=gte100". 
        Please see below for details.
        If nothing passed to **kwargs, return input dataframe.

        Paramters
        ---------
        sort_by: str
            column name which in output dataframe is sorted by values in.
        ascending: bool (default: True)
        **kwargs:
            key is for column, value is for filtering values (items)
            You can use indicators below for filtering way.
            "gt" for ">"
            "gte" for ">="
            "lt" for "<"
            "lte" for "<="
            "ne" for "!="
            "c/" for "contains"
            "" for "=="
            If you pass tuple to value, this function search and filter 
            items recursively.

        Dependencies
        ------------
        pandas
        re
        '''
        res_df = self.df
        def f(res_df, k, v):
            if isinstance(v, str):
                if v == '*':
                    pass
                elif re.search('^gt\d+', v):
                    v = float(re.search('^gt(\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] > v]
                elif re.search('^gte\d+', v):
                    v = float(re.search('^gte(\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] >= v]
                elif re.search('^lt\d+', v):
                    v = float(re.search('^lt(\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] < v]
                elif re.search('^lte\d+', v):
                    v = float(re.search('lte(\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] <= v]
                elif re.search('^ne\d+', v):
                    v = float(re.search('ne(\d+\.*\d*)$', v).group(1))
                    res_df = res_df[res_df[k] != v]
                elif re.search('^c\/', v):
                    v = re.search('^c\/(.+)\/$', v).group(1)
                    res_df = res_df[res_df[k].str.contains(v)]
                elif re.search('^nc\/', v):
                    v = re.search('^nc\/(.+)\/$', v).group(1)
                    res_df = res_df[~res_df[k].str.contains(v)]
                else:
                    res_df = res_df[res_df[k] == v]
            elif isinstance(v, list):
                res_df = res_df[res_df[k].isin(v)]
            else:
                res_df = res_df[res_df[k] == v]

            return res_df

        for k,v in kwargs.items():
            if isinstance(v, tuple): # "and"
                res_df_list = []
                for i in v:
                    tmp_res_df = f(res_df, k, i)
                    res_df = pd.merge(res_df, tmp_res_df, how='inner')
            elif isinstance(v, (str, list)): # list means "or" condition
                res_df = f(res_df, k, v)
            elif isinstance(v, int):
                res_df = res_df[res_df[k] == v]
        if sort_by:
            res_df.sort_values(by=sort_by, ascending=ascending, inplace=True)
            
        return res_df

    def groupby(self, **kwargs):
        return self.df.groupby(**kwargs)

    def head(self, *kwargs):
        return self.df.head(*kwargs)

    def tail(self, *kwargs):
        return self.df.tail(*kwargs)

    def __len__(self):
        return len(self.df.index)

    def __iter__(self):
        return self.df.iterrows()

    def __getitem__(self, key):
        if key == '*':
            return self.df
        else:
            return self.df.loc[key, :]

    def __repr__(self):
        return '<{name}: {desc} ({size} records)>'.format(
            name=type(self).__name__, desc=self.description, size=self.__len__()
        )

class ConvertCoordinates(Database):
    def __init__(self, df, version1, version2, description=''):
        super().__init__(df, description) # use __init__() of parent class
        self.version1, self.version2 = version1, version2

        # Check if all columns are contined in DataFrame
        self.check_columns()
        self.check_same_len()

    def check_columns(self):
        assert f'v{self.version1}_chr' in self.df.columns, \
            f'Column not found: "v{self.version1}_chr" was expected.'
        assert f'v{self.version1}_start' in self.df.columns, \
            f'Column not found: "v{self.version1}_start" was expected.'
        assert f'v{self.version1}_end' in self.df.columns, \
            f'Column not found: "v{self.version1}_end" was expected.'

        assert f'v{self.version2}_chr' in self.df.columns, \
            f'Column not found: "v{self.version2}_chr" was expected.'
        assert f'v{self.version2}_start' in self.df.columns, \
            f'Column not found: "v{self.version2}_start" was expected.'
        assert f'v{self.version2}_end' in self.df.columns, \
            f'Column not found: "v{self.version2}_end" was expected.'

    def check_same_len(self):
        v1_start = f'v{self.version1}_start'
        v1_end = f'v{self.version1}_end'
        v2_start = f'v{self.version2}_start'
        v2_end = f'v{self.version2}_end'

        assert self.df.apply(lambda x: x[v1_end] - x[v1_start] ==\
                x[v2_end] - x[v2_start], axis=1).all()

    def get_rows(self, version, chromosome, start, end, sort_by='', ascending=True):
        """ Returns a subset of DataFrame where a segment include input range 
        completely. Currently this function does not take care of partial matches.
        """
        filt_kw = {
            f'v{version}_chr': chromosome,
            f'v{version}_start': f'lte{start}',
            f'v{version}_end': f'gte{end}',
        }
        return self.filter(sort_by, ascending, **filt_kw)

    def get_another_version(self, version):
        if version == self.version1:
            return self.version2
        elif version == self.version2:
            return self.version1
        else:
            raise Exception(f'Please input either of version {self.version1} '\
                'or {self.version2}.')

    @staticmethod
    def to_PairedGenomicRanges(row, version1, version2, name=None):
        keys = [version1, version2]
        ranges = [
            GenomicRange(
                row[f'v{version1}_chr'], 
                row[f'v{version1}_start'],
                row[f'v{version1}_end']
            ),
            GenomicRange(
                row[f'v{version2}_chr'], 
                row[f'v{version2}_start'],
                row[f'v{version2}_end']
            )
        ]
        is_inversion = True if row['strand'] == '-' else False

        return PairedGenomicRanges(keys, ranges, is_inversion, name)

    def to_list_of_PairedGenomicRanges(self):
        return self.df.apply(lambda x: 
            self.to_PairedGenomicRanges(x, self.version1, self.version2), 
            axis=1).tolist()

    def convert_coordinate(self, query_version, query=None,
                           query_chr='', query_start=None, query_end=None):
        # Parse input query
        if isinstance(query, str):
            query_coord = GenomicRange.from_str(query)
        elif isinstance(query, GenomicRange):
            query_coord = query
        elif query_chr and query_start and query_end:
            query_coord = GenomicRange(
                query_chr, query_start, query_end)
        else:
            raise Exception('Please input query or query_chr, query_start and'\
                ' query_end.')
        
        # Find query coordinate in DataFrame
        conv_table = self.get_rows(
            query_version, query_coord.chromosome, 
            query_coord.start, query_coord.end)

        # Get another version
        match_version = self.get_another_version(query_version)
        
        # If query is not found, return -9
        if len(conv_table) == 0:
            return PairedGenomicRanges(
                keys=[query_version, match_version], 
                ranges=[query_coord, GenomicRange(-9, -9, -9)], 
                is_inversion=-9
            )
            
        # If query range is found in the multiple rows, return -8
        elif len(conv_table) > 1:
            return PairedGenomicRanges(
                keys=[query_version, match_version], 
                ranges=[query_coord, GenomicRange(-8, -8, -8)], 
                is_inversion=-9
            )

        paired = self.to_PairedGenomicRanges(
            conv_table.iloc[0], self.version1, self.version2, conv_table.index)

        # Return converted coordinates
        return paired.convert_range(query_version, query_coord)

    @staticmethod
    def from_query_strs_to_query_coords(query_strs):
        return list(map(GenomicRange.from_str, query_strs))
        
    def convert_coordinates(self, query_version, query_coords):
        for query_coord in query_coords:
            yield self.convert_coordinate(query_version, query_coord)

    def from_querys_to_DataFrame(self, query_version, query_strs):
        query_coords = self.from_query_strs_to_query_coords(query_strs)

        concat_list = [
            pair.to_Series() for parir in 
            self.convert_coordinate(query_version, query_coords)
        ]

        return pd.concat(concat_list)
                
    def __repr__(self):
        return '<{name}: {desc} (versions {v1} and {v2}; {size} records)>'.format(
            name=type(self).__name__, desc=self.description, 
            v1=self.version1, v2=self.version2,
            size=self.__len__()
        )
