import re
import pandas as pd
from collections import namedtuple
from collections.abc import Mapping

class Database(Mapping):
    """ This class inherits Mapping class. __iter__, __getitem__ and __len__ 
    functions are overwritten. This is a base class of SFS class. """
    def __init__(self, df, description=''):
        self.df = df
        self.description = description
        
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
            df: DataFrame (pandas)
                input dataframe
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
    
GenomicCoordinate = namedtuple(
    'GenomicCoordinate', ['chromosome', 'start', 'end', 'version']
)
class ConvertCoordinates(Database):
    def __init__(self, df, version1, version2, description=''):
        super().__init__(df, description) # use __init__() of parent class
        self.version1, self.version2 = version1, version2

        # Check if all columns are contined in DataFrame
        self.check_columns()
        
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
    
    def get_dmel_coordinates(self, query_version, ref_version, query):
        # Parse input query
        if isinstance(query, str):
            query_coord = self.query_parser(query, query_version)
        elif isinstance(query, GenomicCoordinate):
            query_coord = query
        
        # Find query coordinate in DataFrame
        filt_kw = {
            f'v{query_version}_chr': query_coord.chromosome,
            f'v{query_version}_start': f'lte{query_coord.start}',
            f'v{query_version}_end': f'gte{query_coord.end}',
        }
        conv_dicts = self.filter(**filt_kw).T.to_dict()
        
        # If query is not found, return -9
        if len(conv_dicts) == 0:
            return GenomicCoordinate('', -9, -9, -9)
        
        # Collect results
        results = []
        
        for i, conv_dict in conv_dicts.items():
            # Get position of the query coordinates
            # Assumes continuous coordinates in both versions
            ix1 = query_coord.start - conv_dict[f'v{query_version}_start']
            ix2 = query_coord.end - conv_dict[f'v{query_version}_start']

            # Refer to the other version
            ref_range = range(conv_dict[f'v{ref_version}_start'], 
                              conv_dict[f'v{ref_version}_end']+1)
            s2 = ref_range[ix1]
            e2 = ref_range[ix2]
            
            result_coord = GenomicCoordinate(
                conv_dict[f'v{ref_version}_chr'], s2, e2, ref_version)
            
            results.append(result_coord)
            
        if len(results) == 1:
            return query_coord, results[0]
        
        else:
            return query_coord, results

    def recursively_get_dmel_coordinates(self, query_version, ref_version, querys):
        query_coords = [self.query_parser(q) for q in querys]
        out_d = dict.fromkeys(query_coords)

        for query_coord in query_coords:
            _, result = self.get_dmel_coordinates(
                query_version, ref_version, query_coord)
            out_d[query_coord] = result

        return out_d
    
    def query_parser(self, query, version):
        chrname = query.split(':')[0]
        start = int(query.split(':')[1].split('..')[0])
        end = int(query.split(':')[1].split('..')[1])
        
        return GenomicCoordinate(chrname, start, end, version)
    
    def __repr__(self):
        return '<{name}: {desc} (versions {v1} and {v2}; {size} records)>'.format(
            name=type(self).__name__, desc=self.description, 
            v1=self.version1, v2=self.version2,
            size=self.__len__()
        )
