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
    """ Create an instance for coordinate conversion between annoatation 
    versions. """
    def __init__(self, df, description=''):
        super().__init__(df, description) # use __init__() of parent class
    
    def get_dmel_coordinates(self, query_version, ref_version, query):
        # Parse input query
        query_coord = self.query_parser(query, query_version)
        
        # Find query coordinate in DataFrame
        filt_kw = {
            f'v{query_version}_chr': query_coord.chromosome,
            f'v{query_version}_start': f'lte{query_coord.start}',
            f'v{query_version}_end': f'gte{query_coord.end}',
        }
        # Assumes only one pair will hit
        conv_dict = self.filter(**filt_kw).iloc[0].to_dict()
        
        # Get position of the query coordinates
        # Assumes continuous coordinates in both versions
        ix1 = query_coord.start - conv_dict[f'v{query_version}_start']
        ix2 = query_coord.end - conv_dict[f'v{query_version}_end']
        
        # Refer to the other version
        s2 = range(conv_dict[f'v{ref_version}_start'], conv_dict[f'v{ref_version}_end'])[ix1]
        e2 = range(conv_dict[f'v{ref_version}_start'], conv_dict[f'v{ref_version}_end'])[ix2]
        
        return GenomicCoordinate(
            conv_dict[f'v{ref_version}_chr'], s2, e2, ref_version)
    
    def query_parser(self, query, version):
        chrname = query.split(':')[0]
        start = int(query.split(':')[1].split('..')[0])
        end = int(query.split(':')[1].split('..')[1])
        
        return GenomicCoordinate(chrname, start, end, version)
