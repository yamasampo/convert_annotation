import re
import pandas as pd
from collections import namedtuple
from collections.abc import Mapping

class GenomicRange(object):
	def __init__(self, chromosome, start, end):
		self.chromosome = chromosome
		self.start = start
		self.end = end

	@staticmethod
	def from_str(gencoord_str):
		""" Returns GenomicRange object for input string.
		Parameter   
		---------
		gencoord_str: str
			gencoord_str has to be in the format of 
			"{chr}:{start}..{start}:{version}"
		
		Return
		------
		GenomicRange object
		"""
		chrname = gencoord_str.split(':')[0]
		start = int(gencoord_str.split(':')[1].split('..')[0])
		end = int(gencoord_str.split(':')[1].split('..')[1])

		return GenomicRange(chrname, start, end)

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

	def __repr__(self):
		return '{name}(chromosome={chromosome}, start={start}, end={end})'\
			.format(
				name=type(self).__name__, chromosome=self.chromosome, 
				start=self.start, end=self.end
			)

class PairedGenomicRanges(object):
	def __init__(self, query, match, description=''):
		self.query = query
		self.match = match
		self.description = description

	def __eq__(self, compare):
		if isinstance(compare, PairedGenomicRanges):
			if self.query != compare.query:
				return False
			elif self.match != compare.match:
				return False
			# If compare object passed all comparisons, return True
			return True

		return False


	def __repr__(self):
		if self.description:
			return '{name}: {desc} (query={query}, match={match})'\
				.format(
					name=type(self).__name__, 
					query=self.query, match=self.match, desc=self.description
				)

		else:
			return '{name}(query={query}, match={match})'\
				.format(
					name=type(self).__name__, 
					query=self.query, match=self.match
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
    
	def get_dmel_coordinates(self, query_version, ref_version, query=None,
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
		
		# If query is not found, return -9
		if len(conv_table) == 0:
			return PairedGenomicRanges(
				query_coord, 
				GenomicRange(-9, -9, -9)
			)
		# If query range is found in the multiple rows, return -8
		elif len(conv_table) > 1:
			return PairedGenomicRanges(
				query_coord, 
				GenomicRange(-8, -8, -8)
			)

		conv_dict = conv_table.iloc[0].to_dict()
		# Get position of the query coordinates
		# Assumes continuous coordinates in both versions
		ix1 = query_coord.start - conv_dict[f'v{query_version}_start']
		ix2 = query_coord.end - conv_dict[f'v{query_version}_start']

		# Refer to the other version
		ref_range = range(conv_dict[f'v{ref_version}_start'], 
						  conv_dict[f'v{ref_version}_end']+1)
		try:
			s2 = ref_range[ix1]
		except IndexError:
			raise Exception('IndexError found: {} '\
				'while length is {}'.format(ix1, len(ref_range)))
		try:
			e2 = ref_range[ix2]
		except IndexError:
			raise Exception('IndexError found: {} '\
				'while length is {}'.format(ix2, len(ref_range)))

		result_coord = GenomicRange(
			conv_dict[f'v{ref_version}_chr'], s2, e2)
			
		return PairedGenomicRanges(query_coord, result_coord)

	def recursively_get_dmel_coordinates(self, query_version, ref_version, querys):
		query_coords = [self.query_parser(q, query_version) for q in querys]
		result_pairs = []

		for query_coord in query_coords:
			result = self.get_dmel_coordinates(
				query_version, ref_version, query_coord)
			result_pairs.append(result)
		
		return result_pairs
		
	def __repr__(self):
		return '<{name}: {desc} (versions {v1} and {v2}; {size} records)>'.format(
			name=type(self).__name__, desc=self.description, 
			v1=self.version1, v2=self.version2,
			size=self.__len__()
		)
