import pandas as pd
from convert_annotation.classes import ConvertCoordinates

def from_CSV_to_ConvertCoordinates(csv_path, version1, version2, description=''):
    """ Read CSV file into ConvertCoordinates object. """
    return ConvertCoordinates(
		pd.read_csv(csv_path), version1, version2, description)

def parse_query_list(path, expect='itemnum: ', avoid=['itemnum', '/*']):
	""" Returns a list of query names that are written in a given file.
	Parameter
	---------
	path: str
		a path to a file of file name list
	
	Return
	------
	query_list: list
		a list that contains all file names listed in a given file. 
	"""
	query_list = []
	with open(path, 'r') as f:
		for line in f:
			if expect:
				if line.startswith(expect):
					exp_itemnum = int(line.rstrip().split(expect)[1])
            
			bad = 0
			for bad_char in avoid:
				if line.startswith(bad_char):
					bad += 1
			if bad > 0:
				continue
			
			query_list.append(line.rstrip())
	if expect:
		if exp_itemnum:
			assert len(query_list) == exp_itemnum, \
				f'Wrong item number: {len(query_list)} observed instead of '\
					'{exp_itemnum}.'
		else:
			raise Exception('Expected number has not been parsed.')
	
	return query_list
