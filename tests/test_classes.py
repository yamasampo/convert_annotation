""" Nose tests for classes GenomicCoordinate, PairedGenomicCoordinate, 
ConvertCoordinates and their instances. """

import pandas as pd
from convert_annotation.classes import GenomicCoordinate, PairedGenomicCoord, \
	ConvertCoordinates

class TestGenomicCoordinate:
	""" Unit tests for GenomicCoordinate. """
	def setup(self):
		self.gencood = GenomicCoordinate('2L', 1, 10, 5)

	def teardown(self):
		pass

	def test_chromosome(self):
		assert self.gencood.chromosome == '2L'

	def test_start(self):
		assert self.gencood.start == 1

	def test_end(self):
		assert self.gencood.end == 10

	def test_version(self):
		assert self.gencood.version == 5

class TestPairedGenomicCoordinate:
	""" Unit tests for PairedGenomicCoordinate. """
	def setup(self):
		self.paired_gencoord = PairedGenomicCoord(
			GenomicCoordinate('2L', 1, 10, 5),
			GenomicCoordinate('2L', 11, 20, 6)
		)
	
	def teardown(self):
		pass
	
	def test_query(self):
		assert self.paired_gencoord.query.chromosome == '2L'
		assert self.paired_gencoord.query.start == 1
		assert self.paired_gencoord.query.end == 10
		assert self.paired_gencoord.query.version == 5

	def test_reference(self):
		assert self.paired_gencoord.reference.chromosome == '2L'
		assert self.paired_gencoord.reference.start == 11
		assert self.paired_gencoord.reference.end == 20
		assert self.paired_gencoord.reference.version == 6

class TestConvertCoordinates:
	""" Unit tests for ConvertCoordinates. """
	def setup(self):
		df = pd.DataFrame(
			{
				'v5_chr': ['2L', '2L', '2L'], 
				'v5_start': [1, 20, 35], 
				'v5_end': [10, 30, 45], 
				'v6_chr': ['2L', '2R', '2R'], 
				'v6_start': [11, 21, 25], 
				'v6_end': [20, 31, 35], 
				'strand': ['+', '+', '+']
			}
		)
		self.cc = ConvertCoordinates(
			df, 5, 6, 'Test ConvertCoordinates object')
		
	def teardown(self):
		pass

	def test_query_parser(self):
		assert self.cc.query_parser('3R:100..200', 5) == \
			GenomicCoordinate('3R', 100, 200, 5)

	def test_get_dmel_coordinates(self):
		# Test 5 to 6 conversion
		assert self.cc.get_dmel_coordinates(5, 6, '2L:1..10') == \
			PairedGenomicCoord(
				GenomicCoordinate(chromosome='2L', start=1, end=10, version=5),
				[GenomicCoordinate(chromosome='2L', start=11, end=20, version=6)]
			)

		# Test 6 to 5 conversion
		assert self.cc.get_dmel_coordinates(6, 5, '2L:11..15') == \
			PairedGenomicCoord(
				GenomicCoordinate(chromosome='2L', start=11, end=15, version=6),
				[GenomicCoordinate(chromosome='2L', start=1, end=5, version=5)]
			)

		# Test chromosoem changes
		assert self.cc.get_dmel_coordinates(5, 6, '2L:23..27') == \
			PairedGenomicCoord(
				GenomicCoordinate(chromosome='2L', start=23, end=27, version=5),
				[GenomicCoordinate(chromosome='2R', start=24, end=28, version=6)]
			)
		
		# Test duplication in one version
		assert self.cc.get_dmel_coordinates(6, 5, '2R:25..30') == \
			PairedGenomicCoord(
				GenomicCoordinate(chromosome='2R', start=25, end=30, version=6),
				[GenomicCoordinate(chromosome='2L', start=24, end=29, version=5),
				 GenomicCoordinate(chromosome='2L', start=35, end=40, version=5)]
			)

		# Test the case where query is not found in the list
		assert self.cc.get_dmel_coordinates(5, 6, '3L:100..200') == \
			PairedGenomicCoord(
				GenomicCoordinate(chromosome='3L', start=100, end=200, version=5),
				[GenomicCoordinate(chromosome='', start=-9, end=-9, version=-9)]
			), \
			'Found {}'.format(self.cc.get_dmel_coordinates(5, 6, '3L:100..200'))

