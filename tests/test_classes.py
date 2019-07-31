""" Nose tests for classes GenomicCoordinate, PairedGenomicCoordinate, 
ConvertCoordinates and their instances. """

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
	""" Unit tests for PairedGenomicCoordinate """
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


