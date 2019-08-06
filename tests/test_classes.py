""" Nose tests for classes GenomicRange, PairedGenomicRanges, 
ConvertCoordinates and their instances. """
from nose.tools import assert_raises
import pandas as pd
from convert_annotation.classes import GenomicRange, PairedGenomicRanges, \
	ConvertCoordinates

class TestGenomicRange:
    """ Unit tests for GenomicRange. """
    def setup(self):
        self.gencood = GenomicRange('2L', 1, 10)

    def teardown(self):
        pass

    def test_eq(self):
        assert self.gencood == GenomicRange('2L', 1, 10)
        assert self.gencood != ('2L', 1, 10)
        assert self.gencood != GenomicRange('3L', 1, 10)
        assert self.gencood != GenomicRange('2L', 2, 10)
        assert self.gencood != GenomicRange('2L', 1, 20)

    def test_chromosome(self):
        assert self.gencood.chromosome == '2L'

    def test_start(self):
        assert self.gencood.start == 1

    def test_end(self):
        assert self.gencood.end == 10

    def test_from_str(self):
    	assert self.gencood.from_str('3R:100..200') == \
            GenomicRange('3R', 100, 200)

    def test_get_index(self):
        assert self.gencood.get_index(3) == 2
        assert_raises(IndexError, self.gencood.get_index, -2)
        assert_raises(IndexError, self.gencood.get_index, 11)
        assert_raises(IndexError, self.gencood.get_index, 0)

class TestPairedGenomicRanges:
    """ Unit tests for PairedGenomicRanges. """
    def setup(self):
        self.paired_gencoord = PairedGenomicRanges(
            GenomicRange('2L', 1, 10),
            GenomicRange('2L', 11, 20)
        )

    def teardown(self):
        pass

    def test_eq(self):
        assert self.paired_gencoord == PairedGenomicRanges(
            GenomicRange('2L', 1, 10),
            GenomicRange('2L', 11, 20)
        )
        assert self.paired_gencoord != (('2L', 1, 10), ('2L', 11, 20))
        assert self.paired_gencoord != PairedGenomicRanges(
            GenomicRange('3L', 1, 10),
            GenomicRange('2L', 11, 20)
        )
        assert self.paired_gencoord != PairedGenomicRanges(
            GenomicRange('2L', 1, 10),
            GenomicRange('2L', 101, 110)
        )

    def test_query(self):
        assert self.paired_gencoord.query.chromosome == '2L'
        assert self.paired_gencoord.query.start == 1
        assert self.paired_gencoord.query.end == 10

    def test_match(self):
        assert self.paired_gencoord.match.chromosome == '2L'
        assert self.paired_gencoord.match.start == 11
        assert self.paired_gencoord.match.end == 20

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

    def test_get_dmel_coordinates(self):
        # Test 5 to 6 conversion
        assert self.cc.get_dmel_coordinates(5, 6, '2L:1..10') == \
            PairedGenomicRanges(
                GenomicRange(chromosome='2L', start=1, end=10),
                GenomicRange(chromosome='2L', start=11, end=20)
            )

        # Another way to input query coordinates: GenomicRange object
        assert self.cc.get_dmel_coordinates(
            5, 6, 
            GenomicRange(chromosome='2L', start=1, end=10)) == \
            PairedGenomicRanges(
                GenomicRange(chromosome='2L', start=1, end=10),
                GenomicRange(chromosome='2L', start=11, end=20)
            )

        # Another way to input query coordinates: specify each property
        assert self.cc.get_dmel_coordinates(
            5, 6, query_chr='2L', query_start=1, query_end=10) == \
            PairedGenomicRanges(
                GenomicRange(chromosome='2L', start=1, end=10),
                GenomicRange(chromosome='2L', start=11, end=20)
            )

        # Test 6 to 5 conversion
        assert self.cc.get_dmel_coordinates(6, 5, '2L:11..15') == \
            PairedGenomicRanges(
                GenomicRange(chromosome='2L', start=11, end=15),
                GenomicRange(chromosome='2L', start=1, end=5)
            )

        # Test chromosoem changes
        assert self.cc.get_dmel_coordinates(5, 6, '2L:23..27') == \
            PairedGenomicRanges(
                GenomicRange(chromosome='2L', start=23, end=27),
                GenomicRange(chromosome='2R', start=24, end=28)
            )
        
        # Test duplication in one version
        assert self.cc.get_dmel_coordinates(6, 5, '2R:25..30') == \
            PairedGenomicRanges(
                GenomicRange(chromosome='2R', start=25, end=30),
                GenomicRange(chromosome=-8, start=-8, end=-8)
            )

        # Test the case where query is not found in the list
        assert self.cc.get_dmel_coordinates(5, 6, '3L:100..200') == \
            PairedGenomicRanges(
                GenomicRange(chromosome='3L', start=100, end=200),
                GenomicRange(chromosome=-9, start=-9, end=-9)
            ), \
            'Found {}'.format(self.cc.get_dmel_coordinates(5, 6, '3L:100..200'))

