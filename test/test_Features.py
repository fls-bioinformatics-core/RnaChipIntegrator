#
#     test_Features.py: unit tests for Features module
#     Copyright (C) University of Manchester 2011-5 Peter Briggs

from common import *
from rnachipintegrator.Features import Feature
from rnachipintegrator.Features import FeatureSet
import unittest

class TestFeature(unittest.TestCase):

    def setUp(self):
        # Set up Feature objects to be used in the tests
        # Forward strand example
        self.rna_data = \
            Feature('CG9130-RB','chr3L','1252012','1255989','+')
        # Reverse strand example
        self.rna_data_2 = \
            Feature('CG13051-RA','chr3L','16257914','16258166','-')

    def test_properties(self):
        self.assertEqual(self.rna_data.chrom,'chr3L')
        self.assertEqual(self.rna_data_2.chrom,'chr3L')

    def test__eq__(self):
        self.assertEqual(self.rna_data,Feature('CG9130-RB',
                                               'chr3L',
                                               '1252012',
                                               '1255989','+'))
        self.assertNotEqual(self.rna_data,self.rna_data_2)

    def test_contains_position(self):
        position = 1253000
        self.assertTrue(self.rna_data.containsPosition(position),
                        "Transcript should contain position")
        position = 4250000
        self.assertFalse(self.rna_data.containsPosition(position),
                         "Transcript should not contain position")
        position = 10000
        self.assertFalse(self.rna_data.containsPosition(position),
                         "Transcript should not contain position")

    def test_get_closest_edge_distance(self):
        # Single position
        position = 1200000
        # Single reference position
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position),
                         abs(self.rna_data.start-position),
                         "Incorrect closest edge distance #1")
        position = 1300000
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position),
                         abs(self.rna_data.end-position),
                         "Incorrect closest edge distance #2")

    def test_get_closest_edge_distance_outside_region(self):
        # Reference region (2 positions)
        position1 = 1100000
        position2 = 1200000
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2),
                         abs(self.rna_data.start-position2),
                         "Incorrect closest edge distance (region #1)")
        position1 = 1300000
        position2 = 1400000
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2),
                         abs(self.rna_data.end-position1),
                         "Incorrect closest edge distance (region #2)")

    def test_get_closest_edge_distance_partially_inside_region(self):
        # Partially inside reference region
        position1 = 1200000
        position2 = 1255000
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2),
                         abs(self.rna_data.end-position2),
                         "Incorrect closest edge distance (inside region #1)")
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2,
                                                                zero_inside_region=True),
                         0,
                         "Incorrect closest edge distance (inside region #2)")

    def test_get_closest_edge_distance_completely_inside_region(self):
        # Completely inside reference region
        position1 = 1250000
        position2 = 1300000
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2),
                         abs(self.rna_data.start-position1),
                         "Incorrect closest edge distance (inside region #3)")
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2,
                                                                zero_inside_region=True),
                         0,
                         "Incorrect closest edge distance (inside region #4)")

    def test_get_promoter_region(self):
        leading = 10000
        trailing = 2500
        promoter = self.rna_data.getPromoterRegion(leading,trailing)
        self.assertEqual(promoter,
                         (self.rna_data.getTSS()-leading,
                          self.rna_data.getTSS()+trailing),
                         "Incorrect promoter region for + strand example")
        promoter = self.rna_data_2.getPromoterRegion(leading,trailing)
        self.assertEqual(promoter,
                         (self.rna_data_2.getTSS()+leading,
                          self.rna_data_2.getTSS()-trailing),
                         "Incorrect promoter region for - strand example")

class TestFeatureSet(unittest.TestCase):

    def setUp(self):
        # Create input files for tests
        create_test_file('Transcripts-ex1.txt',transcripts_ex1)
        create_test_file('Transcripts-ex2.txt',transcripts_ex2)
        create_test_file('Transcripts-ex2a.txt',transcripts_ex2a)

    def tearDown(self):
        # Remove input files
        delete_test_file('Transcripts-ex1.txt')
        delete_test_file('Transcripts-ex2.txt')

    def test_populate_from_list_of_features(self):
        features = FeatureSet(
            features_list=(
                Feature('CG31973','chr2L',25402,59243,'-'),
                Feature('CG2674-RC','chr2L',107926,114433,'+'),
                Feature('CG2674-RE','chr2L',106903,114433,'+'),
                Feature('CG2674-RA','chr2L',107760,114433,'+')))
        self.assertEqual(features[0],
                         Feature('CG31973','chr2L',25402,59243,'-'))
        self.assertEqual(features[1],
                         Feature('CG2674-RC','chr2L',107926,114433,'+'))
        self.assertEqual(features[2],
                         Feature('CG2674-RE','chr2L',106903,114433,'+'))
        self.assertEqual(features[3],
                         Feature('CG2674-RA','chr2L',107760,114433,'+'))

    def test_reading_in_RNAseq_data(self):
        rna_seq = FeatureSet('Transcripts-ex1.txt')
        self.assertEqual(len(rna_seq),10,
                         "Wrong number of lines from RNA-seq file")
        self.assertTrue(rna_seq.isFlagged(),
                        "Data should be flagged")

    def test_filter_on_chromosome(self):
        rna_seq = FeatureSet('Transcripts-ex1.txt')
        rna_chr = rna_seq.filterByChr('chr3LHet')
        self.assertEqual(len(rna_chr),1,
                         "Wrong number of chromosomes")
        for rna_data in rna_chr:
            self.assertEqual(rna_data.chrom,'chr3LHet',
                             "Wrong chromosome filtered")

    def test_filter_on_strand(self):
        rna_seq = FeatureSet('Transcripts-ex1.txt')
        rna_plus = rna_seq.filterByStrand('+')
        self.assertEqual(len(rna_plus),5,
                         "Wrong number of + strands")
        rna_minus = rna_seq.filterByStrand('-')
        self.assertEqual(len(rna_minus),5,
                         "Wrong number of - strands")

    def test_filter_on_flag(self):
        rna_seq = FeatureSet('Transcripts-ex1.txt')
        rna_flagged = rna_seq.filterByFlag(1)
        self.assertEqual(len(rna_flagged),4,
                         "Wrong number of flagged data lines")

    def test_getTSS(self):
        rna_seq = FeatureSet('Transcripts-ex1.txt')
        rna_plus = rna_seq.filterByStrand('+')
        for rna_data in rna_plus:
            self.assertTrue((rna_data.strand == '+' and
                             rna_data.start == rna_data.getTSS()),
                            "Incorrect TSS on + strand")
        rna_minus = rna_seq.filterByStrand('-')
        for rna_data in rna_minus:
            self.assertTrue((rna_data.strand == '-' and
                             rna_data.end == rna_data.getTSS()),
                            "Incorrect TSS on - strand")

    def test_filter_on_TSS(self):
        rna_seq = FeatureSet('Transcripts-ex1.txt')
        lower,upper = 5000000,10000000
        rna_tss = rna_seq.filterByTSS(upper,lower)
        self.assertEqual(len(rna_tss),3,
                         "Wrong number of transcripts filtered on TSS")
        for rna_data in rna_tss:
            self.assertTrue((rna_data.getTSS() >= lower and
                             rna_data.getTSS() <= upper),
                            "Transcript outside range")

    def test_sort_by_distance(self):
        rna_sort = FeatureSet('Transcripts-ex1.txt')
        position = 4250000
        # Do sort on distance
        # Sort is done in place, so assignment is not required
        # however the sort function should return a reference to
        # the initial object
        result = rna_sort.sortByDistanceFrom(position)
        self.assertEqual(result,rna_sort,
                         "Returned object doesn't match subject")
        # Check that each distance is greater than the previous one
        last_rna_data = None
        for rna_data in rna_sort:
            if not last_rna_data:
                last_rna_data = rna_data
        else:
            self.assertTrue((abs(rna_data.getTSS() - position) >=
                             abs(last_rna_data.getTSS() - position)),
                             "Sort by distance failed")

    def test_sort_by_closest_distance_to_edge(self):
        rna_sort = FeatureSet('Transcripts-ex1.txt')
        position = 4250000
        # Do sort
        # Sort is done in place, so assignment is not required
        # however the sort function should return a reference to
        # the initial object
        result = rna_sort.sortByClosestEdgeTo(position)
        self.assertEqual(result,rna_sort,
                         "Returned object doesn't match subject")
        # Check that the closest distances are in ascending order
        last_rna_data = None
        for rna_data in rna_sort:
            if not last_rna_data:
                last_rna_data = rna_data
            else:
                self.assertTrue((min(abs(rna_data.getTSS() - position),
                                     abs(rna_data.getTES() - position)) >=
                                 min(abs(last_rna_data.getTSS() - position),
                                     abs(last_rna_data.getTES() - position))),
                                "Sort by closest distance to edge failed")

    def test_sort_by_closest_TSS_to_edge(self):
        rna_sort = FeatureSet('Transcripts-ex1.txt')
        position = (16000000,17500000)
        # Do sort
        # Sort is done in place, so assignment is not required
        # however the sort function should return a reference to
        # the initial object
        result = rna_sort.sortByClosestTSSTo(*position)
        self.assertEqual(result,rna_sort,
                         "Returned object doesn't match subject")
        # Check that the closest distances are in ascending order
        last_rna_data = None
        for rna_data in rna_sort:
            if not last_rna_data:
                last_rna_data = rna_data
            else:
                self.assertTrue((min(abs(rna_data.getTSS() - position[0]),
                                     abs(rna_data.getTSS() - position[1])) >=
                                 min(abs(last_rna_data.getTSS() - position[0]),
                                     abs(last_rna_data.getTSS() - position[1]))),
                                "Sort by closest TSS to edge failed")
        

    def test_reading_bad_file_scientific_notation(self):
        self.assertRaises(Exception,FeatureSet,'Transcripts-ex2.txt')

    def test_reading_bad_file_end_lower_than_start(self):
        self.assertRaises(Exception,FeatureSet,'Transcripts-ex2a.txt')

    def test_get_item(self):
        features = FeatureSet('Transcripts-ex1.txt')
        feature = features[2]
        self.assertEqual(feature,Feature('CG32847-RB',
                                         'chr3L',
                                         '15114722',
                                         '15115217',
                                         '+'))

    def test_get_slice(self):
        features = FeatureSet('Transcripts-ex1.txt')
        features_slice = features[1:3]
        self.assertTrue(isinstance(features_slice,FeatureSet))
        self.assertEqual(len(features_slice),2)
        self.assertEqual(features[1],features_slice[0])
        self.assertEqual(features[2],features_slice[1])

    def test__eq__(self):
        # Check equality of FeatureSets
        feature_set1 = FeatureSet()
        feature_set2 = FeatureSet()
        # Empty feature sets
        self.assertEqual(feature_set1,feature_set2)
        # Populate
        feature_set1.addFeature(Feature('CG1000','chr1','1','2','+'))
        feature_set2.addFeature(Feature('CG1000','chr1','1','2','+'))
        self.assertEqual(feature_set1,feature_set2)
        # Add second
        feature_set1.addFeature(Feature('CG2000','chr1','1','2','+'))
        self.assertNotEqual(feature_set1,feature_set2)
        feature_set2.addFeature(Feature('CG2000','chr1','1','2','+'))
        self.assertEqual(feature_set1,feature_set2)
        # Add third
        feature_set1.addFeature(Feature('CG2001','chr2',3,4,'-'))
        feature_set2.addFeature(Feature('CG2002','chr2',3,5,'+'))
        self.assertNotEqual(feature_set1,feature_set2)
