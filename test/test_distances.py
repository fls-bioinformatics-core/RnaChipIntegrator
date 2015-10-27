#
#     test_analysis.py: unit tests for analysis module
#     Copyright (C) University of Manchester 2011-5 Peter Briggs

from common import *
from rnachipintegrator.Peaks import Peak
from rnachipintegrator.Features import Feature
from rnachipintegrator.distances import regions_overlap
from rnachipintegrator.distances import closestDistanceToRegion
from rnachipintegrator.distances import edge_distances
from rnachipintegrator.distances import tss_distances
from rnachipintegrator.distances import direction
import unittest

########################################################################
#
# TestRegionsOverlap
#
#########################################################################

class TestRegionsOverlap(unittest.TestCase):
    
    def test_regions_dont_overlap(self):
        region1 = (1000,2000)
        region2 = (3000,4000)
        self.assertFalse(regions_overlap(region1,region2),
                         "Regions should not overlap")
        self.assertFalse(regions_overlap(region2,region1),
                         "Regions should not overlap (reversed)")

    def test_regions_partially_overlap(self):
        region1 = (1000,2000)
        region2 = (1500,2500)
        self.assertTrue(regions_overlap(region1,region2),
                        "Regions should overlap")
        self.assertTrue(regions_overlap(region2,region1),
                        "Regions should overlap (reversed)")

    def test_one_region_inside_the_other(self):
        region1 = (1000,2000)
        region2 = (1250,1750)
        self.assertTrue(regions_overlap(region1,region2),
                        "Regions should overlap")
        self.assertTrue(regions_overlap(region2,region1),
                        "Regions should overlap (reversed)")

    def test_regions_overlap_limit_orders_dont_matter(self):
        # Check that limits of regions can be in any order
        self.assertTrue(regions_overlap((0,10),(5,15)))
        self.assertTrue(regions_overlap((10,0),(5,15)))
        self.assertTrue(regions_overlap((0,10),(15,5)))
        self.assertTrue(regions_overlap((10,0),(15,5)))

########################################################################
#
# TestClosestDistanceToRegion
#
#########################################################################

class TestClosestDistanceToRegion(unittest.TestCase):

    def test_get_distance_to_point(self):
        # Single position
        reference = 1240000
        position = 1200000
        # Single reference position
        self.assertEqual(closestDistanceToRegion(reference,position),
                         abs(reference-position),
                         "Incorrect distance #1")
        position = 1300000
        self.assertEqual(closestDistanceToRegion(reference,position),
                         abs(reference-position),
                         "Incorrect distance #2")

    def test_get_closest_distance_outside_region(self):
        # Outside reference region (2 positions)
        reference = 1240000
        position1 = 1100000
        position2 = 1200000
        self.assertEqual(closestDistanceToRegion(reference,position1,position2),
                         abs(reference-position2),
                         "Incorrect closest distance (outside region #1)")
        position1 = 1300000
        position2 = 1400000
        self.assertEqual(closestDistanceToRegion(reference,position1,position2),
                         abs(reference-position1),
                         "Incorrect closest distance (outside region #2)")

    def test_get_closest_distance_inside_region(self):
        # Inside reference region
        reference = 1240000
        position1 = 1200000
        position2 = 1255000
        self.assertEqual(closestDistanceToRegion(reference,position1,position2),
                         abs(reference-position2),
                         "Incorrect closest distance (inside region #1)")
        self.assertEqual(closestDistanceToRegion(reference,position1,position2,
                                                 zero_inside_region=True),
                         0,
                         "Incorrect closest distance (inside region #2)")

########################################################################
#
# TestEdgeDistancesFunction
#
#########################################################################

class TestEdgeDistancesFunction(unittest.TestCase):
    def test_distances_peak_before_feature(self):
        self.assertEqual(edge_distances(Peak('chr1','100','200'),
                                        Feature('NM1','chr1','250','400','+')),
                         (50,200))
    def test_distances_feature_before_peak(self):
        self.assertEqual(edge_distances(Peak('chr1','250','400'),
                                   Feature('NM2','chr1','100','200','+')),
                         (50,150))
    def test_distances_peak_overlaps_feature_start(self):
        self.assertEqual(edge_distances(Peak('chr1','100','250'),
                                   Feature('NM3','chr1','200','400','+')),
                         (0,100))
    def test_distances_feature_overlaps_peak_start(self):
        self.assertEqual(edge_distances(Peak('chr1','250','350'),
                                   Feature('NM4','chr1','300','400','+')),
                         (0,50))
    def test_distances_peak_contains_feature(self):
        self.assertEqual(edge_distances(Peak('chr1','250','450'),
                                   Feature('NM5','chr1','300','400','+')),
                         (0,0))
    def test_distances_feature_contains_peak(self):
        self.assertEqual(edge_distances(Peak('chr1','250','350'),
                                   Feature('NM6','chr1','100','400','+')),
                         (0,0))

########################################################################
#
# TestTSSDistancesFunction
#
#########################################################################

class TestTSSDistancesFunction(unittest.TestCase):
    def test_tss_distances_peak_before_TSS(self):
        self.assertEqual(tss_distances(Peak('chr1','100','200'),
                                       Feature('NM1','chr1','250','400','+')),
                         (50,150))
        self.assertEqual(tss_distances(Peak('chr1','100','200'),
                                       Feature('NM1','chr1','250','400','-')),
                         (200,300))
    def test_tss_distances_TSS_before_peak(self):
        self.assertEqual(tss_distances(Peak('chr1','250','400'),
                                       Feature('NM2','chr1','100','200','+')),
                         (150,300))
        self.assertEqual(tss_distances(Peak('chr1','250','400'),
                                       Feature('NM2','chr1','100','200','-')),
                         (50,200))
    def test_tss_distances_peak_contains_TSS(self):
        self.assertEqual(tss_distances(Peak('chr1','100','250'),
                                       Feature('NM3','chr1','200','400','+')),
                         (0,75))
        self.assertEqual(tss_distances(Peak('chr1','250','350'),
                                       Feature('NM4','chr1','200','300','-')),
                         (0,50))

########################################################################
#
# TestDirectionFunction
#
#########################################################################

from rnachipintegrator.distances import UPSTREAM,DOWNSTREAM,OVERLAP

class TestDirectionFunction(unittest.TestCase):
    # Positive strand, peak relative to feature
    def test_peak_upstream_from_feature(self):
        peak = Peak('chr1','100','200')
        feature = Feature('NM1','chr1','250','400','+')
        self.assertEqual(direction(feature,peak),UPSTREAM)
    def test_peak_downstream_from_feature(self):
        peak = Peak('chr1','450','550')
        feature = Feature('NM1','chr1','250','400','+')
        self.assertEqual(direction(feature,peak),DOWNSTREAM)
    def test_peak_upstream_from_feature_partial_overlap(self):
        peak = Peak('chr1','100','300')
        feature = Feature('NM1','chr1','250','400','+')
        self.assertEqual(direction(feature,peak),UPSTREAM)
    def test_peak_downstream_from_feature_partial_overlap(self):
        peak = Peak('chr1','350','550')
        feature = Feature('NM1','chr1','250','400','+')
        self.assertEqual(direction(feature,peak),DOWNSTREAM)
    def test_peak_full_overlap_feature(self):
        peak = Peak('chr1','200','450')
        feature = Feature('NM1','chr1','250','400','+')
        self.assertEqual(direction(feature,peak),OVERLAP)
    # Negative strand, peak relative to feature
    def test_peak_upstream_from_feature_neg_strand(self):
        peak = Peak('chr1','100','200')
        feature = Feature('NM1','chr1','250','400','-')
        self.assertEqual(direction(feature,peak),DOWNSTREAM)
    def test_peak_downstream_from_feature_neg_strand(self):
        peak = Peak('chr1','450','550')
        feature = Feature('NM1','chr1','250','400','-')
        self.assertEqual(direction(feature,peak),UPSTREAM)
    def test_peak_upstream_from_feature_neg_strand_partial_overlap(self):
        peak = Peak('chr1','100','300')
        feature = Feature('NM1','chr1','250','400','-')
        self.assertEqual(direction(feature,peak),DOWNSTREAM)
    def test_peak_downstream_from_feature_neg_strand_partial_overlap(self):
        peak = Peak('chr1','350','550')
        feature = Feature('NM1','chr1','250','400','-')
        self.assertEqual(direction(feature,peak),UPSTREAM)
    def test_peak_full_overlap_feature_neg_strand(self):
        peak = Peak('chr1','200','450')
        feature = Feature('NM1','chr1','250','400','-')
        self.assertEqual(direction(feature,peak),OVERLAP)
    # Positive strand, feature relative to peak
    def test_feature_upstream_from_peak(self):
        peak = Peak('chr1','450','550')
        feature = Feature('NM1','chr1','250','400','+')
        self.assertEqual(direction(peak,feature),UPSTREAM)
    def test_feature_downstream_from_peak(self):
        peak = Peak('chr1','100','200')
        feature = Feature('NM1','chr1','250','400','+')
        self.assertEqual(direction(peak,feature),DOWNSTREAM)
    def test_feature_upstream_from_peak_partial_overlap(self):
        peak = Peak('chr1','350','550')
        feature = Feature('NM1','chr1','250','400','+')
        self.assertEqual(direction(peak,feature),UPSTREAM)
    def test_feature_downstream_from_peak_partial_overlap(self):
        peak = Peak('chr1','100','300')
        feature = Feature('NM1','chr1','250','400','+')
        self.assertEqual(direction(peak,feature),DOWNSTREAM)
    def test_feature_full_overlap_peak(self):
        peak = Peak('chr1','200','450')
        feature = Feature('NM1','chr1','250','400','+')
        self.assertEqual(direction(peak,feature),OVERLAP)
    # Negative strand, feature relative to peak
    def test_feature_upstream_from_peak_neg_strand(self):
        peak = Peak('chr1','100','200')
        feature = Feature('NM1','chr1','250','400','-')
        self.assertEqual(direction(peak,feature),UPSTREAM)
    def test_feature_downstream_from_peak_neg_strand(self):
        peak = Peak('chr1','450','550')
        feature = Feature('NM1','chr1','250','400','-')
        self.assertEqual(direction(peak,feature),DOWNSTREAM)
    def test_feature_upstream_from_peak_neg_strand_partial_overlap(self):
        peak = Peak('chr1','100','300')
        feature = Feature('NM1','chr1','250','400','-')
        self.assertEqual(direction(peak,feature),UPSTREAM)
    def test_feature_downstream_from_peak_neg_strand_partial_overlap(self):
        peak = Peak('chr1','350','550')
        feature = Feature('NM1','chr1','250','400','-')
        self.assertEqual(direction(peak,feature),DOWNSTREAM)
    def test_feature_full_overlap_peak_neg_strand(self):
        peak = Peak('chr1','200','450')
        feature = Feature('NM1','chr1','250','400','-')
        self.assertEqual(direction(peak,feature),OVERLAP)
