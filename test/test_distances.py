#
#     test_analysis.py: unit tests for analysis module
#     Copyright (C) University of Manchester 2011-5 Peter Briggs

from common import *
from rnachipintegrator.distances import regions_overlap
from rnachipintegrator.distances import closestDistanceToRegion
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
