#
#     test_Peaks.py: unit tests for Peaks module
#     Copyright (C) University of Manchester 2011-5 Peter Briggs

from common import *
from rnachipintegrator.Peaks import PeakSet
from rnachipintegrator.Peaks import Peak
from rnachipintegrator.Features import FeatureSet
from rnachipintegrator.distances import GetNearestTranscriptToPeak
import unittest

class TestPeakSet(unittest.TestCase):

    def setUp(self):
        # Create input files for tests
        create_test_file('ChIP_peaks-ex1.txt',chip_peaks_ex1)
        create_test_file('ChIP_peaks-ex2.txt',chip_peaks_ex2)
        create_test_file('ChIP_peaks-ex5.txt',chip_peaks_ex5)

    def tearDown(self):
        # Remove input files
        delete_test_file('ChIP_peaks-ex1.txt')
        delete_test_file('ChIP_peaks-ex2.txt')
        delete_test_file('ChIP_peaks-ex5.txt')

    def test_reading_in_ChIPseq_data(self):
        peaks = PeakSet('ChIP_peaks-ex1.txt')
        self.assertEqual(len(peaks),5,
                         "Wrong number of lines read from ChIP-seq file")

    def test_is_summit_data(self):
        peaks = PeakSet('ChIP_peaks-ex1.txt')
        self.assertTrue(peaks.isSummit(),"ChIP data are summits")
        peaks = PeakSet('ChIP_peaks-ex5.txt')
        self.assertFalse(peaks.isSummit(),"ChIP data are not summits")

    def test_filter_on_chromosome(self):
        peaks_chr = PeakSet('ChIP_peaks-ex2.txt')
        chromosome = 'chr2L'
        self.assertEqual(len(peaks_chr),10,
                         "Wrong number of lines from ChIP-seq file")
        peaks_chr = peaks_chr.filterByChr(chromosome)
        self.assertEqual(len(peaks_chr), 5,
                         "Wrong number of lines from ChIP-seq after chr filter")
        for peak in peaks_chr:
            self.assertTrue((peak.chrom == chromosome),
                            "Wrong chromosome name filtered by chr")

    def test_filter_on_peak_position(self):
        peaks = PeakSet('ChIP_peaks-ex1.txt')
        lower,upper = 12000000,15000000
        peaks_pos = peaks.filterByPosition(upper,lower)
        self.assertEqual(len(peaks_pos),2,
                         "Wrong number of peaks filtered")
        for peak in peaks_pos:
            print str(peak)
            self.assertTrue((peak.start >= lower and
                             peak.start <= upper),
                            "Peak should have been filtered out")

    def test_sort_by_distance_from(self):
        peaks_sort = PeakSet('ChIP_peaks-ex1.txt')
        position = 12000000
        # Do the sorting
        peaks_sort.sortByDistanceFrom(position)
        # Check that each distance is greater than the previous one
        last_peak = None
        for peak in peaks_sort:
            if not last_peak:
                last_peak = peak
            else:
                self.assertTrue((abs(peak.start - position) >=
                                 abs(last_peak.start - position)),
                                "Sort by distance failed")

    def test_ChIP_peak_inside_region(self):
        peaks = PeakSet('ChIP_peaks-ex1.txt')
        lower,upper = 4252000,4254000
        self.assertTrue(peaks[0].insideRegion(upper,lower),
                        "ChIP peak should be in region")
        self.assertTrue(peaks[0].insideRegion(lower,upper),
                        "ChIP peak should be in region (reversed limits)")
        upper,lower = 4252000,4242000
        self.assertFalse(peaks[0].insideRegion(upper,lower),
                         "ChIP peak should not be inside region")
        self.assertFalse(peaks[0].insideRegion(lower,upper),
                         "ChIP peak should not be inside region (reversed limits")

    def test_handling_quoted_chr(self):
        peaks = PeakSet('ChIP_peaks-ex5.txt')
        self.assertEqual(len(peaks),5,
                         "Wrong number of lines read from ChIP-seq file")
        peaks_chr = peaks.filterByChr("chr4")
        self.assertEqual(len(peaks_chr),2,
                         "Wrong number of ChIP-seq records filtered")

    def test_get_item(self):
        peaks = PeakSet('ChIP_peaks-ex1.txt')
        peak = peaks[2]
        self.assertEqual(peak,Peak('chr3L','12139192','12139193'))

    def test_get_slice(self):
        peaks = PeakSet('ChIP_peaks-ex1.txt')
        peaks_slice = peaks[1:3]
        self.assertTrue(isinstance(peaks_slice,PeakSet))
        self.assertEqual(len(peaks_slice),2)
        self.assertEqual(peaks[1],peaks_slice[0])
        self.assertEqual(peaks[2],peaks_slice[1])

class TestPeak(unittest.TestCase):

    def test__eq__(self):
        self.assertTrue(Peak('chr2L','66811','66812') == 
                        Peak('chr2L','66811','66812'))
        self.assertFalse(Peak('chr2L','66811','66812') == 
                         Peak('chr2L','249177','605951'))

class TestFeatureSetWithChIPSeqData(unittest.TestCase):
    def setUp(self):
        # Create input files for tests
        create_test_file('Transcripts-ex1.txt',transcripts_ex1)
        create_test_file('ChIP_peaks-ex1.txt',chip_peaks_ex1)

    def tearDown(self):
        # Remove input files
        delete_test_file('Transcripts-ex1.txt')
        delete_test_file('ChIP_peaks-ex1.txt')

    def test_closest_transcript_to_peak(self):
        features = FeatureSet('Transcripts-ex1.txt')
        feature1 = features[1]
        feature2 = features[2]
        peaks = PeakSet('ChIP_peaks-ex1.txt')
        peak = peaks[0]
        nearest = GetNearestTranscriptToPeak(feature1,feature2,peak)
        # 2nd feature should be closer than first
        self.assertEqual(nearest,feature1,
                         "Wrong transcript selected as nearest")
        # test when only one is set
        nearest = GetNearestTranscriptToPeak(None,feature2,peak)
        self.assertEqual(nearest,feature2,
                         "Wrong transcript selected as nearest")
