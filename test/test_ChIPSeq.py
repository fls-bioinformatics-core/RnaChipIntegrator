#
#     test_ChIPSeq.py: unit tests for ChIPSeq module
#     Copyright (C) University of Manchester 2011-5 Peter Briggs

from common import *
from rnachipintegrator.ChIPSeq import ChIPSeqData
from rnachipintegrator.RNASeq import RNASeqData
from rnachipintegrator.distances import GetNearestTranscriptToPeak
import unittest

class TestChIPSeqData(unittest.TestCase):

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
        chip_seq = ChIPSeqData('ChIP_peaks-ex1.txt')
        self.assertEqual(len(chip_seq),5,
                         "Wrong number of lines read from ChIP-seq file")

    def test_is_summit_data(self):
        chip_seq = ChIPSeqData('ChIP_peaks-ex1.txt')
        self.assertTrue(chip_seq.isSummit(),
                        "ChIP data are summits")
        chip_seq = ChIPSeqData('ChIP_peaks-ex5.txt')
        self.assertFalse(chip_seq.isSummit(),
                        "ChIP data are not summits")

    def test_filter_on_chromosome(self):
        chip_chr = ChIPSeqData('ChIP_peaks-ex2.txt')
        chromosome = 'chr2L'
        self.assertEqual(len(chip_chr),10,
                         "Wrong number of lines from ChIP-seq file")
        chip_chr = chip_chr.filterByChr(chromosome)
        self.assertEqual(len(chip_chr), 5,
                         "Wrong number of lines from ChIP-seq after chr filter")
        for chip_data in chip_chr:
            self.assertTrue((chip_data.chr == chromosome),
                            "Wrong chromosome name filtered by chr")

    def test_filter_on_peak_position(self):
        chip_seq = ChIPSeqData('ChIP_peaks-ex1.txt')
        lower,upper = 12000000,15000000
        chip_pos = chip_seq.filterByPosition(upper,lower)
        self.assertEqual(len(chip_pos),2,
                         "Wrong number of peaks filtered")
        for chip_data in chip_pos:
            self.assertTrue((chip_data.start >= lower and
                             chip_data.start <= upper),
                            "Peak should have been filtered out")

    def test_sort_by_distance_from(self):
        chip_sort = ChIPSeqData('ChIP_peaks-ex1.txt')
        position = 12000000
        # Do the sorting
        chip_sort.sortByDistanceFrom(position)
        # Check that each distance is greater than the previous one
        last_chip_data = None
        for chip_data in chip_sort:
            if not last_chip_data:
                last_chip_data = chip_data
            else:
                self.assertTrue((abs(chip_data.start - position) >=
                                 abs(last_chip_data.start - position)),
                                "Sort by distance failed")

    def test_ChIP_peak_inside_region(self):
        chip_seq = ChIPSeqData('ChIP_peaks-ex1.txt')
        lower,upper = 4252000,4254000
        self.assertTrue(chip_seq[0].insideRegion(upper,lower),
                        "ChIP peak should be in region")
        self.assertTrue(chip_seq[0].insideRegion(lower,upper),
                        "ChIP peak should be in region (reversed limits)")
        upper,lower = 4252000,4242000
        self.assertFalse(chip_seq[0].insideRegion(upper,lower),
                         "ChIP peak should not be inside region")
        self.assertFalse(chip_seq[0].insideRegion(lower,upper),
                         "ChIP peak should not be inside region (reversed limits")

    def test_handling_quoted_chr(self):
        chip_seq = ChIPSeqData('ChIP_peaks-ex5.txt')
        self.assertEqual(len(chip_seq),5,
                         "Wrong number of lines read from ChIP-seq file")
        chip_chr = chip_seq.filterByChr("chr4")
        self.assertEqual(len(chip_chr),2,
                         "Wrong number of ChIP-seq records filtered")

class TestRNASeqDataWithChIPSeqData(unittest.TestCase):
    def setUp(self):
        # Create input files for tests
        create_test_file('Transcripts-ex1.txt',transcripts_ex1)
        create_test_file('ChIP_peaks-ex1.txt',chip_peaks_ex1)

    def tearDown(self):
        # Remove input files
        delete_test_file('Transcripts-ex1.txt')
        delete_test_file('ChIP_peaks-ex1.txt')

    def test_closest_transcript_to_peak(self):
        rna_seq = RNASeqData('Transcripts-ex1.txt')
        rna_data1 = rna_seq[1]
        rna_data2 = rna_seq[2]
        chip_seq = ChIPSeqData('ChIP_peaks-ex1.txt')
        chip_peak = chip_seq[0]
        nearest = GetNearestTranscriptToPeak(rna_data1,rna_data2,chip_peak)
        # 2nd peak should be closer than first
        self.assertEqual(nearest,rna_data1,
                         "Wrong transcript selected as nearest")
        # test when only one is set
        nearest = GetNearestTranscriptToPeak(None,rna_data2,chip_peak)
        self.assertEqual(nearest,rna_data2,
                         "Wrong transcript selected as nearest")
