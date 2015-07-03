#
#     test_analysis.py: unit tests for analysis module
#     Copyright (C) University of Manchester 2011-5 Peter Briggs

from common import *
from rnachipintegrator.analysis import AnalyseNearestTSSToSummits
from rnachipintegrator.analysis import AnalyseNearestTranscriptsToPeakEdges
from rnachipintegrator.analysis import AnalyseNearestPeaksToTranscripts
from rnachipintegrator.ChIPSeq import ChIPSeqData
from rnachipintegrator.RNASeq import RNASeqData
import unittest

########################################################################
#
# TestAnalyseNearestTSSToSummits
#
#########################################################################

class TestAnalyseNearestTSSToSummits(unittest.TestCase):

    def setUp(self):
        # Create input files for tests
        create_test_file('Transcripts-ex3.txt',transcripts_ex3)
        create_test_file('ChIP_peaks-ex3.txt',chip_peaks_ex3)

    def tearDown(self):
        # Remove input files
        delete_test_file('Transcripts-ex3.txt')
        delete_test_file('ChIP_peaks-ex3.txt')

    def test_AnalyseNearestTSSToSummits(self):
        # Initialise data
        rna_seq = RNASeqData('Transcripts-ex3.txt')
        chip_seq = ChIPSeqData('ChIP_peaks-ex3.txt')
        max_distance = 130000
        # Run the analysis
        results = AnalyseNearestTSSToSummits(chip_seq,
                                             rna_seq,
                                             max_distance)
        # Verify the results are as expected
        self.assertEqual(len(results),
                         nearest_transcript_to_peak_ex3.count('\n'),
                         "Didn't get expected number of results")
        for i in range(len(results)):
            # Get the expected line from the reference output
            expected = nearest_transcript_to_peak_ex3.split('\n')[i+1]
            # Build a line from the actual results
            items = []
            for field in ('chr','start','geneID','nearest','TSS',
                          'distance_to_TSS','distance_to_TES',
                          'strand','in_the_gene','transcripts_inbetween',
                          'transcript_ids_inbetween'):
                items.append(str(results[i][field]))
            actual = '\t'.join(items)
            self.assertEqual(actual,expected,
                             "Result differs for line %s\nRef:%s\nAct:%s" % \
                                 (i,expected,actual))

########################################################################
#
# TestAnalyseNearestTranscriptsToPeakEdges
#
#########################################################################

class TestAnalyseNearestTranscriptsToPeakEdges(unittest.TestCase):

    def setUp(self):
        # Create input files for tests
        create_test_file('Transcripts-ex3.txt',transcripts_ex3)
        create_test_file('ChIP_peaks_binding_region-ex3.txt',
                         chip_peaks_binding_region_ex3)

    def tearDown(self):
        # Remove input files
        delete_test_file('Transcripts-ex3.txt')
        delete_test_file('ChIP_peaks_binding_region-ex3.txt')

    def test_AnalyseNearestTranscriptsToPeakEdges(self):
        # Initialise data
        rna_seq = RNASeqData('Transcripts-ex3.txt')
        chip_seq = ChIPSeqData('ChIP_peaks_binding_region-ex3.txt')
        promoter_region = (10000,2500)
        max_closest=4
        # Remove the flag from the gene data
        self.assertTrue(rna_seq.isFlagged(),"initial gene data should be flagged")
        for data in rna_seq:
            data.flag = None
        self.assertFalse(rna_seq.isFlagged(),"failed to remove flag on gene data")
        # Run the analysis
        results = AnalyseNearestTranscriptsToPeakEdges(chip_seq,
                                                       rna_seq,
                                                       promoter_region,
                                                       max_closest)
        # Verify the results are as expected
        self.assertEqual(len(results),
                         nearest_transcript_to_peak_edge_ex3.count('\n'),
                         "Didn't get expected number of results")
        for i in range(len(results)):
            # Get the expected line from the reference output
            expected = nearest_transcript_to_peak_edge_ex3.split('\n')[i+1]
            # Build a line from the actual results
            items = []
            for field in ('chr','start','end','geneID','strand','TSS','TES',
                          'dist_closest_edge','dist_TSS','dist_TES',
                          'overlap_transcript','overlap_promoter'):
                items.append(str(results[i][field]))
            actual = '\t'.join(items)
            self.assertEqual(actual,expected,
                             "Result differs for line %s\nRef:%s\nAct:%s" % \
                                 (i+1,expected,actual))

    def test_AnalyseNearestTranscriptsToPeakEdges_with_diff_expression(self):
        # Initialise data
        rna_seq = RNASeqData('Transcripts-ex3.txt')
        chip_seq = ChIPSeqData('ChIP_peaks_binding_region-ex3.txt')
        promoter_region = (10000,2500)
        max_closest=4
        # Check the differential expression flag on the gene data
        self.assertTrue(rna_seq.isFlagged(),"gene data should be flagged")
        # Run the analysis
        results = AnalyseNearestTranscriptsToPeakEdges(chip_seq,
                                                       rna_seq,
                                                       promoter_region,
                                                       max_closest)
        # Verify the results are as expected
        self.assertEqual(len(results),
                         nearest_transcript_to_peak_edge_diff_expression_ex3.count('\n'),
                         "Didn't get expected number of results")
        for i in range(len(results)):
            # Get the expected line from the reference output
            expected = nearest_transcript_to_peak_edge_diff_expression_ex3.split('\n')[i+1]
            # Build a line from the actual results
            items = []
            for field in ('chr','start','end','geneID','strand','TSS','TES',
                          'dist_closest_edge','dist_TSS','dist_TES',
                          'overlap_transcript','overlap_promoter'):
                items.append(str(results[i][field]))
            actual = '\t'.join(items)
            self.assertEqual(actual,expected,
                             "Result differs for line %s\nRef:%s\nAct:%s" % \
                                 (i+1,expected,actual))

    def test_AnalyseNearestTranscriptTSSToPeakEdges(self):
        # Initialise data
        rna_seq = RNASeqData('Transcripts-ex3.txt')
        chip_seq = ChIPSeqData('ChIP_peaks_binding_region-ex3.txt')
        promoter_region = (10000,2500)
        max_closest=4
        # Remove the flag from the gene data
        self.assertTrue(rna_seq.isFlagged(),"initial gene data should be flagged")
        for data in rna_seq:
            data.flag = None
        self.assertFalse(rna_seq.isFlagged(),"failed to remove flag on gene data")
        # Run the analysis
        results = AnalyseNearestTranscriptsToPeakEdges(chip_seq,
                                                       rna_seq,
                                                       promoter_region,
                                                       max_closest,
                                                       TSS_only=True)
        # Verify the results are as expected
        self.assertEqual(len(results),
                         nearest_transcript_tss_to_peak_edge_ex3.count('\n'),
                         "Didn't get expected number of results")
        for i in range(len(results)):
            # Get the expected line from the reference output
            expected = nearest_transcript_tss_to_peak_edge_ex3.split('\n')[i+1]
            # Build a line from the actual results
            items = []
            for field in ('chr','start','end','geneID','strand','TSS','TES',
                          'dist_closest_edge','dist_TSS','dist_TES',
                          'overlap_transcript','overlap_promoter'):
                items.append(str(results[i][field]))
            actual = '\t'.join(items)
            self.assertEqual(actual,expected,
                             "Result differs for line %s\nRef:%s\nAct:%s" % \
                                 (i+1,expected,actual))

    def test_AnalyseNearestTranscriptTSSToPeakEdges_with_diff_expression(self):
        # Initialise data
        rna_seq = RNASeqData('Transcripts-ex3.txt')
        chip_seq = ChIPSeqData('ChIP_peaks_binding_region-ex3.txt')
        promoter_region = (10000,2500)
        max_closest=4
        # Check the differential expression flag on the gene data
        self.assertTrue(rna_seq.isFlagged(),"gene data should be flagged")
        # Run the analysis
        results = AnalyseNearestTranscriptsToPeakEdges(chip_seq,
                                                       rna_seq,
                                                       promoter_region,
                                                       max_closest,
                                                       TSS_only=True)
        # Verify the results are as expected
        self.assertEqual(len(results),
                         nearest_transcript_tss_to_peak_edge_diff_expression_ex3.count('\n'),
                         "Didn't get expected number of results")
        for i in range(len(results)):
            # Get the expected line from the reference output
            expected = nearest_transcript_tss_to_peak_edge_diff_expression_ex3.split('\n')[i+1]
            # Build a line from the actual results
            items = []
            for field in ('chr','start','end','geneID','strand','TSS','TES',
                          'dist_closest_edge','dist_TSS','dist_TES',
                          'overlap_transcript','overlap_promoter'):
                items.append(str(results[i][field]))
            actual = '\t'.join(items)
            self.assertEqual(actual,expected,
                             "Result differs for line %s\nRef:%s\nAct:%s" % \
                                 (i+1,expected,actual))

########################################################################
#
# TestAnalyseNearestPeaksToTranscripts
#
#########################################################################

class TestAnalyseNearestPeaksToTranscripts(unittest.TestCase):

    def setUp(self):
        # Create input files for tests
        create_test_file('Transcripts-ex4.txt',transcripts_ex4)
        create_test_file('ChIP_peaks-ex4.txt',chip_peaks_ex4)

    def tearDown(self):
        # Remove input files
        delete_test_file('Transcripts-ex4.txt')
        delete_test_file('ChIP_peaks-ex4.txt')

    def test_AnalyseNearestTranscriptsToPeaks(self):
        # Initialise data
        rna_seq = RNASeqData('Transcripts-ex4.txt')
        chip_seq = ChIPSeqData('ChIP_peaks-ex4.txt')
        window_width = 20000
        # Run the analysis
        results = AnalyseNearestPeaksToTranscripts(rna_seq,
                                                   chip_seq,
                                                   window_width)
        # Verify the results are as expected
        self.assertEqual(len(results),
                         nearest_peak_to_transcript_ex4.count('\n'),
                         "Didn't get expected number of results")
        for i in range(len(results)):
            # Get the expected line from the reference output
            expected = nearest_peak_to_transcript_ex4.split('\n')[i+1]
            # Build a line from the actual results
            items = []
            for field in ('geneID','chr_RNA','start','end','strand',
                          'differentially_expressed','number_of_peaks',
                          'chr_ChIP_1','summit_1','distance_1',
                          'chr_ChIP_2','summit_2','distance_2',
                          'chr_ChIP_3','summit_3','distance_3'):
                try:
                    items.append(str(results[i][field]))
                except KeyError:
                    items.append('---')
            actual = '\t'.join(items)
            self.assertEqual(actual,expected,
                             "Result differs for line %s\nRef:'%s'\nAct:'%s'" \
                                 % (i,expected,actual))
