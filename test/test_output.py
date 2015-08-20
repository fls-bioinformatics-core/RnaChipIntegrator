#
#     test_ouptut.py: unit tests for output module
#     Copyright (C) University of Manchester 2011-5 Peter Briggs

import unittest
import rnachipintegrator.output as output
from rnachipintegrator.output import AnalysisReporter
from rnachipintegrator.Peaks import Peak,PeakSet
from rnachipintegrator.Features import Feature,FeatureSet

class TestAnalysisReporter(unittest.TestCase):

    def test_make_header_single_line(self):
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=('peak.chr',
                                      'peak.start',
                                      'peak.end',
                                      'number_of_results',
                                      'list(feature.id)'))
        self.assertEqual(ap.make_header(0),
                         "peak.chr\tpeak.start\tpeak.end\tnumber_of_results")
        self.assertEqual(ap.make_header(2),
                         "peak.chr\tpeak.start\tpeak.end\tnumber_of_results\tfeature.id_1\tfeature.id_2")

    def test_make_header_multi_line(self):
        ap = AnalysisReporter(output.MULTI_LINE,
                              fields=('feature.id',
                                      'feature.chr',
                                      'feature.start',
                                      'feature.end',
                                      'peak.chr',
                                      'peak.start',
                                      'peak.end'))
        self.assertEqual(ap.make_header(0),
                         "feature.id\tfeature.chr\tfeature.start\tfeature.end\tpeak.chr\tpeak.start\tpeak.end")
        self.assertEqual(ap.make_header(2),
                         "feature.id\tfeature.chr\tfeature.start\tfeature.end\tpeak.chr\tpeak.start\tpeak.end")

    def test_report_nearest_features_single_line(self):
        # Set up some test data
        peak = Peak('chr2L',66811,66812)
        features = FeatureSet(
            features_list=(
                Feature('CG31973','chr2L',25402,59243,'-'),
                Feature('CG2674-RE','chr2L',106903,114433,'+'),
                Feature('CG2674-RC','chr2L',107926,114433,'+')))
        # Expected lines
        expected = (
            "chr2L\t66811\t66812\t3\t"
            "CG31973\t7568\t7568\t"
            "CG2674-RE\t40091\t40091\t"
            "CG2674-RC\t41114\t41114",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=('peak.chr',
                                      'peak.start',
                                      'peak.end',
                                      'number_of_results',
                                      'list(feature.id,'
                                      'dist_closest,dist_TSS)'))
        # Check that output matches
        for line,expected_line in zip(ap.report_nearest_features(peak,
                                                                 features),
                                      expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_features_multi_line(self):
        # Set up some test data
        peak = Peak('chr2L',66811,66812)
        features = FeatureSet(
            features_list=(
                Feature('CG31973','chr2L',25402,59243,'-'),
                Feature('CG2674-RE','chr2L',106903,114433,'+'),
                Feature('CG2674-RC','chr2L',107926,114433,'+')))
        # Expected lines
        expected = (
            "chr2L\t66811\t66812\t1 of 3\tCG31973\t7568\t7568",
            "chr2L\t66811\t66812\t2 of 3\tCG2674-RE\t40091\t40091",
            "chr2L\t66811\t66812\t3 of 3\tCG2674-RC\t41114\t41114",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.MULTI_LINE,
                              fields=('peak.chr',
                                      'peak.start',
                                      'peak.end',
                                      'order',
                                      'feature.id',
                                      'dist_closest',
                                      'dist_TSS'))
        # Check that output matches
        for line,expected_line in zip(ap.report_nearest_features(peak,
                                                                 features),
                                      expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_peaks_single_line(self):
        # Set up some test data
        feature = Feature('CG31973','chr2L','25402','59243','-')
        peaks = PeakSet(
            peaks_list=(
                Peak('chr2L','66711','66911'),
                Peak('chr2L','249077','249277'),
                Peak('chr2L','605850','606050')))
        # Expected lines
        expected = (
            "CG31973\t3\t"
            "chr2L\t66711\t66911\t7468\t7468\t"
            "chr2L\t249077\t249277\t189834\t189834\t"
            "chr2L\t605850\t606050\t546607\t546607",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=('feature.id',
                                      'number_of_results',
                                      'list(peak.chr,peak.start,peak.end,'
                                      'dist_closest,dist_TSS)'))
        # Check that output matches
        for line,expected_line in zip(ap.report_nearest_peaks(feature,
                                                              peaks),
                                      expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_peaks_multi_line(self):
        # Set up some test data
        feature = Feature('CG31973','chr2L','25402','59243','-')
        peaks = PeakSet(
            peaks_list=(
                Peak('chr2L','66711','66911'),
                Peak('chr2L','249077','249277'),
                Peak('chr2L','605850','606050')))
        # Expected lines
        expected = (
            "CG31973\t1 of 3\tchr2L\t66711\t66911\t7468\t7468",
            "CG31973\t2 of 3\tchr2L\t249077\t249277\t189834\t189834",
            "CG31973\t3 of 3\tchr2L\t605850\t606050\t546607\t546607",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.MULTI_LINE,
                              fields=('feature.id',
                                      'order',
                                      'peak.chr',
                                      'peak.start',
                                      'peak.end',
                                      'dist_closest',
                                      'dist_TSS'))
        # Check that output matches
        for line,expected_line in zip(ap.report_nearest_peaks(feature,
                                                              peaks),
                                      expected):
            self.assertEqual(line,expected_line)
