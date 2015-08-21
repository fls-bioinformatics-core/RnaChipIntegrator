#
#     test_ouptut.py: unit tests for output module
#     Copyright (C) University of Manchester 2011-5 Peter Briggs

import unittest
from itertools import izip_longest
import rnachipintegrator.output as output
from rnachipintegrator.output import AnalysisReporter,AnalysisReportWriter
from rnachipintegrator.Peaks import Peak,PeakSet
from rnachipintegrator.Features import Feature,FeatureSet

class TestAnalysisReporterHeader(unittest.TestCase):

    def test_make_header_single_line(self):
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=('peak.chr',
                                      'peak.start',
                                      'peak.end',
                                      'number_of_results',
                                      'list(feature.id)'))
        self.assertEqual(ap.make_header(),
                         "peak.chr\tpeak.start\tpeak.end\tnumber_of_results")
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=('peak.chr',
                                      'peak.start',
                                      'peak.end',
                                      'number_of_results',
                                      'list(feature.id)'),
                              max_hits=2)
        self.assertEqual(ap.make_header(),
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
        self.assertEqual(ap.make_header(),
                         "feature.id\tfeature.chr\tfeature.start\tfeature.end\tpeak.chr\tpeak.start\tpeak.end")
        ap = AnalysisReporter(output.MULTI_LINE,
                              fields=('feature.id',
                                      'feature.chr',
                                      'feature.start',
                                      'feature.end',
                                      'peak.chr',
                                      'peak.start',
                                      'peak.end'),
                              max_hits=2)
        self.assertEqual(ap.make_header(),
                         "feature.id\tfeature.chr\tfeature.start\tfeature.end\tpeak.chr\tpeak.start\tpeak.end")

class TestAnalysisReporterNearestFeatures(unittest.TestCase):

    def setUp(self):
        # Set up some test data
        self.peak = Peak('chr2L',66811,66812)
        self.features = FeatureSet(
            features_list=(
                Feature('CG31973','chr2L',25402,59243,'-'),
                Feature('CG2674-RE','chr2L',106903,114433,'+'),
                Feature('CG2674-RC','chr2L',107926,114433,'+')))
        self.single_line_fields = ('peak.chr',
                                   'peak.start',
                                   'peak.end',
                                   'number_of_results',
                                   'list(feature.id,'
                                   'dist_closest,dist_TSS)')
        self.multi_line_fields = ('peak.chr',
                                  'peak.start',
                                  'peak.end',
                                  'order',
                                  'feature.id',
                                  'dist_closest',
                                  'dist_TSS')

    def test_report_nearest_features_single_line(self):
        # Expected lines
        expected = (
            "chr2L\t66811\t66812\t3\t"
            "CG31973\t7568\t7568\t"
            "CG2674-RE\t40091\t40091\t"
            "CG2674-RC\t41114\t41114",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=self.single_line_fields)
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_features(self.peak,
                                           self.features),
                expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_features_single_line_max_hits(self):
        # Expected lines
        expected = (
            "chr2L\t66811\t66812\t2\t"
            "CG31973\t7568\t7568\t"
            "CG2674-RE\t40091\t40091",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=self.single_line_fields,
                              max_hits=2)
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_features(self.peak,
                                           self.features),
                expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_features_single_line_pad(self):
        # Expected lines
        expected = (
            "chr2L\t66811\t66812\t3\t"
            "CG31973\t7568\t7568\t"
            "CG2674-RE\t40091\t40091\t"
            "CG2674-RC\t41114\t41114\t"
            ".\t.\t.",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=self.single_line_fields,
                              max_hits=4,
                              pad=True,
                              null_placeholder='.')
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_features(self.peak,
                                           self.features),
                expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_features_multi_line(self):
        # Expected lines
        expected = (
            "chr2L\t66811\t66812\t1 of 3\tCG31973\t7568\t7568",
            "chr2L\t66811\t66812\t2 of 3\tCG2674-RE\t40091\t40091",
            "chr2L\t66811\t66812\t3 of 3\tCG2674-RC\t41114\t41114",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.MULTI_LINE,
                              fields=self.multi_line_fields)
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_features(self.peak,
                                           self.features),
                expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_features_multi_line_max_hits(self):
        # Expected lines
        expected = (
            "chr2L\t66811\t66812\t1 of 2\tCG31973\t7568\t7568",
            "chr2L\t66811\t66812\t2 of 2\tCG2674-RE\t40091\t40091",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.MULTI_LINE,
                              fields=self.multi_line_fields,
                              max_hits=2)
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_features(self.peak,
                                           self.features),
                expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_features_multi_line_pad(self):
        # Expected lines
        expected = (
            "chr2L\t66811\t66812\t1 of 3\tCG31973\t7568\t7568",
            "chr2L\t66811\t66812\t2 of 3\tCG2674-RE\t40091\t40091",
            "chr2L\t66811\t66812\t3 of 3\tCG2674-RC\t41114\t41114",
            "chr2L\t66811\t66812\t.\t.\t.\t.",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.MULTI_LINE,
                              fields=self.multi_line_fields,
                              max_hits=4,
                              pad=True,
                              null_placeholder='.')
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_features(self.peak,
                                           self.features),
                expected):
            self.assertEqual(line,expected_line)

class TestAnalysisReporterNearestPeaks(unittest.TestCase):

    def setUp(self):
        # Set up some test data
        self.feature = Feature('CG31973','chr2L','25402','59243','-')
        self.peaks = PeakSet(
            peaks_list=(
                Peak('chr2L','66711','66911'),
                Peak('chr2L','249077','249277'),
                Peak('chr2L','605850','606050')))
        self.single_line_fields = ('feature.id',
                                   'number_of_results',
                                   'list(peak.chr,peak.start,peak.end,'
                                   'dist_closest,dist_TSS)')
        self.multi_line_fields = ('feature.id',
                                  'order',
                                  'peak.chr',
                                  'peak.start',
                                  'peak.end',
                                  'dist_closest',
                                  'dist_TSS')

    def test_report_nearest_peaks_single_line(self):
        # Expected lines
        expected = (
            "CG31973\t3\t"
            "chr2L\t66711\t66911\t7468\t7468\t"
            "chr2L\t249077\t249277\t189834\t189834\t"
            "chr2L\t605850\t606050\t546607\t546607",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=self.single_line_fields)
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_peaks(self.feature,
                                        self.peaks),
                expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_peaks_single_line_max_hits(self):
        # Expected lines
        expected = (
            "CG31973\t2\t"
            "chr2L\t66711\t66911\t7468\t7468\t"
            "chr2L\t249077\t249277\t189834\t189834",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=self.single_line_fields,
                              max_hits=2,
                              null_placeholder='.')
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_peaks(self.feature,
                                        self.peaks),
                expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_peaks_single_line_pad(self):
        # Expected lines
        expected = (
            "CG31973\t3\t"
            "chr2L\t66711\t66911\t7468\t7468\t"
            "chr2L\t249077\t249277\t189834\t189834\t"
            "chr2L\t605850\t606050\t546607\t546607\t"
            ".\t.\t.\t.\t.",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=self.single_line_fields,
                              max_hits=4,
                              pad=True)
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_peaks(self.feature,
                                        self.peaks),
                expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_peaks_multi_line(self):
        # Expected lines
        expected = (
            "CG31973\t1 of 3\tchr2L\t66711\t66911\t7468\t7468",
            "CG31973\t2 of 3\tchr2L\t249077\t249277\t189834\t189834",
            "CG31973\t3 of 3\tchr2L\t605850\t606050\t546607\t546607",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.MULTI_LINE,
                              fields=self.multi_line_fields)
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_peaks(self.feature,
                                        self.peaks),
                expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_peaks_multi_line_max_hits(self):
        # Expected lines
        expected = (
            "CG31973\t1 of 2\tchr2L\t66711\t66911\t7468\t7468",
            "CG31973\t2 of 2\tchr2L\t249077\t249277\t189834\t189834",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.MULTI_LINE,
                              fields=self.multi_line_fields,
                              max_hits=2)
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_peaks(self.feature,
                                        self.peaks),
                expected):
            self.assertEqual(line,expected_line)

    def test_report_nearest_peaks_multi_line_pad(self):
        # Expected lines
        expected = (
            "CG31973\t1 of 3\tchr2L\t66711\t66911\t7468\t7468",
            "CG31973\t2 of 3\tchr2L\t249077\t249277\t189834\t189834",
            "CG31973\t3 of 3\tchr2L\t605850\t606050\t546607\t546607",
            "CG31973\t.\t.\t.\t.\t.\t.",
        )
        # Set up to report some stuff
        ap = AnalysisReporter(output.MULTI_LINE,
                              fields=self.multi_line_fields,
                              max_hits=4,
                              pad=True,
                              null_placeholder='.')
        # Check that output matches
        for line,expected_line in izip_longest(
                ap.report_nearest_peaks(self.feature,
                                        self.peaks),
                expected):
            self.assertEqual(line,expected_line)

import tempfile
class TestAnalysisReportWriter(unittest.TestCase):

    def test_write_features(self):
        # Set up some test data
        peak = Peak('chr2L',66811,66812)
        features = FeatureSet(
            features_list=(
                Feature('CG31973','chr2L',25402,59243,'-'),
                Feature('CG2674-RE','chr2L',106903,114433,'+'),
                Feature('CG2674-RC','chr2L',107926,114433,'+')))
        # Temp output file
        fp,outfile = tempfile.mkstemp()
        # Write peaks to file
        ap = AnalysisReportWriter(output.MULTI_LINE,
                                  fields=('peak.chr',
                                          'peak.start',
                                          'peak.end',
                                          'order',
                                          'feature.id',
                                          'dist_closest',
                                          'dist_TSS'),
                                  outfile=outfile)
        ap.write_nearest_features(peak,features)
        ap.close()
        # Expected and actual output
        expected_output = \
            "#peak.chr\tpeak.start\tpeak.end\torder\tfeature.id\tdist_closest\tdist_TSS\n" \
            "chr2L\t66811\t66812\t1 of 3\tCG31973\t7568\t7568\n" \
            "chr2L\t66811\t66812\t2 of 3\tCG2674-RE\t40091\t40091\n" \
            "chr2L\t66811\t66812\t3 of 3\tCG2674-RC\t41114\t41114\n"
        actual_output = open(outfile,'r').read()
        # Check that output matches
        self.assertEqual(expected_output,actual_output)

    def test_write_features_summary(self):
        # Set up some test data
        peak = Peak('chr2L',66811,66812)
        features = FeatureSet(
            features_list=(
                Feature('CG31973','chr2L',25402,59243,'-'),
                Feature('CG2674-RE','chr2L',106903,114433,'+'),
                Feature('CG2674-RC','chr2L',107926,114433,'+')))
        # Temp output file
        fp,summary = tempfile.mkstemp()
        # Write peaks to file
        ap = AnalysisReportWriter(output.MULTI_LINE,
                                  fields=('peak.chr',
                                          'peak.start',
                                          'peak.end',
                                          'order',
                                          'feature.id',
                                          'dist_closest',
                                          'dist_TSS'),
                                  summary=summary)
        ap.write_nearest_features(peak,features)
        ap.close()
        # Expected and actual output
        expected_output = \
            "#peak.chr\tpeak.start\tpeak.end\torder\tfeature.id\tdist_closest\tdist_TSS\n" \
            "chr2L\t66811\t66812\t1 of 3\tCG31973\t7568\t7568\n"
        actual_output = open(summary,'r').read()
        # Check that output matches
        self.assertEqual(expected_output,actual_output)

    def test_write_peaks(self):
        # Set up some test data
        feature = Feature('CG31973','chr2L','25402','59243','-')
        peaks = PeakSet(
            peaks_list=(
                Peak('chr2L','66711','66911'),
                Peak('chr2L','249077','249277'),
                Peak('chr2L','605850','606050')))
        # Temp output file
        fp,outfile = tempfile.mkstemp()
        # Write peaks to file
        ap = AnalysisReportWriter(output.MULTI_LINE,
                                  fields=('feature.id',
                                          'order',
                                          'peak.chr',
                                          'peak.start',
                                          'peak.end',
                                          'dist_closest',
                                          'dist_TSS'),
                                  outfile=outfile)
        ap.write_nearest_peaks(feature,peaks)
        ap.close()
        # Expected and actual output
        expected_output = \
            "#feature.id\torder\tpeak.chr\tpeak.start\tpeak.end\tdist_closest\tdist_TSS\n" \
            "CG31973\t1 of 3\tchr2L\t66711\t66911\t7468\t7468\n" \
            "CG31973\t2 of 3\tchr2L\t249077\t249277\t189834\t189834\n" \
            "CG31973\t3 of 3\tchr2L\t605850\t606050\t546607\t546607\n"
        actual_output = open(outfile,'r').read()
        # Check that output matches
        self.assertEqual(expected_output,actual_output)

    def test_write_peaks_summary(self):
        # Set up some test data
        feature = Feature('CG31973','chr2L','25402','59243','-')
        peaks = PeakSet(
            peaks_list=(
                Peak('chr2L','66711','66911'),
                Peak('chr2L','249077','249277'),
                Peak('chr2L','605850','606050')))
        # Temp output files
        fp,summary= tempfile.mkstemp()
        # Write peaks to file
        ap = AnalysisReportWriter(output.MULTI_LINE,
                                  fields=('feature.id',
                                          'order',
                                          'peak.chr',
                                          'peak.start',
                                          'peak.end',
                                          'dist_closest',
                                          'dist_TSS'),
                                  summary=summary)
        ap.write_nearest_peaks(feature,peaks)
        ap.close()
        # Expected and actual output
        expected_output = \
            "#feature.id\torder\tpeak.chr\tpeak.start\tpeak.end\tdist_closest\tdist_TSS\n" \
            "CG31973\t1 of 3\tchr2L\t66711\t66911\t7468\t7468\n"
        actual_output = open(summary,'r').read()
        # Check that output matches
        self.assertEqual(expected_output,actual_output)
