#
#     test_ouptut.py: unit tests for output module
#     Copyright (C) University of Manchester 2011-5 Peter Briggs

import unittest
from itertools import izip_longest
import rnachipintegrator.output as output
from rnachipintegrator.output import AnalysisReporter,AnalysisReportWriter
from rnachipintegrator.output import describe_fields
from rnachipintegrator.output import update_text
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

    def test_make_header_single_line_custom_feature_type(self):
        ap = AnalysisReporter(output.SINGLE_LINE,
                              fields=('peak.chr',
                                      'peak.start',
                                      'peak.end',
                                      'number_of_results',
                                      'list(feature.id)'),
                              feature_type='gene',
                              max_hits=2)
        self.assertEqual(ap.make_header(),
                         "peak.chr\tpeak.start\tpeak.end\tnumber_of_results\tgene.id_1\tgene.id_2")

    def test_make_header_multi_line_custom_feature_type(self):
        ap = AnalysisReporter(output.MULTI_LINE,
                              fields=('feature.id',
                                      'feature.chr',
                                      'feature.start',
                                      'feature.end',
                                      'peak.chr',
                                      'peak.start',
                                      'peak.end'),
                              feature_type='gene',
                              max_hits=2)
        self.assertEqual(ap.make_header(),
                         "gene.id\tgene.chr\tgene.start\tgene.end\tpeak.chr\tpeak.start\tpeak.end")

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

class TestDescribeFieldsFunction(unittest.TestCase):

    def test_describe_unqualified_fields(self):
        desc = describe_fields(('chr','start','end','id',
                                'strand','TSS','TES',
                                'differentially_expressed'))
        expected_fields = dict()
        for x in ('chr','start','end','id',
                  'strand','TSS','TES',
                  'differentially_expressed'):
            expected_fields[x] = output.FIELDS[x].\
                                 replace("<FEATURE>","feature")
        self.assertEqual(desc[0],('chr',expected_fields['chr']))
        self.assertEqual(desc[1],('start',expected_fields['start']))
        self.assertEqual(desc[2],('end',expected_fields['end']))
        self.assertEqual(desc[3],('id',expected_fields['id']))
        self.assertEqual(desc[4],('strand',expected_fields['strand']))
        self.assertEqual(desc[5],('TSS',expected_fields['TSS']))
        self.assertEqual(desc[6],('TES',expected_fields['TES']))

    def test_describe_peak_fields(self):
        desc = describe_fields(('peak.chr','peak.start','peak.end'))
        self.assertEqual(desc[0],('peak.chr',output.FIELDS['peak.chr']))
        self.assertEqual(desc[1],('peak.start',output.FIELDS['peak.start']))
        self.assertEqual(desc[2],('peak.end',output.FIELDS['peak.end']))

    def test_describe_feature_fields(self):
        desc = describe_fields(('feature.chr','feature.id',
                                'feature.start','feature.end',
                                'feature.strand',
                                'feature.TSS','feature.TES'))
        expected_fields = dict()
        for x in ('feature.chr','feature.id',
                  'feature.start','feature.end',
                  'feature.strand',
                  'feature.TSS','feature.TES'):
            expected_fields[x] = output.FIELDS[x].\
                                 replace("<FEATURE>","feature")
        self.assertEqual(desc[0],('feature.chr',expected_fields['feature.chr']))
        self.assertEqual(desc[1],('feature.id',expected_fields['feature.id']))
        self.assertEqual(desc[2],('feature.start',
                                  expected_fields['feature.start']))
        self.assertEqual(desc[3],('feature.end',expected_fields['feature.end']))
        self.assertEqual(desc[4],('feature.strand',
                                  expected_fields['feature.strand']))
        self.assertEqual(desc[5],('feature.TSS',expected_fields['feature.TSS']))
        self.assertEqual(desc[6],('feature.TES',expected_fields['feature.TES']))

    def test_describe_derived_fields(self):
        desc = describe_fields(('dist_closest',
                                'dist_TSS','dist_TES',
                                'direction',
                                'overlap_feature',
                                'overlap_promoter',
                                'in_the_feature',
                                'order','number_of_results'))
        expected_fields = dict()
        for x in ('dist_closest',
                  'dist_TSS','dist_TES',
                  'direction',
                  'overlap_feature',
                  'overlap_promoter',
                  'in_the_feature',
                  'order','number_of_results'):
            expected_fields[x] = output.FIELDS[x].\
                                 replace("<FEATURE>","feature").\
                                 replace("<SOURCE>","source").\
                                 replace("<TARGET>","target")
        self.assertEqual(desc[0],('dist_closest',
                                  expected_fields['dist_closest']))
        self.assertEqual(desc[1],('dist_TSS',expected_fields['dist_TSS']))
        self.assertEqual(desc[2],('dist_TES',expected_fields['dist_TES']))
        self.assertEqual(desc[3],('direction',expected_fields['direction']))
        self.assertEqual(desc[4],('overlap_feature',
                                  expected_fields['overlap_feature']))
        self.assertEqual(desc[5],('overlap_promoter',
                                  expected_fields['overlap_promoter']))
        self.assertEqual(desc[6],('in_the_feature',
                                  expected_fields['in_the_feature']))
        self.assertEqual(desc[7],('order',expected_fields['order']))
        self.assertEqual(desc[8],('number_of_results',
                                  expected_fields['number_of_results']))

    def test_describe_listed_fields_for_peaks(self):
        desc = describe_fields(('id','list(chr,start,end,dist_TSS)'))
        expected_fields = dict()
        for x in ('id','chr','start','end','dist_TSS'):
            expected_fields[x] = output.FIELDS[x].\
                                 replace("<FEATURE>","feature")
        self.assertEqual(desc[0],('id',expected_fields['id']))
        self.assertEqual(desc[1],(('For each hit:',)))
        self.assertEqual(desc[2],('chr_#',expected_fields['chr']))
        self.assertEqual(desc[3],('start_#',expected_fields['start']))
        self.assertEqual(desc[4],('end_#',expected_fields['end']))
        self.assertEqual(desc[5],('dist_TSS_#',expected_fields['dist_TSS']))

    def test_describe_listed_fields_for_features(self):
        desc = describe_fields(('chr','start','list(id,TSS,strand)'))
        expected_fields = dict()
        for x in ('chr','start','id','TSS','strand'):
            expected_fields[x] = output.FIELDS[x].\
                                 replace("<FEATURE>","feature")
        self.assertEqual(desc[0],('chr',expected_fields['chr']))
        self.assertEqual(desc[1],('start',expected_fields['start']))
        self.assertEqual(desc[2],(('For each hit:',)))
        self.assertEqual(desc[3],('id_#',expected_fields['id']))
        self.assertEqual(desc[4],('TSS_#',expected_fields['TSS']))
        self.assertEqual(desc[5],('strand_#',expected_fields['strand']))

    def test_update_placeholders_in_fields(self):
        desc = describe_fields(('chr','start','list(id,TSS,strand,direction)'),
                               feature="gene",
                               source="peak",
                               target="gene")
        expected_fields = dict()
        for x in ('chr','start','id','TSS','strand','direction'):
            expected_fields[x] = output.FIELDS[x].\
                                 replace("<FEATURE>","gene").\
                                 replace("<SOURCE>","peak").\
                                 replace("<TARGET>","gene")
        self.assertEqual(desc[0],('chr',expected_fields['chr']))
        self.assertEqual(desc[1],('start',expected_fields['start']))
        self.assertEqual(desc[2],(('For each hit:',)))
        self.assertEqual(desc[3],('id_#',expected_fields['id']))
        self.assertEqual(desc[4],('TSS_#',expected_fields['TSS']))
        self.assertEqual(desc[5],('strand_#',expected_fields['strand']))
        self.assertEqual(desc[6],('direction_#',
                                  expected_fields['direction']))

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

class TestUpdateTextFunction(unittest.TestCase):

    def test_with_no_placeholders(self):
        self.assertEqual(update_text("Hello you!"),
                         "Hello you!")

    def test_with_placeholder(self):
        self.assertEqual(update_text("Hello <NAME>!",NAME="Joe"),
                         "Hello Joe!")

    def test_with_placeholder_no_value_supplied(self):
        self.assertEqual(update_text("Hello <NAME>!"),
                         "Hello <NAME>!")

    def test_with_placeholder_value_is_None(self):
        self.assertEqual(update_text("Hello <NAME>!",NAME=None),
                         "Hello <NAME>!")

    def test_with_multiple_placeholders(self):
        self.assertEqual(update_text(
            "Get <TRANSPORT_MODE> from <HERE> to <THERE>",
            TRANSPORT_MODE="bus",
            HERE="home",
            THERE="town"),
                         "Get bus from home to town")
