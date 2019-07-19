#
#     test_ouptut.py: unit tests for cli module
#     Copyright (C) University of Manchester 2019 Peter Briggs

import unittest
import tempfile
import shutil
import os
from rnachipintegrator.cli import CLI
from rnachipintegrator.cli import OutputFiles
from rnachipintegrator.cli import AnalysisParams
from rnachipintegrator.cli import BatchNamer
from rnachipintegrator.cli import read_feature_file
from rnachipintegrator.cli import read_peak_file

class TestCLI(unittest.TestCase):
    def test_get_version(self):
        cli = CLI("Test",version="v1.0.2")
        self.assertEqual(cli.get_version(),"v1.0.2")

    def test_add_option(self):
        cli = CLI("Test")
        cli.add_option('--attempts',
                       dest='attempts',
                       default=0,
                       type=int,
                       help="Number of attempts")
        options,args = cli.parse_args(list())
        self.assertEqual(options.attempts,0)
        options,args = cli.parse_args(['--attempts=11',])
        self.assertEqual(options.attempts,11)

    def test_add_option_group(self):
        cli = CLI("Test")
        cli.add_option_group("Example options")
        cli.add_option('--attempts',
                       dest='attempts',
                       default=0,
                       type=int,
                       help="Number of attempts",
                       group="Example options")
        options,args = cli.parse_args(list())
        self.assertEqual(options.attempts,0)
        options,args = cli.parse_args(['--attempts=11',])
        self.assertEqual(options.attempts,11)

    def test_add_edge_option(self):
        cli = CLI("Test")
        cli.add_edge_option()
        options,args = cli.parse_args(list())
        self.assertEqual(options.edge,'tss')
        options,args = cli.parse_args(['--edge=tss',])
        self.assertEqual(options.edge,'tss')
        options,args = cli.parse_args(['--edge=both',])
        self.assertEqual(options.edge,'both')

    def test_add_only_de_option(self):
        cli = CLI("Test")
        cli.add_only_de_option()
        options,args = cli.parse_args(list())
        self.assertFalse(options.only_diff_expressed)
        options,args = cli.parse_args(['--only-DE',])
        self.assertTrue(options.only_diff_expressed)

    def test_add_number_option(self):
        cli = CLI("Test")
        cli.add_number_option()
        options,args = cli.parse_args(list())
        self.assertEqual(options.max_closest,None)
        options,args = cli.parse_args(['--number=9',])
        self.assertEqual(options.max_closest,9)

    def test_add_promoter_region_option(self):
        cli = CLI("Test")
        cli.add_promoter_region_option()
        options,args = cli.parse_args(list())
        self.assertEqual(options.promoter_region,"1000,100")
        options,args = cli.parse_args(['--promoter_region=1500,200',])
        self.assertEqual(options.promoter_region,"1500,200")

    def test_name_option(self):
        cli = CLI("Test")
        cli.add_name_option()
        options,args = cli.parse_args(list())
        self.assertEqual(options.name,None)
        options,args = cli.parse_args(['--name=test',])
        self.assertEqual(options.name,"test")

    def test_add_compact_option(self):
        cli = CLI("Test")
        cli.add_compact_option()
        options,args = cli.parse_args(list())
        self.assertFalse(options.compact)
        options,args = cli.parse_args(['--compact',])
        self.assertTrue(options.compact)

    def test_add_summary_option(self):
        cli = CLI("Test")
        cli.add_summary_option()
        options,args = cli.parse_args(list())
        self.assertFalse(options.summary)
        options,args = cli.parse_args(['--summary',])
        self.assertTrue(options.summary)

    def test_add_pad_option(self):
        cli = CLI("Test")
        cli.add_pad_option()
        options,args = cli.parse_args(list())
        self.assertFalse(options.pad_output)
        options,args = cli.parse_args(['--pad',])
        self.assertTrue(options.pad_output)

    def test_add_xlsx_option(self):
        cli = CLI("Test")
        cli.add_xlsx_option()
        options,args = cli.parse_args(list())
        self.assertFalse(options.xlsx_output)
        options,args = cli.parse_args(['--xlsx',])
        self.assertTrue(options.xlsx_output)

    def test_add_analyses_option(self):
        cli = CLI("Test")
        cli.add_analyses_option()
        options,args = cli.parse_args(list())
        self.assertEqual(options.analyses,'all')
        options,args = cli.parse_args(['--analyses=all',])
        self.assertEqual(options.analyses,'all')
        options,args = cli.parse_args(['--analyses=peak_centric',])
        self.assertEqual(options.analyses,'peak_centric')
        options,args = cli.parse_args(['--analyses=gene_centric',])
        self.assertEqual(options.analyses,'gene_centric')

    def test_feature_option(self):
        cli = CLI("Test")
        cli.add_feature_option()
        options,args = cli.parse_args(list())
        self.assertEqual(options.feature_type,'gene')
        options,args = cli.parse_args(['--feature=transcript',])
        self.assertEqual(options.feature_type,"transcript")

    def test_peak_id_option(self):
        cli = CLI("Test")
        cli.add_peak_id_option()
        options,args = cli.parse_args(list())
        self.assertEqual(options.peak_id,None)
        options,args = cli.parse_args(['--peak_id=4',])
        self.assertEqual(options.peak_id,4)

    def test_add_peak_cols_option(self):
        cli = CLI("Test")
        cli.add_peak_cols_option()
        options,args = cli.parse_args(list())
        self.assertEqual(options.peak_cols,None)
        options,args = cli.parse_args(['--peak_cols=3,1,2',])
        self.assertEqual(options.peak_cols,"3,1,2")

class TestOutputFiles(unittest.TestCase):

    def test_properties(self):
        outputs = OutputFiles("test")
        self.assertEqual(outputs.peak_centric_out,
                         "test_peak_centric.txt")
        self.assertEqual(outputs.gene_centric_out,
                         "test_gene_centric.txt")
        self.assertEqual(outputs.peak_centric_summary,
                         "test_peak_centric_summary.txt")
        self.assertEqual(outputs.gene_centric_summary,
                         "test_gene_centric_summary.txt")
        self.assertEqual(outputs.xlsx_out,
                         "test.xlsx")
        self.assertEqual(outputs.files,
                         ("test_peak_centric.txt",
                          "test_gene_centric.txt",
                          "test_peak_centric_summary.txt",
                          "test_gene_centric_summary.txt",
                          "test.xlsx"))

class TestAnalysisParams(unittest.TestCase):

    def test_defaults(self):
        params = AnalysisParams()
        self.assertEqual(params.peaks,None)
        self.assertEqual(params.genes,None)
        self.assertEqual(params.cutoff,None)
        self.assertFalse(params.tss_only)
        self.assertFalse(params.only_differentially_expressed)

    def test_set_values(self):
        params = AnalysisParams(peaks="peaks",
                                genes="genes",
                                cutoff=100000,
                                tss_only=True,
                                only_differentially_expressed=True)
        self.assertEqual(params.peaks,"peaks")
        self.assertEqual(params.genes,"genes")
        self.assertEqual(params.cutoff,100000)
        self.assertTrue(params.tss_only)
        self.assertTrue(params.only_differentially_expressed)

class TestBatchNamer(unittest.TestCase):

    def test_no_batching(self):
        self.assertEqual(BatchNamer(
            "BASE",
            peak_files=('peaks1',),
            gene_files=('genes1',),
            cutoffs=(100000,)).get_name('peaks1',
                                        'genes1',
                                        100000),
                         "BASE")

    def test_multiple_peaks_only(self):
        self.assertEqual(BatchNamer(
            "BASE",
            peak_files=('peaks1','peaks2'),
            gene_files=('genes1',),
            cutoffs=(100000,)).get_name('peaks1',
                                        'genes1',
                                        100000),
                         "BASE_peaks1")

    def test_multiple_genes_only(self):
        self.assertEqual(BatchNamer(
            "BASE",
            peak_files=('peaks1',),
            gene_files=('genes1','genes2'),
            cutoffs=(100000,)).get_name('peaks1',
                                        'genes1',
                                        100000),
                         "BASE_genes1")

    def test_multiple_cutoffs_only(self):
        self.assertEqual(BatchNamer(
            "BASE",
            peak_files=('peaks1',),
            gene_files=('genes1',),
            cutoffs=(100000,150000)).get_name('peaks1',
                                              'genes1',
                                              100000),
                         "BASE_d100000")

    def test_pairs_of_multiples(self):
        self.assertEqual(BatchNamer(
            "BASE",
            peak_files=('peaks1','peaks2',),
            gene_files=('genes1',),
            cutoffs=(100000,150000)).get_name('peaks1',
                                              'genes1',
                                              100000),
                         "BASE_peaks1_d100000")
        self.assertEqual(BatchNamer(
            "BASE",
            peak_files=('peaks1',),
            gene_files=('genes1','genes2'),
            cutoffs=(100000,150000)).get_name('peaks1',
                                              'genes1',
                                              100000),
                         "BASE_genes1_d100000")
        self.assertEqual(BatchNamer(
            "BASE",
            peak_files=('peaks1','peaks2',),
            gene_files=('genes1','genes2'),
            cutoffs=(100000,)).get_name('peaks1',
                                        'genes1',
                                        100000),
                         "BASE_peaks1_genes1")

    def test_all_multiples(self):
        self.assertEqual(BatchNamer(
            "BASE",
            peak_files=('peaks1','peaks2',),
            gene_files=('genes1','genes2'),
            cutoffs=(100000,150000)).get_name('peaks1',
                                              'genes1',
                                              100000),
                         "BASE_peaks1_genes1_d100000")

class TestReadFeatureFile(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if os.path.exists(self.wd):
            shutil.rmtree(self.wd)

    def test_read_feature_file(self):
        feature_data = \
"""
AK030377_A330023F24Rik	chr1	196781953	196826186	-	1
AK080193_A530079E22Rik	chr1	89401837	89403491	-	1
AK082264_C230030N03Rik	chr1	34735043	34781084	+	1
BC006931_AI597479	chr1	43153807	43172843	+	1
"""
        feature_file = os.path.join(self.wd,"features.txt")
        open(feature_file,'w').write(feature_data)
        features = read_feature_file(feature_file)
        self.assertEqual(len(features),4)

    def test_read_feature_file_with_header(self):
        feature_data = \
"""RefSeq_Gene Symbol	chr1	start	stop	strand	diff_exp
AK030377_A330023F24Rik	chr1	196781953	196826186	-	1
AK080193_A530079E22Rik	chr1	89401837	89403491	-	1
AK082264_C230030N03Rik	chr1	34735043	34781084	+	1
BC006931_AI597479	chr1	43153807	43172843	+	1
"""
        feature_file = os.path.join(self.wd,"features.txt")
        open(feature_file,'w').write(feature_data)
        features = read_feature_file(feature_file)
        self.assertEqual(len(features),4)

class TestReadPeakFile(unittest.TestCase):

    def setUp(self):
        self.wd = tempfile.mkdtemp()

    def tearDown(self):
        if os.path.exists(self.wd):
            shutil.rmtree(self.wd)

    def test_read_peak_file(self):
        peak_data = \
"""chr1	9619046	9619167
chr1	9619175	9619382
chr1	10617233	10617437
chr1	13829227	13829564
chr1	17030405	17030545
chr1	17030405	17030545
chr1	24609562	24609623
chr1	24609562	24609623
"""
        peak_file = os.path.join(self.wd,"peaks.txt")
        open(peak_file,'w').write(peak_data)
        peaks = read_peak_file(peak_file)
        self.assertEqual(len(peaks),8)

    def test_read_peak_file_with_header(self):
        peak_data = \
"""peak	s	e
chr1	9619046	9619167
chr1	9619175	9619382
chr1	10617233	10617437
chr1	13829227	13829564
chr1	17030405	17030545
chr1	17030405	17030545
chr1	24609562	24609623
chr1	24609562	24609623
"""
        peak_file = os.path.join(self.wd,"peaks.txt")
        open(peak_file,'w').write(peak_data)
        peaks = read_peak_file(peak_file)
        self.assertEqual(len(peaks),8)

    def test_read_peak_file_cols(self):
        peak_data = \
"""9619046	9619167	chr1
9619175	9619382	chr1
10617233	10617437	chr1
13829227	13829564	chr1
17030405	17030545	chr1
17030405	17030545	chr1
24609562	24609623	chr1
24609562	24609623	chr1
"""
        peak_file = os.path.join(self.wd,"peaks.txt")
        open(peak_file,'w').write(peak_data)
        peaks = read_peak_file(peak_file,peak_cols=(3,1,2))
        self.assertEqual(len(peaks),8)

    def test_read_peak_file_ids(self):
        peak_data = \
"""chr1	9619046	9619167	P001
chr1	9619175	9619382	P002
chr1	10617233	10617437	P003
chr1	13829227	13829564	P004
chr1	17030405	17030545	P005
chr1	17030405	17030545	P006
chr1	24609562	24609623	P007
chr1	24609562	24609623	P008
"""
        peak_file = os.path.join(self.wd,"peaks.txt")
        open(peak_file,'w').write(peak_data)
        peaks = read_peak_file(peak_file,peak_id_col=4)
        self.assertEqual(len(peaks),8)
        for peak in peaks:
            self.assertNotEqual(peak.id,None)
