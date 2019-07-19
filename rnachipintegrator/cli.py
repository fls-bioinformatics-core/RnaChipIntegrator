#!/usr/bin/env python
#
#     RnaChipIntegrator.py: analyse genomic features (genes) with peak data
#     Copyright (C) University of Manchester 2011-19 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
#     This code is free software; you can redistribute it and/or modify it
#     under the terms of the Artistic License 2.0 (see the file LICENSE
#     included with the distribution).
#
########################################################################
#
# RnaChipIntegrator.py
#
#########################################################################

"""
RnaChipIntegrator.py

Analyse genomic features (genes) with peak data.
"""

#######################################################################
# Imports
#######################################################################

import sys
import os
import argparse
from .Features import FeatureSet
from .Peaks import PeakSet
import analysis
import output
import xls_output
import logging
from multiprocessing import Pool

from . import get_version
__version__ = get_version()

logging.getLogger().setLevel(logging.WARNING)
logging.basicConfig(format='%(levelname)s: %(message)s')

#######################################################################
# Data
#######################################################################

# Default values
class _DEFAULTS(object):
    PROMOTER_REGION = (1000,100)
    CUTOFF = 1000000
    MAX_CLOSEST = None
    PAD_OUTPUT = False
    FEATURE_TYPE = 'gene'

# Front matter
_PROGRAM_INFO = """
Find nearest peaks to genes (and vice versa)

University of Manchester
Faculty of Biology Medicine and Health
Bioinformatics Core Facility
Authors: Peter Briggs, Ian Donaldson and Leo Zeef

If you use this program in your published work then please cite:

   Briggs PJ, Donaldson IJ, Zeef LAH. RnaChipIntegrator
   (version %s). Available at:
   https://github.com/fls-bioinformatics-core/RnaChipIntegrator
""" % __version__

#######################################################################
# Classes
#######################################################################

class CLI(object):
    """
    Utility class for building command line parser

    The CLI class provides a simple way to construct
    command line parsers for RnaChipIntegrator programs.
    It wraps an OptionParser instance with convenience
    methods for adding options and groups, alongside a
    set of pre-implemented options.

    Example usage:

    >>> p = CLI("Example utility")
    >>> p.add_cutoff_option()
    >>> p.add_option('--force',action='store_true',dest='force')
    >>> p.parse_args()
    >>> print(options.cutoff)
    >>> print(options.force)

    """
    def __init__(self,usage,version=None,description=None):
        """
        Create a new CLI instance

        Arguments:
          usage (str): usage string for the utility
          version (str): optional, set the version
            for the utility which is displayed by the
            '--version' option. If not supplied then
            defaults to the program name plus the
            version number of the package.
          description (str): optional, text
            describing the utility which is displayed
            by the '--help' option
        """
        self._prog = os.path.basename(sys.argv[0])
        if version is None:
            self._version = "%(prog)s "+__version__
        else:
            self._version = version
        self._version = self._version % { 'prog': self._prog }
        self.parser = argparse.ArgumentParser(
            prog=self._prog,
            usage=usage,
            version=self._version,
            description=description)
        self.option_groups = dict()

    def get_version(self):
        """
        Return the version string from the parser
        """
        return self._version

    def parse_args(self,args=None):
        """
        Parse an argument list

        Arguments:
          args (list): optional, list of options
            to parse (otherwise defaults to
            command line supplied to the utility)

        Returns:
          Tuple: tuple (options,arguments).
        """
        if args is None:
            args = sys.argv[1:]
        return self.parser.parse_known_args(args)

    def error(self,*args):
        """
        Raise an error from within the parser

        Wrapper for the 'error' method of the underlying
        parser object

        Arguments:
          args (list): arguments to be passed to the
            'error' method.
        """
        return self.parser.error(*args)

    def add_option_group(self,name):
        """
        Add an option group to the parser

        Arguments:
          name (str): name/title text for the group
        """
        group = self.parser.add_argument_group(name)
        self.option_groups[name] = group

    def add_option(self,*args,**kws):
        """
        Add an option to the parser

        Passes arguments and keywords to the
        optparse.add_option method.

        If keywords contain the 'group' keyword
        then the option is added to the named
        group (must have been created by a call
        to 'add_option_group').

        Arguments:
          args (list): any argument accepted by
            optparse.add_option
          kws (mapping): any keyword accepted by
            optparse.add_option
        """
        if 'group' in kws:
            group = kws['group']
            del(kws['group'])
        else:
            group = None
        if group is not None:
            p = self.option_groups[group]
        else:
            p = self.parser
        p.add_argument(*args,**kws)

    def add_edge_option(self,group=None):
        """
        Add --edge option to the parser

        Acceptable values are 'tss' (default) or 'both'.
        The value is accessed via the 'edge' property of
        the parser options.
        """
        self.add_option('--edge',
                        action='store',dest="edge",
                        choices=('tss','both'),
                        default='tss',
                        help="Gene edges to consider when calculating "
                        "distances between genes and peaks, either: "
                        "'tss' (default: only use gene TSS) or 'both' "
                        "(use whichever of TSS or TES gives shortest "
                        "distance)",
                        group=group)

    def add_only_de_option(self,group=None):
        """
        Add --only-DE option to the parser

        Default value is False. The value is accessed via
        the 'only_diff_expressed' property of the parser
        options.
        """
        self.add_option('--only-DE',
                        action='store_true',dest='only_diff_expressed',
                        default=False,
                        help="Only use genes flagged as differentially "
                        "expressed in analyses (input gene data must "
                        "include DE flags)",
                        group=group)

    def add_number_option(self,group=None):
        """
        Add --number option to the parser

        Default value is None, otherwise it is an integer.
        The value is accessed via the 'max_closest' property
        of the parser options.
        """
        self.add_option('--number',
                        action='store',dest='max_closest',
                        type=int,default=_DEFAULTS.MAX_CLOSEST,
                        help="Filter results after applying --cutoff[s] "
                        "to report only the nearest MAX_CLOSEST number "
                        "of pairs for each gene/peak from the analyses "
                        "(default is to report all results)",
                        group=group)

    def add_promoter_region_option(self,group=None):
        """
        Add --promoter_region option to the parser

        Default value is "1000,100". The value is accessed
        via the 'promoter_region' property of the parser
        options.
        """
        self.add_option('--promoter_region',
                        action="store",dest="promoter_region",
                        default="%d,%d" % _DEFAULTS.PROMOTER_REGION,
                        help="Define promoter region with respect to "
                        "gene TSS in the form UPSTREAM,DOWNSTREAM "
                        "(default -%d to %dbp of TSS)" %
                        _DEFAULTS.PROMOTER_REGION,
                        group=group)

    def add_name_option(self,group=None):
        """
        Add --name to the parser

        Default value is None. The value is accessed
        via the 'name' property of the parser options.
        """
        self.add_option('--name',
                        action='store',dest='name',
                        default=None,
                        help="Set basename for output files",
                        group=group)

    def add_compact_option(self,group=None):
        """
        Add --compact option to the parser

        Default value is False. The value is accessed via
        the 'compact' property of the parser options.
        """
        self.add_option('--compact',
                        action='store_true',dest='compact',
                        default=False,
                        help="Output all hits for each peak or gene on a "
                        "single line (cannot be used with --summary)",
                        group=group)

    def add_summary_option(self,group=None):
        """
        Add --summary option to the parser

        Default value is False. The value is accessed via
        the 'summary' property of the parser options.
        """
        self.add_option('--summary',
                        action='store_true',dest='summary',
                        default=False,
                        help="Output 'summary' for each analysis, "
                        "consisting of only the top hit for each peak "
                        "or gene (cannot be used with --compact)",
                        group=group)

    def add_pad_option(self,group=None):
        """
        Add --pad option to the parser

        Default value is False. The value is accessed via
        the 'pad_output' property of the parser options.
        """
        self.add_option('--pad',
                        action="store_true",dest="pad_output",
                        help="Where less than MAX_CLOSEST hits are found, "
                        "pad output with blanks to ensure that MAX_CLOSEST "
                        "hits are still reported (nb --pad is implied for "
                        "--compact)",
                        group=group)

    def add_xlsx_option(self,group=None):
        """
        Add --xlsx option to the parser

        Default value is False. The value is accessed via
        the 'xlsx_output' property of the parser options.
        """
        self.add_option('--xlsx',
                        action="store_true",dest="xlsx_output",
                        help="Output XLSX spreadsheet with results",
                        group=group)

    def add_analyses_option(self,group=None):
        """
        Add --analyses option to the parser

        Acceptable values are 'all' (default), 'peak_centric',
        or 'gene_centric'. The value is accessed via the
        'analyses' property of the parser options.
        """
        self.add_option('--analyses',
                        action='store',dest="analyses",
                        choices=('all','gene_centric','peak_centric',),
                        default="all",
                        help="Select which analyses to run: can be one "
                        "of 'all' (default, runs all available "
                        "analyses), 'peak_centric' or 'gene_centric'",
                        group=group)

    def add_feature_option(self,group=None):
        """
        Add --feature option to the parser


        Default value is 'gene'. The value is accessed
        via the 'feature_type' property of the parser
        options.
        """
        self.add_option('--feature',
                        action="store",dest="feature_type",
                        default=_DEFAULTS.FEATURE_TYPE,
                        help="Rename '%s' to FEATURE_TYPE in output "
                        "(e.g. 'transcript' etc)"
                        % _DEFAULTS.FEATURE_TYPE,
                        group=group)

    def add_peak_id_option(self,group=None):
        """
        Add --peak_id option to the parser

        Default value is None. The value is accessed
        via the 'peak_id' property of the parser options.
        """
        self.add_option('--peak_id',
                        action="store",dest="peak_id",
                        type=int,
                        help="Column to use as an ID for each peak "
                        "from the input peak file (first column is "
                        "column 1). If specified then IDs will be "
                        "transferred to the output files when peaks "
                        "are reported",
                        group=group)

    def add_peak_cols_option(self,group=None):
        """
        Add --peak_cols option to the parser

        Default value is None. The value is accessed
        via the 'peak_cols' property of the parser options.
        """
        self.add_option('--peak_cols',
                        action="store",dest="peak_cols",
                        help="List of 3 column indices (e.g. '1,4,5') "
                        "indicating columns to use for chromosome, "
                        "start and end from the input peak file (if not "
                        "first three columns)",
                        group=group)

class OutputFiles(object):
    """
    Utility class for handling output file naming/clean up

    Given a basename, provides methods to fetch
    the name for various output files.

    Also provides a method to automatically clean
    up (i.e. remove) files which would otherwise
    be overwritten by analysis.
    """
    # Set up output file names
    def __init__(self,basename):
        """
        Create new OutputFiles instance

        Arguments:
          basename (str): basename to generate
            file names using
        """
        self._basename = str(basename)

    @property
    def peak_centric_out(self):
        """
        Output file for peak-centric analyses
        """
        return self._basename+"_peak_centric.txt"

    @property
    def gene_centric_out(self):
        """
        Output file for gene-centric analyses
        """
        return self._basename+"_gene_centric.txt"

    @property
    def peak_centric_summary(self):
        """
        Output file for summary of peak-centric analyses
        """
        return self._basename+"_peak_centric_summary.txt"

    @property
    def gene_centric_summary(self):
        """
        Output file for summary of gene-centric analyses
        """
        return self._basename+"_gene_centric_summary.txt"

    @property
    def xlsx_out(self):
        """
        Output XLSX file name
        """
        return self._basename+".xlsx"

    @property
    def files(self):
        """
        List of generated file names
        """
        return (self.peak_centric_out,
                self.gene_centric_out,
                self.peak_centric_summary,
                self.gene_centric_summary,
                self.xlsx_out)

    def remove_files(self):
        """
        Remove output files which already exist
        """
        print("Removing pre-existing output files for '%s'" %
              self._basename)
        for f in self.files:
            if os.path.isfile(f):
                print("\tRemoving %s" % f)
                os.remove(f)
        print("")

class AnalysisParams(object):
    """
    Class for passing around analysis parameters

    Class that enables parameters for analyses to
    be stored (at instantiation) and accessed via
    object properties, e.g.

    >>> params = AnalysisParams(peaks=peakset,genes=geneset)
    >>> params.peaks

    The following parameters are available:

    genes: FeatureSet instance
    peaks: PeakSet instance
    cutoff: cutoff distance
    tss_only: whether to use TSS or both TSS and TES
    only_differentially_expressed: whether to restrict
      genes to those flagged as differenitally
      expressed
    """
    def __init__(self,**kws):
        """
        Create AnalysisParams instance

        Arguments:
          kws (mapping): parameters and their
           associated values
        """
        self._params = dict(
            genes=None,
            peaks=None,
            cutoff=None,
            tss_only=False,
            only_differentially_expressed=False)
        for key in kws:
            if key not in self._params:
                raise AttributeError("'%s' object has no attribute "
                                     "'%s'" % (self.__name__,key))
            self._params[key] = kws[key]
    @property
    def genes(self):
        return self._params['genes']
    @property
    def peaks(self):
        return self._params['peaks']
    @property
    def cutoff(self):
        return self._params['cutoff']
    @property
    def tss_only(self):
        return self._params['tss_only']
    @property
    def only_differentially_expressed(self):
        return self._params['only_differentially_expressed']

class BatchNamer(object):
    """
    Class for generating output names for batch mode

    Names will be of the form:

    {BASENAME}[_{PEAK_FILE}][_{GENE_FILE}][_d{CUTOFF}]

    where elements in square braces are only included if
    there are multiple values across all analyses.

    """
    def __init__(self,basename,peak_files,gene_files,cutoffs):
        """
        Create new Namer instance

        Arguments:
          basename (str): basename for the generated names
          peak_files (list): list of peak files across all
            analyses
          gene_files (list): list of gene files across all
            analyses
          cutoffs (list): list of cutoffs across all analyses
        """
        self._basename = basename
        self._multi_peaks = len(peak_files) > 1
        self._multi_genes = len(gene_files) > 1
        self._multi_cutoffs = len(cutoffs) > 1
        self._delim = '_'

    def get_name(self,peak_file,gene_file,cutoff):
        """
        Return a name based on the supplied arguments

        Arguments:
          peak_file (str): name/path for peak file
          gene_file (str): name/path for gene file
          cutoff (int): cutoff distance used
        """
        peak_file = os.path.splitext(os.path.basename(peak_file))[0]
        gene_file = os.path.splitext(os.path.basename(gene_file))[0]
        return "%s%s%s%s" % (self._basename,
                             ("%s%s" % (self._delim,peak_file)
                              if self._multi_peaks else ""),
                             ("%s%s" % (self._delim,gene_file)
                              if self._multi_genes else ""),
                             ("%sd%s" % (self._delim,cutoff)
                              if self._multi_cutoffs else ""))

#######################################################################
# Functions
#######################################################################

def read_feature_file(feature_file,feature_type='gene'):
    """
    Read in feature data from file

    Arguments:
      feature_file (str): path to file to read data from
      feature_type (str): optional, name for feature type
        (default: 'gene')

    Returns:
      FeatureSet: features read in from the file.
    """

    print("Reading in data from '%s'" % feature_file)
    try:
        features = FeatureSet(feature_file)
    except Exception as ex:
        logging.fatal("Failed to read in %s data: %s" % (feature_type,
                                                         ex))
        print("Please fix errors in input file before running again")
        sys.exit(1)
    if not len(features):
        logging.fatal("No %s data read in" % feature_type)
        sys.exit(1)
    print("%d %s records read in" % (len(features),
                                     feature_type))
    if features.isFlagged():
        print("\tData include differential expression flag")
        print("\t%d records flagged as differentially expressed" %
              len(features.filterByFlag(1)))
    print("")
    return features

def read_peak_file(peak_file,peak_cols=None,peak_id_col=None):
    """
    Read in peak data from file

    Arguments:
      peak_file (str): path to file to read data from
      peak_cols (tuple): optional, tuple of column numbers
        (starting from 1) for the chromosome, start and end
        columns
      peak_id_col (int): optional, number of the column
        to read IDs from (starting from 1)

    Returns:
      PeakSet: peaks read in from the file.
    """
    # Read in peak data
    print("Reading in data from '%s'" % peak_file)
    print("Using columns %s from peaks file as chrom, start, end" %
          (peak_cols,))
    if peak_id_col is not None:
        peak_id_col = int(peak_id_col)
        print("Using column %s from peaks file as peak ID" % peak_id_col)
    try:
        peaks = PeakSet(peak_file,columns=peak_cols,
                        id_column=peak_id_col)
    except Exception as ex:
        logging.fatal("Failed to read peak data (%s)" % ex)
        print("Please fix errors in input file before running again")
        sys.exit(1)
    if not len(peaks):
        logging.error("No peak data read in")
        sys.exit(1)
    print("%d peak records read in (%s)" % (len(peaks),
                                            'summits' if peaks.isSummit()
                                            else 'regions'))
    print("")
    return peaks

def find_nearest_features(params):
    """
    Wrapper for 'find_nearest_features' using AnalysisParams

    Arguments:
      params (AnalysisParams): params controlling how the
        nearest genes are determined

    Yields:
      Tuple: tuple of (peak,nearest_genes,params)
    """
    for peak,genes in analysis.find_nearest_features(
            params.peaks,
            params.genes,
            tss_only=params.tss_only,
            distance=params.cutoff,
            only_differentially_expressed=\
            params.only_differentially_expressed):
        yield (peak,genes,params)

def find_nearest_peaks(params):
    """
    Wrapper for 'find_nearest_peaks' using AnalysisParams

    Arguments:
      params (AnalysisParams): params controlling how the
        nearest peaks are determined

    Yields:
      Tuple: tuple of (gene,nearest_peaks,params)
    """
    for gene,peaks in analysis.find_nearest_peaks(
            params.genes,
            params.peaks,
            tss_only=params.tss_only,
            distance=params.cutoff,
            only_differentially_expressed=\
            params.only_differentially_expressed):
        yield (gene,peaks,params)

def find_nearest_features_as_list(params):
    """
    Wrapper to fetch results from 'find_nearest_features' as a list
    """
    return list(find_nearest_features(params))

def find_nearest_peaks_as_list(params):
    """
    Wrapper to fetch results from 'find_nearest_peaks' as a list
    """
    return list(find_nearest_peaks(params))

#######################################################################
# Main program
#######################################################################

def main(args=None):
    """
    Implements the 'RnaChipIntegrator' CLI
    """

    # Command line arguments
    if args is None:
        args = sys.argv[1:]

    p = CLI(usage="%(prog)s [options] [GENES PEAKS]",
            description=
            "Analyse GENES (any set of genes or genomic "
            "features) against PEAKS (a set of regions) "
            "and report nearest genes to each peak (and "
            "vice versa)")

    p.add_option_group("Analysis options")
    p.add_option('--cutoff',action='store',
                 dest='max_distance',default=None,
                 help="Maximum distance allowed between peaks "
                 "and genes before being omitted from the "
                 "analyses (default %dbp; set to zero for no "
                 "cutoff, use --cutoffs instead to specify "
                 "multiple distances)" % _DEFAULTS.CUTOFF,
                 group="Analysis options")
    p.add_edge_option(group="Analysis options")
    p.add_only_de_option(group="Analysis options")

    p.add_option_group("Reporting options")
    p.add_number_option(group="Reporting options")
    p.add_promoter_region_option(group="Reporting options")

    p.add_option_group("Output options")
    p.add_name_option(group="Output options")
    p.add_compact_option(group="Output options")
    p.add_summary_option(group="Output options")
    p.add_pad_option(group="Output options")
    p.add_xlsx_option(group="Output options")

    p.add_option_group("Batch options")
    p.add_option('--cutoffs',action='store',dest='cutoffs',
                 default=None,
                 help="Comma-separated list of one or more "
                 "maximum distances allowed between peaks "
                 "and genes (bp). An analysis will be "
                 "performed for each GENES-PEAKS pair at "
                 "each cutoff distance (default %dbp; set "
                 "to zero for no cutoff NB cannot be used "
                 "in conjunction with the --cutoff option)"
                 % _DEFAULTS.CUTOFF,
                 group="Batch options")
    p.add_option('--genes',action='store',dest="genes",nargs="+",
                 metavar="GENES_FILE",
                 help="Specify multiple genes files (if used then "
                 "peaks file(s) must be specified using --peaks "
                 "option)",
                 group="Batch options")
    p.add_option('--peaks',action='store',dest="peaks",nargs="+",
                 metavar="PEAKS_FILE",
                 help="Specify multiple peaks files (if used then "
                 "genes file(s) must be specified using --genes "
                 "option)",
                 group="Batch options")
    p.add_option('-n','--nprocessors',action='store',
                 type=int,dest='nprocs',default=1,
                 help="Number of processors/cores to run the "
                 "program using (default: 1)",
                 group="Batch options")
    p.add_option('--split-outputs',action='store_true',
                 dest='multiple_outputs',
                 help="In batch mode write results of each "
                 "analysis to separate file (default is to "
                 "write all results to single file)",
                 group="Batch options")

    p.add_option_group("Advanced options")
    p.add_analyses_option(group="Advanced options")
    p.add_feature_option(group="Advanced options")
    p.add_peak_cols_option(group="Advanced options")
    p.add_peak_id_option(group="Advanced options")

    # Process command line
    options,args = p.parse_args()

    # Input files
    peak_files = []
    gene_files = []
    if options.peaks:
        peak_files = [f for f in options.peaks]
    if options.genes:
        gene_files = [f for f in options.genes]
    if peak_files and gene_files:
        if len(args) > 0:
            p.error("too many arguments: files already supplied via "
                    "--peaks and --genes")
    elif peak_files:
        p.error("must supply GENES file(s) via --genes when using "
                "--peaks option")
    elif gene_files:
        p.error("must supply PEAKS file(s) via --peaks when using "
                "--genes option")
    elif len(args) == 2:
        gene_files.append(args[0])
        peak_files.append(args[1])
    else:
        p.error("need to supply genes and peaks via command line "
                "or via --genes and --peaks options")

    # Report version and authors
    print(p.get_version())
    print(_PROGRAM_INFO)

    # Process cutoffs
    if options.max_distance and options.cutoffs:
        p.error("can't use --cutoff and --cutoffs together")
    if options.max_distance:
        try:
            cutoffs = [int(options.max_distance),]
        except ValueError:
            p.error("bad cutoff value: '%s'" % options.max_distance)
    elif options.cutoffs:
        cutoffs = list()
        for cutoff in str(options.cutoffs).split(','):
            try:
                cutoffs.append(int(cutoff))
            except ValueError:
                p.error("Bad cutoff value: '%s'" % cutoff)
    else:
        cutoffs = [_DEFAULTS.CUTOFF,]
    cutoffs.sort()

    # Deal with zero cutoff distance meaning 'no cutoff'
    cutoffs = [d if d != 0 else None for d in cutoffs]

    # Gene edge to use
    if options.edge == 'tss':
        tss_only = True
    else:
        tss_only = False

    # Promoter region
    promoter = (abs(int(options.promoter_region.split(',')[0])),
                abs(int(options.promoter_region.split(',')[1])))

    # Columns for peak data
    if options.peak_cols is None:
        peak_cols = (1,2,3)
    else:
        try:
            peak_cols = tuple([int(x)
                               for x in options.peak_cols.split(',')])
        except Exception as ex:
            logging.fatal("Bad column assignment for --peak_cols")
            sys.exit(1)

    # Use IDs for peaks?
    try:
        peak_id_col = int(options.peak_id)
    except TypeError:
        peak_id_col = None

    # Set up reporting
    if options.compact:
        mode = output.SINGLE_LINE
        placeholder = '.'
        if peak_id_col is None:
            peak_fields = ('peak.chr','peak.start','peak.end',
                           'list(feature.id,strand,TSS,TES,dist_closest,'
                           'dist_TSS,dist_TES,direction,overlap_feature,'
                           'overlap_promoter)')
            gene_fields = ('feature.id','feature.chr','feature.start',
                           'feature.end','feature.strand',
                           'list(chr,start,end,dist_closest,dist_TSS,'
                           'direction,in_the_feature)')
        else:
            peak_fields = ('peak.id','peak.chr','peak.start','peak.end',
                           'list(feature.id,strand,TSS,TES,dist_closest,'
                           'dist_TSS,dist_TES,direction,overlap_feature,'
                           'overlap_promoter)')
            gene_fields = ('feature.id','feature.chr','feature.start',
                           'feature.end','feature.strand',
                           'list(peak.id,chr,start,end,dist_closest,dist_TSS,'
                           'direction,in_the_feature)')
        if options.summary:
            options.summary = False
            logging.error("--summary option not compatible with --compact")
            sys.exit(1)
    else:
        mode = output.MULTI_LINE
        placeholder = '---'
        if peak_id_col is None:
            peak_fields = ('peak.chr','peak.start','peak.end',
                           'feature.id','strand','TSS','TES',
                           'dist_closest','dist_TSS','dist_TES',
                           'direction','overlap_feature','overlap_promoter')
            gene_fields = ('feature.id','feature.chr','feature.start',
                           'feature.end','feature.strand',
                           'chr','start','end','dist_closest','dist_TSS',
                           'direction','in_the_feature')
        else:
            peak_fields = ('peak.id','peak.chr','peak.start','peak.end',
                           'feature.id','strand','TSS','TES',
                           'dist_closest','dist_TSS','dist_TES',
                           'direction','overlap_feature','overlap_promoter')
            gene_fields = ('feature.id','feature.chr','feature.start',
                           'feature.end','feature.strand',
                           'peak.id','chr','start','end',
                           'dist_closest','dist_TSS','direction',
                           'in_the_feature')

    # Split output files
    if options.multiple_outputs:
        multiple_outputs = (len(gene_files) > 1) or \
                           (len(peak_files) > 1) or \
                           (len(cutoffs) > 1)
    else:
        multiple_outputs = False

    # Update fields for batch mode
    if len(cutoffs) > 1:
        peak_fields = tuple(['cutoff']+list(peak_fields))
        gene_fields = tuple(['cutoff']+list(gene_fields))
    if len(gene_files) > 1:
        peak_fields = tuple(['feature_file']+list(peak_fields))
        gene_fields = tuple(['feature_file']+list(gene_fields))
    if len(peak_files) > 1:
        peak_fields = tuple(['peak_file']+list(peak_fields))
        gene_fields = tuple(['peak_file']+list(gene_fields))

    # Analyses to run
    peak_centric = (options.analyses in ("all","peak_centric",))
    gene_centric = (options.analyses in ("all","gene_centric",))

    # Report inputs
    print("Genes files    : %s" % gene_files[0])
    for gene_file in gene_files[1:]:
        print("                 %s" % gene_file)
    print("Peaks files    : %s" % peak_files[0])
    for peak_file in peak_files[1:]:
        print("                 %s" % peak_file)
    print("Cutoffs (bp)   : %s" % ','.join([str(d) if d is not None
                                            else "no cutoff"
                                            for d in cutoffs]))
    print("Edge           : %s" % ('TSS only' if tss_only
                                   else 'TSS or TES'))
    print("DE only        : %s" % ('yes' if options.only_diff_expressed
                                   else 'no'))
    print("Nprocs         : %s" % options.nprocs)
    print("Max no. of hits: %s" % ('All' if options.max_closest is None
                                   else "%d" % options.max_closest))
    print("Promoter region: -%d to %d (bp from TSS)" % promoter)
    print("Feature type   : %s" % options.feature_type)
    print("")
    print("Analyses:")
    print("- Peak-centric: %s" % ('yes' if peak_centric else 'no'))
    print("- Gene-centric: %s" % ('yes' if gene_centric else 'no'))

    # Read in gene data
    gene_lists = dict()
    for gene_file in gene_files:
        genes = read_feature_file(gene_file)
        if options.only_diff_expressed and not genes.isFlagged():
            logging.fatal("--only-DE flag needs input genes flagged as "
                          "differentially expressed")
            sys.exit(1)
        gene_lists[gene_file] = genes

    # Read in peak data
    peak_lists = dict()
    for peak_file in peak_files:
        peak_lists[peak_file] = read_peak_file(peak_file,
                                               peak_cols=peak_cols,
                                               peak_id_col=peak_id_col)

    # Output files
    if options.name is not None:
        basename = options.name
    else:
        basename = os.path.splitext(os.path.basename(gene_files[0]))[0]
    if multiple_outputs:
        # Split outputs into multiple files
        namer = BatchNamer(basename,peak_files,gene_files,cutoffs)
        for peak_file in peak_files:
            for gene_file in gene_files:
                for cutoff in cutoffs:
                    name = namer.get_name(peak_file,gene_file,cutoff)
                    # Remove existing files
                    OutputFiles(name).remove_files()
    else:
        # All outputs into same files
        OutputFiles(basename).remove_files()

    # Assemble inputs over cutoffs and peak and gene lists
    # The same parameters can be used for both peak- and
    # gene-centric analyses
    analysis_params = []
    for peaks in peak_files:
        for genes in gene_files:
            for cutoff in cutoffs:
                analysis_params.append(
                    AnalysisParams(
                        genes=gene_lists[genes],
                        peaks=peak_lists[peaks],
                        cutoff=cutoff,
                        tss_only=tss_only,
                        only_differentially_expressed=
                        options.only_diff_expressed
                    ))

    # Run the peak-centric analyses
    if peak_centric:
        print("**** Peak-centric analysis: nearest genes to each peak ****")
        # Run the analyses
        if options.nprocs > 1:
            # Multiple cores
            pool = Pool(options.nprocs)
            # Version of multiprocessing which can also
            # handle ctrl-C terminating the program
            # See http://bryceboe.com/2010/08/26/python-multiprocessing-and-keyboardinterrupt/
            p = pool.map_async(find_nearest_features_as_list,
                               analysis_params)
            try:
                results = p.get(0xFFFF)
            except KeyboardInterrupt:
                print("KeyboardInterrupt")
                sys.exit(1)
        else:
            # Single core
            results = map(lambda p: list(find_nearest_features(p)),
                          analysis_params)
        # Output the results
        peak_centric_outputs = dict()
        peak_centric_summary = dict()
        append_outputs = False
        for result in results:
            # Deal with output files
            if multiple_outputs:
                peak,nearest_genes,params = result[0]
                name = namer.get_name(params.peaks.source_file,
                                      params.genes.source_file,
                                      params.cutoff)
            else:
                name = basename
            outputs = OutputFiles(name)
            # Set up reporter
            reporter = output.AnalysisReportWriter(
                mode,peak_fields,
                promoter_region=promoter,
                null_placeholder=placeholder,
                max_hits=options.max_closest,
                pad=options.pad_output,
                outfile=outputs.peak_centric_out,
                summary=(outputs.peak_centric_summary
                         if options.summary else None),
                feature_type=options.feature_type,
                append=append_outputs)
            # Write results
            for peak,nearest_genes,params in result:
                reporter.write_nearest_features(
                    peak,nearest_genes,
                    peak_file=params.peaks.source_file,
                    feature_file=params.genes.source_file,
                    cutoff=params.cutoff)
            if not append_outputs:
                print("Results written to %s" % outputs.peak_centric_out)
                peak_centric_outputs[name] = outputs.peak_centric_out
                if options.summary:
                    print("Summary written to %s" %
                          outputs.peak_centric_summary)
                    peak_centric_summary[name] = outputs.peak_centric_summary
            print("")
            if not multiple_outputs and not append_outputs:
                append_outputs = True
            # Close reporter
            reporter.close()

    # Run the gene/feature-centric analyses
    if gene_centric:
        print("**** Gene-centric analysis: nearest peaks to each gene ****")
        # Run the analyses
        if options.nprocs > 1:
            # Multiple cores
            pool = Pool(options.nprocs)
            # Version of multiprocessing which can also
            # handle ctrl-C terminating the program
            # See http://bryceboe.com/2010/08/26/python-multiprocessing-and-keyboardinterrupt/
            p = pool.map_async(find_nearest_peaks_as_list,
                               analysis_params)
            try:
                results = p.get(0xFFFF)
            except KeyboardInterrupt:
                print("KeyboardInterrupt")
                sys.exit(1)
        else:
            # Single core
            results = map(lambda p: list(find_nearest_peaks(p)),
                          analysis_params)
        # Output the results
        gene_centric_outputs = dict()
        gene_centric_summary = dict()
        append_outputs = False
        for result in results:
            # Deal with output files
            if multiple_outputs:
                gene,nearest_peaks,params = result[0]
                name = namer.get_name(params.peaks.source_file,
                                      params.genes.source_file,
                                      params.cutoff)
            else:
                name = basename
            outputs = OutputFiles(name)
            # Set up reporter
            reporter = output.AnalysisReportWriter(
                mode,gene_fields,
                null_placeholder=placeholder,
                max_hits=options.max_closest,
                pad=options.pad_output,
                outfile=outputs.gene_centric_out,
                summary=(outputs.gene_centric_summary
                         if options.summary else None),
                feature_type=options.feature_type,
                append=append_outputs)
            # Write results
            for gene,nearest_peaks,params in result:
                reporter.write_nearest_features(
                    gene,nearest_peaks,
                    peak_file=params.peaks.source_file,
                    feature_file=params.genes.source_file,
                    cutoff=params.cutoff)
            if not append_outputs:
                print("Results written to %s" % outputs.gene_centric_out)
                gene_centric_outputs[name] = outputs.gene_centric_out
                if options.summary:
                    print("Summary written to %s" %
                          outputs.gene_centric_summary)
                    gene_centric_summary[name] = outputs.gene_centric_summary
            print("")
            if not multiple_outputs and not append_outputs:
                append_outputs = True
            # Close reporter
            reporter.close()

    # Make XLSX file
    if options.xlsx_output:
        print("**** Writing XLSX file ****")
        xlsx = xls_output.XLSX(OutputFiles(basename).xlsx_out,
                               p.get_version(),
                               options.feature_type)
        # Write the settings
        xlsx.append_to_notes("Input %ss file\t%s" % (options.feature_type,
                                                     gene_files[0]))
        for gene_file in gene_files[1:]:
            xlsx.append_to_notes("\t%s" % gene_file)
        xlsx.append_to_notes("Input peaks file\t%s" % peak_files[0])
        for peak_file in peak_files[1:]:
            xlsx.append_to_notes("\t%s" % peak_file)
        cutoff_distances = [str(d) if d != 0 else 'no cutoff'
                            for d in cutoffs]
        xlsx.append_to_notes("Cutoff distances (bp)\t%s" %
                             ','.join([str(d) for d in cutoff_distances]))
        xlsx.append_to_notes("Maximum no. of hits to report\t%s"
                             % ('All' if options.max_closest is None
                                else "%d" % options.max_closest))
        xlsx.append_to_notes("Promoter region (bp from TSS)\t-%d to %d" %
                             promoter)
        if tss_only:
            xlsx.append_to_notes("Distances calculated from\tTSS only")
        else:
            xlsx.append_to_notes("Distances calculated from\tTSS or TES")
        xlsx.append_to_notes("Only use differentially expressed %ss\t%s" %
                             (options.feature_type,
                              "Yes" if options.only_diff_expressed
                              else "No"))
        # Add features to peaks
        if peak_centric:
            names = sorted(peak_centric_outputs.keys())
            xlsx.write_peak_centric(peak_fields)
            if options.summary:
                xlsx.append_to_notes("\n'Peak-centric (summary)' lists the "
                                     "'top' result (i.e. closest peak/%s "
                                     "pair) for each peak" %
                                     options.feature_type)
            for i,name in enumerate(names):
                title = "Peak-centric"
                if multiple_outputs:
                    title += " #%d" % (i+1)
                xlsx.add_result_sheet(title,peak_centric_outputs[name])
                if options.summary:
                    title = "Peak-centric (summary)"
                    if multiple_outputs:
                        title += " #%d" % (i+1)
                    xlsx.add_result_sheet(title,
                                          peak_centric_summary[name])
        # Add peaks to features
        if gene_centric:
            names = sorted(gene_centric_outputs.keys())
            xlsx.write_feature_centric(gene_fields)
            if options.summary:
                xlsx.append_to_notes("\n'%s-centric (summary)' lists the "
                                     "'top' result (i.e. closest %s/peak "
                                     "pair) for each %s" %
                                     (options.feature_type.title(),
                                      options.feature_type,
                                      options.feature_type))
            for i,name in enumerate(names):
                title = "%s-centric" % options.feature_type.title()
                if multiple_outputs:
                    title += " #%d" % (i+1)
                xlsx.add_result_sheet(title,gene_centric_outputs[name])
                if options.summary:
                    title = "%s-centric (summary)" % \
                            options.feature_type.title()
                    if multiple_outputs:
                        title += " #%d" % (i+1)
                    xlsx.add_result_sheet(title,
                                          gene_centric_summary[name])
        xlsx.write()
        print("Wrote %s" % xlsx.xlsx_file)
        print("")

    # Finished
    print("Done")
