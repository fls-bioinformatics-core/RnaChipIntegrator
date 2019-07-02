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
import optparse
from .Features import FeatureSet
from .Peaks import PeakSet
import analysis
import output
import xls_output
import logging
logging.getLogger().setLevel(logging.WARNING)
logging.basicConfig(format='%(levelname)s: %(message)s')

from . import get_version
__version__ = get_version()

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
    >>> print options.cutoff
    >>> print options.force

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
        if version is None:
            version = "%prog "+__version__
        self.parser = optparse.OptionParser(
            usage=usage,
            version=version,
            description=description)
        self.option_groups = dict()

    def get_version(self):
        """
        Return the version string from the parser
        """
        return self.parser.get_version()

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
        return self.parser.parse_args(args)

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
        group = optparse.OptionGroup(self.parser,name)
        self.parser.add_option_group(group)
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
        p.add_option(*args,**kws)

    def add_edge_option(self,group=None):
        """
        Add --edge option to the parser

        Acceptable values are 'tss' (default) or 'both'.
        The value is accessed via the 'edge' property of
        the parser options.
        """
        self.add_option('--edge',
                        action='store',dest="edge",
                        type="choice",choices=('tss','both'),
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
                        type='int',default=_DEFAULTS.MAX_CLOSEST,
                        help="Filter results after applying --cutoff "
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
                        type='choice',default="all",
                        choices=('all','gene_centric','peak_centric',),
                        help="Select which analyses to run: can be one "
                        "of 'all' (default, runs all available "
                        "analyses), 'peak_centric' or 'gene_centric'")

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
                        type='int',
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
        print "Removing pre-existing output files"
        for f in self.files:
            if os.path.isfile(f):
                print "\tRemoving %s" % f
                os.remove(f)
        print

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

    print "Reading in data from '%s'" % feature_file
    try:
        features = FeatureSet(feature_file)
    except Exception as ex:
        logging.fatal("Failed to read in %s data: %s" % (feature_type,
                                                         ex))
        print "Please fix errors in input file before running again"
        sys.exit(1)
    if not len(features):
        logging.fatal("No %s data read in" % feature_type)
        sys.exit(1)
    print "%d %s records read in" % (len(features),
                                     feature_type)
    if features.isFlagged():
        print "\tData include differential expression flag"
        print "\t%d records flagged as differentially expressed" % \
            len(features.filterByFlag(1))
    print
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
    print "Reading in data from '%s'" % peak_file
    print "Using columns %s from peaks file as chrom, start, end" % \
        (peak_cols,)
    if peak_id_col is not None:
        peak_id_col = int(peak_id_col)
        print "Using column %s from peaks file as peak ID" % peak_id_col
    try:
        peaks = PeakSet(peak_file,columns=peak_cols,
                        id_column=peak_id_col)
    except Exception as ex:
        logging.fatal("Failed to read peak data (%s)" % ex)
        print "Please fix errors in input file before running again"
        sys.exit(1)
    if not len(peaks):
        logging.error("No peak data read in")
        sys.exit(1)
    print "%d peak records read in (%s)" % (len(peaks),
                                            'summits' if peaks.isSummit()
                                            else 'regions')
    print
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

    p = CLI(usage="%prog [options] GENES PEAKS",
            description=
            "Analyse GENES (any set of genes or genomic "
            "features) against one or more sets of PEAKS "
            "(a set of regions) and report nearest genes "
            "to each peak (and vice versa)")

    # Analysis options
    p.add_option_group("Analysis options")
    p.add_option('--cutoff',action='store',dest='max_distance',
                 type='int',default=_DEFAULTS.CUTOFF,
                 help="Maximum distance allowed between peaks "
                 "and genes before being omitted from the "
                 "analyses (default %dbp; set to zero for no "
                 "cutoff)" % _DEFAULTS.CUTOFF,
                 group="Analysis options")
    p.add_edge_option(group="Analysis options")
    p.add_only_de_option(group="Analysis options")

    # Reporting options
    p.add_option_group("Reporting options")
    p.add_number_option(group="Reporting options")
    p.add_promoter_region_option(group="Reporting options")

    # Output options
    p.add_option_group("Output options")
    p.add_name_option(group="Output options")
    p.add_compact_option(group="Output options")
    p.add_summary_option(group="Output options")
    p.add_pad_option(group="Output options")
    p.add_xlsx_option(group="Output options")

    # Advanced options
    p.add_option_group("Advanced options")
    p.add_option('--analyses',action='store',dest="analyses",
                 type='choice',default="all",
                 choices=('all','gene_centric','peak_centric',),
                 help="Select which analyses to run: can be one "
                 "of 'all' (default, runs all available "
                 "analyses), 'peak_centric' or 'gene_centric'",
                 group="Advanced options")
    p.add_feature_option(group="Advanced options")
    p.add_peak_cols_option(group="Advanced options")
    p.add_peak_id_option(group="Advanced options")

    # Process command line
    options,args = p.parse_args()

    # Input files
    if len(args) != 2:
        p.error("need to supply 2 files (genes and peaks)")
    gene_file,peak_file = args

    # Report version and authors
    print p.get_version()
    print _PROGRAM_INFO

    # Promoter region
    promoter = (abs(int(options.promoter_region.split(',')[0])),
                abs(int(options.promoter_region.split(',')[1])))

    # Reporting options
    max_distance = options.max_distance
    if max_distance <= 0:
        max_distance = None
    max_closest = options.max_closest

    # Gene edge to use
    if options.edge == 'tss':
        tss_only = True
    else:
        tss_only = False

    # Feature type
    if options.feature_type is None:
        feature_type = 'gene'
    else:
        feature_type = options.feature_type

    # Columns to extract from input peaks file
    if options.peak_cols is None:
        peak_cols = (1,2,3)
    else:
        try:
            peak_cols = tuple([int(x)
                               for x in options.peak_cols.split(',')])
        except Exception as ex:
            p.error("Bad column assignment for --peak_cols")

    # Handle peak IDs, if specified
    if options.peak_id is not None:
        peak_id_col = int(options.peak_id)
    else:
        peak_id_col = None

    # Reporting formats
    if options.compact:
        mode = output.SINGLE_LINE
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
        placeholder = '.'
        if options.summary:
            options.summary = False
            logging.error("--summary option not compatible with --compact")
            sys.exit(1)
    else:
        mode = output.MULTI_LINE
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
                           'peak.id','chr','start','end','dist_closest','dist_TSS',
                           'direction','in_the_feature')
        placeholder = '---'

    # Analyses to run
    peak_centric = (options.analyses in ("all","peak_centric",))
    gene_centric = (options.analyses in ("all","gene_centric",))

    # Report settings
    print "Input genes file: %s" % gene_file
    print "Input peaks file: %s" % peak_file
    print
    if max_distance is not None:
        print "Maximum cutoff distance: %d (bp)" % max_distance
    else:
        print "Maximum cutoff distance: no cutoff"
    print "Maximum no. of hits    : %s" % ('All' if max_closest is None
                                           else "%d" % max_closest)
    print "Promoter region        : -%d to %d (bp from TSS)" % promoter
    print
    print "Analyses:"
    print "- Peak-centric: %s" % ('yes' if peak_centric else 'no')
    print "- Gene-centric: %s" % ('yes' if gene_centric else 'no')
    if tss_only:
        print "Distances will be calculated from gene TSS only"
    else:
        print "Distances will be calculated from nearest of gene TSS or TES"
    print
    print "Genomic features are '%ss'" % feature_type
    print

    # Read in gene data
    genes = read_feature_file(gene_file)
    if options.only_diff_expressed:
        if not genes.isFlagged():
            logging.fatal("--only-DE flag needs input genes flagged as "
                          "differentially expressed")
            sys.exit(1)
        else:
            print
            print "*** Only differentially expressed genes will used ***"
    print

    # Read in peak data
    peaks = read_peak_file(peak_file,
                           peak_cols=peak_cols,
                           peak_id_col=peak_id_col)

    # Set up output files
    if options.name is not None:
        basename = options.name
    else:
        basename = os.path.splitext(os.path.basename(gene_file))[0]
    outputs = OutputFiles(basename)

    # Clean up any pre-existing output files that would otherwise
    # be overwritten
    outputs.remove_files()

    # Do the analyses
    if peak_centric:
        print "**** Peak-centric analysis: nearest genes to each peak ****"
        reporter = output.AnalysisReportWriter(
            mode,peak_fields,
            promoter_region=promoter,
            null_placeholder=placeholder,
            max_hits=max_closest,
            pad=options.pad_output,
            outfile=outputs.peak_centric_out,
            summary=(outputs.peak_centric_summary
                     if options.summary else None),
            feature_type=feature_type)
        params = AnalysisParams(genes=genes,
                                peaks=peaks,
                                cutoff=max_distance,
                                tss_only=tss_only,
                                only_differentially_expressed=
                                options.only_diff_expressed)
        for peak,nearest_genes,prms in find_nearest_features(params):
            reporter.write_nearest_features(peak,nearest_genes)
        reporter.close()
        print "Results written to %s" % outputs.peak_centric_out
        if options.summary:
            print "Summary written to %s" % peak_centric_summary
        print

    if gene_centric:
        print "**** Gene-centric analysis: nearest peaks to each gene ****"
        reporter = output.AnalysisReportWriter(
            mode,gene_fields,
            null_placeholder=placeholder,
            max_hits=max_closest,
            pad=options.pad_output,
            outfile=outputs.gene_centric_out,
            summary=(outputs.gene_centric_summary
                     if options.summary else None),
            feature_type=feature_type)
        params = AnalysisParams(genes=genes,
                                peaks=peaks,
                                cutoff=max_distance,
                                tss_only=tss_only,
                                only_differentially_expressed=
                                options.only_diff_expressed)
        for gene,nearest_peaks,prms in find_nearest_peaks(params):
            reporter.write_nearest_peaks(gene,nearest_peaks)
        reporter.close()
        print "Results written to %s" % outputs.gene_centric_out
        if options.summary:
            print "Summary written to %s" % outputs.gene_centric_summary
        print

    # Make XLSX file
    if options.xlsx_output:
        print "**** Writing XLSX file ****"
        xlsx = xls_output.XLSX(outputs.xlsx_out,
                               p.get_version(),
                               feature_type)
        # Write the settings
        xlsx.append_to_notes("Input %ss file\t%s" % (feature_type,
                                                     gene_file))
        xlsx.append_to_notes("Input peaks file\t%s" % peak_file)
        if max_distance is not None:
            xlsx.append_to_notes("Maximum cutoff distance (bp)\t%d" %
                                 max_distance)
        else:
            xlsx.append_to_notes("Maximum cutoff distance (bp)\tno cutoff")
        xlsx.append_to_notes("Maximum no. of hits to report\t%s"
                             % ('All' if max_closest is None
                                else "%d" % max_closest))
        xlsx.append_to_notes("Promoter region (bp from TSS)\t-%d to %d" %
                             promoter)
        if tss_only:
            xlsx.append_to_notes("Distances calculated from\tTSS only")
        else:
            xlsx.append_to_notes("Distances calculated from\tTSS or TES")
        xlsx.append_to_notes("Only use differentially expressed %ss\t%s" %
                             (feature_type,
                              "Yes" if options.only_diff_expressed else "No"))
        # Add features to peaks
        if peak_centric:
            xlsx.write_peak_centric(peak_fields)
            xlsx.add_result_sheet('Peak-centric',outputs.peak_centric_out)
            if options.summary:
                xlsx.append_to_notes("\n'Peak-centric (summary)' lists the "
                                     "'top' result (i.e. closest peak/%s "
                                     "pair) for each peak" % feature_type)
                xlsx.add_result_sheet('Peak-centric (summary)',
                                      peak_centric_summary)
        # Add peaks to features
        if gene_centric:
            xlsx.write_feature_centric(gene_fields)
            xlsx.add_result_sheet('%s-centric' % feature_type.title(),
                                  outputs.gene_centric_out)
            if options.summary:
                xlsx.append_to_notes("\n'%s-centric (summary)' lists the "
                                     "'top' result (i.e. closest %s/peak "
                                     "pair) for each %s" %
                                     (feature_type.title(),
                                      feature_type,
                                      feature_type))
                xlsx.add_result_sheet('%s-centric (summary)' %
                                      feature_type.title(),
                                      outputs.gene_centric_summary)
        xlsx.write()
        print "Wrote %s" % outputs.xlsx_out
        print

    # Finished
    print "Done"
