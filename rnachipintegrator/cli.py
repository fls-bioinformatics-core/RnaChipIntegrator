#!/usr/bin/env python
#
#     RnaChipIntegrator.py: analyse genomic features (genes) with peak data
#     Copyright (C) University of Manchester 2011-16 Peter Briggs, Leo Zeef
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
import time,datetime
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

DEFAULT_FEATURE_TYPE = 'gene'

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
    
    # Defaults
    promoter = (1000,100)
    max_distance = 1000000
    
    # Parse command line
    p = optparse.OptionParser(usage="%prog [options] GENES PEAKS",
                              version="%prog "+__version__,
                              description=
                              "Analyse GENES (any set of genes or genomic "
                              "features) against PEAKS (a set of regions) "
                              "and report nearest genes to each peak (and "
                              "vice versa)")

    # Analysis options
    analysis_opts = optparse.OptionGroup(p,"Analysis options")
    analysis_opts.add_option('--cutoff',action='store',dest='max_distance',
                             type='int',default=max_distance,
                             help="Maximum distance allowed between peaks "
                             "and genes before being omitted from the "
                             "analyses (default %dbp; set to zero for no "
                             "cutoff)" % max_distance)
    analysis_opts.add_option('--edge',action='store',dest="edge",
                             type="choice",choices=('tss','both'),
                             default='tss',
                             help="Gene edges to consider when calculating "
                             "distances between genes and peaks, either: "
                             "'tss' (default: only use gene TSS) or 'both' "
                             "(use whichever of TSS or TES gives shortest "
                             "distance)")
    analysis_opts.add_option('--only-DE',action='store_true',
                             dest='only_diff_expressed',default=False,
                             help="Only use genes flagged as differentially "
                             "expressed in analyses (input gene data must "
                             "include DE flags)")
    p.add_option_group(analysis_opts)

    # Reporting options
    reporting_opts = optparse.OptionGroup(p,"Reporting options")
    reporting_opts.add_option('--number',action='store',dest='max_closest',
                              type='int',default=None,
                              help="Filter results after applying --cutoff "
                              "to report only the nearest MAX_CLOSEST number "
                              "of pairs for each gene/peak from the analyses "
                              "(default is to report all results)")
    reporting_opts.add_option('--promoter_region',action="store",
                              dest="promoter_region",
                              default="%d,%d" % promoter,
                              help="Define promoter region with respect to "
                              "gene TSS in the form UPSTREAM,DOWNSTREAM "
                              "(default -%d to %dbp of TSS)" %  promoter)
    p.add_option_group(reporting_opts)

    # Output options
    output_opts = optparse.OptionGroup(p,"Output options")
    output_opts.add_option('--name',action='store',dest='name',default=None,
                           help="Set basename for output files")
    output_opts.add_option('--compact',action='store_true',dest='compact',
                           default=False,
                           help="Output all hits for each peak or gene on a "
                           "single line (cannot be used with --summary)")
    output_opts.add_option('--summary',action='store_true',dest='summary',
                           default=False,
                           help="Output 'summary' for each analysis, "
                           "consisting of only the top hit for each peak or "
                           "gene (cannot be used with --compact)")
    output_opts.add_option('--pad',action="store_true",dest="pad_output",
                           help="Where less than MAX_CLOSEST hits are found, "
                           "pad output with blanks to ensure that MAX_CLOSEST "
                           "hits are still reported (nb --pad is implied for "
                           "--compact)")
    output_opts.add_option('--xlsx',action="store_true",dest="xlsx_output",
                           help="Output XLSX spreadsheet with results")
    p.add_option_group(output_opts)

    # Advanced options
    advanced_opts = optparse.OptionGroup(p,"Advanced options")
    advanced_opts.add_option('--feature',action="store",dest="feature_type",
                             help="rename '%s' to FEATURE_TYPE in output (e.g. "
                             "'transcript' etc)" % DEFAULT_FEATURE_TYPE)
    advanced_opts.add_option('--peak_cols',action="store",dest="peak_cols",
                             help="list of 3 column indices (e.g. '1,4,5') "
                             "indicating columns to use for chromosome, "
                             "start and end from the input peak file (if not "
                             "first three columns).")
    p.add_option_group(advanced_opts)

    # Process command line
    options,args = p.parse_args()

    # Input files
    if len(args) != 2:
        p.error("need to supply 2 files (genes and peaks)")
    gene_file,peak_file = args

    # Report version and authors
    p.print_version()
    print
    print "Find nearest peaks to genes (and vice versa)"
    print
    print "University of Manchester"
    print "Faculty of Life Sciences"
    print "Bioinformatics Core Facility"
    print "Authors: Peter Briggs, Ian Donaldson and Leo Zeef"
    print
    print "If you use this program in your published work then please cite:"
    print
    print "   Briggs PJ, Donaldson IJ, Zeef LAH. RnaChipIntegrator"
    print "   (version %s). Available at:" % __version__
    print "   https://github.com/fls-bioinformatics-core/RnaChipIntegrator"
    print

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

    # Reporting formats
    if options.compact:
        mode = output.SINGLE_LINE
        peak_fields = ('peak.chr','peak.start','peak.end',
                       'list(feature.id,strand,TSS,TES,dist_closest,'
                       'dist_TSS,dist_TES,direction,overlap_feature,'
                       'overlap_promoter)')
        gene_fields = ('feature.id','feature.chr','feature.start',
                       'feature.end','feature.strand',
                       'list(chr,start,end,dist_closest,dist_TSS,'
                       'direction,in_the_feature)')
        placeholder = '.'
        if options.summary:
            options.summary = False
            logging.error("--summary option not compatible with --compact")
            sys.exit(1)
    else:
        mode = output.MULTI_LINE
        peak_fields = ('peak.chr','peak.start','peak.end',
                       'feature.id','strand','TSS','TES',
                       'dist_closest','dist_TSS','dist_TES',
                       'direction','overlap_feature','overlap_promoter')
        gene_fields = ('feature.id','feature.chr','feature.start',
                       'feature.end','feature.strand',
                       'chr','start','end','dist_closest','dist_TSS',
                       'direction','in_the_feature')
        placeholder = '---'

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
        except Exception, ex:
            p.error("Bad column assignment for --peak_cols")

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
    if tss_only:
        print "Distances will be calculated from gene TSS only"
    else:
        print "Distances will be calculated from nearest of gene TSS or TES"
    print
    print "Genomic features are '%ss'" % feature_type
    print

    # Read in gene data
    try:
        genes = FeatureSet(gene_file)
    except Exception, ex:
        logging.critical("Failed to read in gene data: %s" % ex)
        print "Please fix errors in input file before running again"
        sys.exit(1)
    if not len(genes):
        logging.error("No gene data read in")
        sys.exit(1)
    print "%d gene records read in" % len(genes)

    # Differential expression handling
    use_differentially_expressed = False
    if genes.isFlagged():
        print "\tGene data include differential expression flag"
        print "\t%d genes flagged as differentially expressed" % \
            len(genes.filterByFlag(1))
        if options.only_diff_expressed:
            print
            print "*** Only differentially expressed genes will used ***"
            use_differentially_expressed = True
    elif options.only_diff_expressed:
        logging.error("--only-DE flag needs input genes flagged as "
                      "differentially expressed")
        sys.exit(1)
    print

    # Read in peak data
    print "Using columns %s from peaks file as chrom, start, end" % \
        (peak_cols,)
    try:
        peaks = PeakSet(peak_file,columns=peak_cols)
    except Exception,ex:
        logging.critical("Failed to read peak data (%s)" % ex)
        print "Please fix errors in input file before running again"
        sys.exit(1)
    if not len(peaks):
        logging.error("No peak data read in")
        sys.exit(1)
    print "%d peak records read in" % len(peaks)
    if peaks.isSummit():
        print "\tPeak data are summits"
    else:
        print "\tPeak data are regions"
    print

    # Set up output files
    if options.name is not None:
        basename = options.name
    else:
        basename = os.path.splitext(os.path.basename(gene_file))[0]
    peak_centric_out = basename+"_peak_centric.txt"
    gene_centric_out = basename+"_gene_centric.txt"
    if options.summary:
        peak_centric_summary = basename+"_peak_centric_summary.txt"
        gene_centric_summary = basename+"_gene_centric_summary.txt"
    else:
        peak_centric_summary = None
        gene_centric_summary = None
    if options.xlsx_output:
        xlsx_out = basename+".xlsx"
    else:
        xlsx_out = None

    # Clean up any pre-existing output files that would otherwise
    # be overwritten
    print "Checking for pre-existing output files"
    for f in (peak_centric_out,peak_centric_summary,
              gene_centric_out,gene_centric_summary,
              xlsx_out):
        if f is not None and os.path.isfile(f):
            print "\tRemoving %s" % f
            os.remove(f)
    print

    # Do the analyses
    print "**** Peak-centric analysis: nearest genes to each peak ****"
    reporter = output.AnalysisReportWriter(mode,peak_fields,
                                           promoter_region=promoter,
                                           null_placeholder=placeholder,
                                           max_hits=max_closest,
                                           pad=options.pad_output,
                                           outfile=peak_centric_out,
                                           summary=peak_centric_summary,
                                           feature_type=feature_type)
    for peak,nearest_genes in analysis.find_nearest_features(
            peaks,genes,tss_only=tss_only,distance=max_distance,
            only_differentially_expressed=use_differentially_expressed):
        reporter.write_nearest_features(peak,nearest_genes)
    reporter.close()
    print "Results written to %s" % peak_centric_out
    if peak_centric_summary:
        print "Summary written to %s" % peak_centric_summary
    print

    print "**** Gene-centric analysis: nearest peaks to each gene ****"
    reporter = output.AnalysisReportWriter(mode,gene_fields,
                                           null_placeholder=placeholder,
                                           max_hits=max_closest,
                                           pad=options.pad_output,
                                           outfile=gene_centric_out,
                                           summary=gene_centric_summary,
                                           feature_type=feature_type)
    for gene,nearest_peaks in analysis.find_nearest_peaks(
            genes,peaks,tss_only=tss_only,distance=max_distance,
            only_differentially_expressed=use_differentially_expressed):
        reporter.write_nearest_peaks(gene,nearest_peaks)
    reporter.close()
    print "Results written to %s" % gene_centric_out
    if gene_centric_summary:
        print "Summary written to %s" % gene_centric_summary
    print

    # Make XLSX file
    if options.xlsx_output:
        print "**** Writing XLSX file ****"
        xlsx = xls_output.XLSX(xlsx_out,p.get_version(),feature_type)
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
                              "Yes" if use_differentially_expressed else "No"))
        # Add features to peaks
        xlsx.write_peak_centric(peak_fields)
        xlsx.add_result_sheet('Peak-centric',peak_centric_out)
        if options.summary:
            xlsx.append_to_notes("\n'Peak-centric (summary)' lists the 'top' "
                                 "result (i.e. closest peak/%s pair) for each "
                                 "peak" % feature_type)
            xlsx.add_result_sheet('Peak-centric (summary)',
                                  peak_centric_summary)
        # Add peaks to features
        xlsx.write_feature_centric(gene_fields)
        xlsx.add_result_sheet('%s-centric' % feature_type.title(),
                              gene_centric_out)
        if options.summary:
            xlsx.append_to_notes("\n'%s-centric (summary)' lists the 'top' "
                                 "result (i.e. closest %s/peak pair) for each "
                                 "%s" % (feature_type.title(),
                                         feature_type,
                                         feature_type))
            xlsx.add_result_sheet('%s-centric (summary)' %
                                  feature_type.title(),
                                  gene_centric_summary)
        xlsx.write()
        print "Wrote %s" % xlsx_out
        print

    # Finished
    print "Done"
