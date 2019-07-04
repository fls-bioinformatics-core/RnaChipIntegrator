#!/usr/bin/env python
#
#     batch.py: implements 'batch' version of RnaChipIntegrator
#     Copyright (C) University of Manchester 2018-2019 Peter Briggs,
#     Leo Zeef & Ian Donaldson
#
#     This code is free software; you can redistribute it and/or modify it
#     under the terms of the Artistic License 2.0 (see the file LICENSE
#     included with the distribution).
#
########################################################################
#
# batch.py
#
#########################################################################

"""
batch.py

Batch version of RnaChipIntegrator for analysing genomic features
(genes) with peak data.
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
from multiprocessing import Pool

from .cli import _DEFAULTS
from .cli import _PROGRAM_INFO
from .cli import CLI
from .cli import OutputFiles
from .cli import AnalysisParams
from .cli import read_feature_file
from .cli import read_peak_file
from .cli import find_nearest_features
from .cli import find_nearest_peaks

from . import get_version
__version__ = get_version()

logging.getLogger().setLevel(logging.WARNING)
logging.basicConfig(format='%(levelname)s: %(message)s')

#######################################################################
# Functions
#######################################################################

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
    Implements the 'RnaChipIntegrator-batch' CLI
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
    p.add_option('--cutoffs',action='store',dest='cutoffs',
                 default=None,
                 help="Comma-separated list of one or more "
                 "maximum distances allowed between peaks "
                 "and genes (bp). An analysis will be "
                 "performed for each GENES-PEAKS pair at "
                 "each cutoff distance (default %dbp; set "
                 "to zero for no cutoff)" % _DEFAULTS.CUTOFF,
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
    print p.get_version()
    print _PROGRAM_INFO

    # Process cutoffs
    if options.cutoffs is None:
        cutoffs = [_DEFAULTS.CUTOFF,]
    else:
        cutoffs = list()
        for cutoff in str(options.cutoffs).split(','):
            try:
                cutoffs.append(int(cutoff))
            except ValueError:
                logging.critical("Bad cutoff value: '%s'" % cutoff)
                sys.exit(1)
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
        peak_id_col = int(options.peak_id) - 1
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
    print "Genes files    : %s" % gene_files[0]
    for gene_file in gene_files[1:]:
        print "                 %s" % gene_file
    print "Peaks files    : %s" % peak_files[0]
    for peak_file in peak_files[1:]:
        print "                 %s" % peak_file
    print "Cutoffs (bp)   : %s" % ','.join([str(d) if d != 0
                                            else "no cutoff"
                                            for d in cutoffs])
    print "Edge           : %s" % ('TSS only' if tss_only
                                   else 'TSS or TES')
    print "DE only        : %s" % ('yes' if options.only_diff_expressed
                                   else 'no')
    print "Nprocs         : %s" % options.nprocs
    print "Max no. of hits: %s" % ('All' if options.max_closest is None
                                   else "%d" % options.max_closest)
    print "Promoter region: -%d to %d (bp from TSS)" % promoter
    print "Feature type   : %s" % options.feature_type
    print
    print "Analyses:"
    print "- Peak-centric: %s" % ('yes' if peak_centric else 'no')
    print "- Gene-centric: %s" % ('yes' if gene_centric else 'no')

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
    outputs = OutputFiles(basename)
    outputs.remove_files()

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

    # Run the analyses
    if peak_centric:
        print "**** Peak-centric analysis: nearest genes to each peak ****"
        # Build reporter
        reporter = output.AnalysisReportWriter(
            mode,peak_fields,
            promoter_region=promoter,
            null_placeholder=placeholder,
            max_hits=options.max_closest,
            pad=options.pad_output,
            outfile=outputs.peak_centric_out,
            summary=(outputs.peak_centric_summary
                     if options.summary else None),
            feature_type=options.feature_type)
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
                print "KeyboardInterrupt"
                sys.exit(1)
        else:
            # Single core
            results = map(lambda p: list(find_nearest_features(p)),
                          analysis_params)
        # Output the results
        for result in results:
            for peak,nearest_genes,params in result:
                reporter.write_nearest_features(
                    peak,nearest_genes,
                    peak_file=params.peaks.source_file,
                    feature_file=params.genes.source_file,
                    cutoff=params.cutoff)
        reporter.close()
        print "Results written to %s" % outputs.peak_centric_out
        if options.summary:
            print "Summary written to %s" % outputs.peak_centric_summary
        print

    # Run the analyses
    if gene_centric:
        print "**** Gene-centric analysis: nearest peaks to each gene ****"
        # Build reporter
        reporter = output.AnalysisReportWriter(
            mode,gene_fields,
            null_placeholder=placeholder,
            max_hits=options.max_closest,
            pad=options.pad_output,
            outfile=outputs.gene_centric_out,
            summary=(outputs.gene_centric_summary
                     if options.summary else None),
            feature_type=options.feature_type)
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
                print "KeyboardInterrupt"
                sys.exit(1)
        else:
            # Single core
            results = map(lambda p: list(find_nearest_peaks(p)),
                          analysis_params)
        # Output the results
        for result in results:
            for gene,nearest_peaks,params in result:
                reporter.write_nearest_features(
                    gene,nearest_peaks,
                    peak_file=params.peaks.source_file,
                    feature_file=params.genes.source_file,
                    cutoff=params.cutoff)
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
            xlsx.write_peak_centric(peak_fields)
            xlsx.add_result_sheet('Peak-centric',outputs.peak_centric_out)
            if options.summary:
                xlsx.append_to_notes("\n'Peak-centric (summary)' lists the "
                                     "'top' result (i.e. closest peak/%s "
                                     "pair) for each peak" %
                                     options.feature_type)
                xlsx.add_result_sheet('Peak-centric (summary)',
                                      outputs.peak_centric_summary)
        # Add peaks to features
        if gene_centric:
            xlsx.write_feature_centric(gene_fields)
            xlsx.add_result_sheet('%s-centric' % options.feature_type.title(),
                                  outputs.gene_centric_out)
            if options.summary:
                xlsx.append_to_notes("\n'%s-centric (summary)' lists the "
                                     "'top' result (i.e. closest %s/peak "
                                     "pair) for each %s" %
                                     (options.feature_type.title(),
                                      options.feature_type,
                                      options.feature_type))
                xlsx.add_result_sheet('%s-centric (summary)' %
                                      options.feature_type.title(),
                                      outputs.gene_centric_summary)
        xlsx.write()
        print "Wrote %s" % outputs.xlsx_out
        print

    # Finished
    print "Done"
