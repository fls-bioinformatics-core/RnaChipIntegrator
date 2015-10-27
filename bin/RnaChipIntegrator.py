#!/usr/bin/env python
#
#     RnaChipIntegrator.py: analyse genomic features with peak data
#     Copyright (C) University of Manchester 2011-15 Peter Briggs, Leo Zeef
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

Analyse genomic features with peak data.

"""

#######################################################################
# Imports
#######################################################################

import sys
import os
import optparse
import time,datetime
from rnachipintegrator.Features import FeatureSet
from rnachipintegrator.Peaks import PeakSet
import rnachipintegrator.analysis as analysis
import rnachipintegrator.output as output
import rnachipintegrator.xls_output as xls_output
import logging
logging.getLogger().setLevel(logging.WARNING)
logging.basicConfig(format='%(levelname)s: %(message)s')

import rnachipintegrator
__version__ = rnachipintegrator.get_version()

#######################################################################
# Main program
#######################################################################

if __name__ == '__main__':
    
    # Defaults
    promoter_region = (1000,100)
    max_distance = 1000000
    max_closest = 4
    
    # Parse command line
    p = optparse.OptionParser(usage="%prog [options] FEATURES PEAKS",
                              version="%prog "+__version__,
                              description=
                              "Analyse FEATURES (any set of genomic features) "
                              "against PEAKS (a set of regions) and report "
                              "nearest features to each peak (and vice versa)")
    p.add_option('--cutoff',action='store',dest='max_distance',
                 type='int',default=max_distance,
                 help="Maximum distance allowed between peaks and "
                 "features before being omitted from the analysis "
                 "(default %dbp; set to zero for no cutoff)" %
                 max_distance)
    p.add_option('--number',action='store',dest='max_closest',
                 type='int',default=max_closest,
                 help="Maximum number of hits to report from the analyses "
                 "(default %d; set to zero to report all hits)" %
                 max_closest)
    p.add_option('--edge',action='store',dest="edge",type="choice",
                 choices=('tss','both'),default='tss',
                 help="Feature edges to consider when calculating distances "
                 "between features and peaks, either: 'tss' (default: only "
                 "use TSS) or 'both' (use whichever of TSS or TES gives "
                 "shortest distance)")
    p.add_option('--promoter_region',action="store",dest="promoter_region",
                 default="%d,%d" % promoter_region,
                 help="Define promoter region with respect to feature TSS "
                 "in the form UPSTREAM,DOWNSTREAM (default -%d to %dbp of "
                 "TSS)" %  promoter_region)
    p.add_option('--only-DE',action='store_true',
                 dest='only_diff_expressed',default=False,
                 help="Only use features flagged as differentially expressed "
                 "in analyses (input feature data must include DE flags)")
    p.add_option('--name',action='store',dest='name',default=None,
                 help="Set basename for output files")
    p.add_option('--compact',action='store_true',dest='compact',default=False,
                 help="Output minimal information in a compact format")
    p.add_option('--summary',action='store_true',dest='summary',default=False,
                 help="Output 'summary' for each analysis, consisting of "
                 "only the top hit for each peak or feature (cannot be used "
                 "with --compact)")
    p.add_option('--pad',action="store_true",dest="pad_output",
                 help="Where less than MAX_CLOSEST hits are found, pad "
                 "output with blanks to ensure that MAX_CLOSEST hits "
                 "are still reported")
    p.add_option('--xls',action="store_true",dest="xls_output",
                 help="Output XLS spreadsheet with results")
    p.add_option('--feature',action="store",dest="feature_type",
                 help="rename features to FEATURE_TYPE in output (e.g. "
                 "'gene', 'transcript' etc)")
    p.add_option('--peak_cols',action="store",dest="peak_cols",
                 help="list of 3 column indices (e.g. '1,4,5') indicating "
                 "columns to use for chromosome, start and end from the "
                 "input peak file (if not first three columns).")
    options,args = p.parse_args()

    # Input files
    if len(args) != 2:
        p.error("need to supply 2 files (features and peaks)")
    feature_file,peak_file = args

    # Report version and authors
    p.print_version()
    print
    print "Find nearest peaks to genomic features (and vice versa)"
    print
    print "University of Manchester"
    print "Faculty of Life Sciences"
    print "Bioinformatics Core Facility"
    print "Authors: Ian Donaldson, Leo Zeef and Peter Briggs"
    print

    # Promoter region
    promoter = (abs(int(options.promoter_region.split(',')[0])),
                abs(int(options.promoter_region.split(',')[1])))

    # Reporting options
    max_distance = options.max_distance
    if max_distance <= 0:
        max_distance = None
    max_closest = options.max_closest
    if max_closest <= 0:
        max_closest = None

    # Feature edge to use
    if options.edge == 'tss':
        tss_only = True
    else:
        tss_only = False

    # Reporting formats
    if options.compact:
        mode = output.SINGLE_LINE
        peak_fields = ('chr','start','end','list(feature.id)')
        feature_fields = ('feature.id','list(chr,start,end,dist_closest)')
        placeholder = '.'
        if options.summary:
            options.summary = False
            logging.error("--summary option not compatible with --compact")
            sys.exit(1)
    else:
        mode = output.MULTI_LINE
        peak_fields = ('chr','start','end',
                       'feature.id','strand','TSS','TES',
                       'dist_closest','dist_TSS','dist_TES',
                       'overlap_feature',
                       'overlap_promoter')
        feature_fields = ('feature.id',
                          'feature.chr','feature.start','feature.end',
                          'feature.strand',
                          'peak.chr','peak.start','peak.end','order',
                          'dist_closest','dist_TSS','dist_TES')
        placeholder = '---'

    # Feature type
    if options.feature_type is None:
        feature_type = 'feature'
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
    print "Input features file: %s" % feature_file
    print "Input peaks file   : %s" % peak_file
    print
    print "Maximum cutoff distance: %d (bp)" % max_distance
    print "Maximum no. of hits    : %d" % max_closest
    print "Promoter region        : -%d to %d (bp from TSS)" % \
        promoter_region
    print
    if tss_only:
        print "Distances will be calculated from feature TSS only"
    else:
        print "Distances will be calculated from nearest of feature TSS or TES"
    print
    print "Features will be referred to as '%ss'" % feature_type
    print

    # Read in feature data
    try:
        features = FeatureSet(feature_file)
    except Exception, ex:
        logging.critical("Failed to read in feature data: %s" % ex)
        print "Please fix errors in input file before running again"
        sys.exit(1)
    if not len(features):
        logging.error("No feature data read in")
        sys.exit(1)
    print "%d feature records read in" % len(features)

    # Differential expression handling
    use_differentially_expressed = False
    if features.isFlagged():
        print "\tFeature data include differential expression flag"
        print "\t%d features flagged as differentially expressed" % \
            len(features.filterByFlag(1))
        if options.only_diff_expressed:
            print
            print "*** Only differentially expressed features will used ***"
            use_differentially_expressed = True
    elif options.only_diff_expressed:
        logging.error("--only-DE flag needs input features flagged as "
                      "differentially expressed")
        sys.exit(1)
    print

    # Read in peak data
    print "Using columns %s from peaks file as chrom, start, end" % \
        (peak_cols,)
    peaks = PeakSet(peak_file,columns=peak_cols)
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
        basename = os.path.splitext(os.path.basename(feature_file))[0]

    # Do the analyses
    print "**** Nearest features to peaks ****"
    outfile = basename+"_features_per_peak.txt"
    if options.summary:
        summary = basename+"_features_per_peak_summary.txt"
    else:
        summary = None
    reporter = output.AnalysisReportWriter(mode,peak_fields,
                                           promoter_region=promoter,
                                           null_placeholder=placeholder,
                                           max_hits=max_closest,
                                           pad=options.pad_output,
                                           outfile=outfile,
                                           summary=summary,
                                           feature_type=options.feature_type)
    for peak,nearest_features in analysis.find_nearest_features(
            peaks,features,tss_only=tss_only,distance=max_distance,
            only_differentially_expressed=use_differentially_expressed):
        reporter.write_nearest_features(peak,nearest_features)
    reporter.close()
    print "Results written to %s" % outfile
    if summary:
        print "Summary written to %s" % summary
    print

    print "**** Nearest peaks to features ****"
    outfile = basename+"_peaks_per_feature.txt"
    if options.summary:
        summary = basename+"_peaks_per_feature_summary.txt"
    else:
        summary = None
    reporter = output.AnalysisReportWriter(mode,feature_fields,
                                           null_placeholder=placeholder,
                                           max_hits=max_closest,
                                           pad=options.pad_output,
                                           outfile=outfile,
                                           summary=summary,
                                           feature_type=options.feature_type)
    for feature,nearest_peaks in analysis.find_nearest_peaks(
            features,peaks,tss_only=tss_only,distance=max_distance,
            only_differentially_expressed=use_differentially_expressed):
        reporter.write_nearest_peaks(feature,nearest_peaks)
    reporter.close()
    print "Results written to %s" % outfile
    if summary:
        print "Summary written to %s" % summary
    print

    # Make XLS file
    if options.xls_output:
        print "**** Writing XLS file ****"
        xls = xls_output.XLS(p.get_version(),feature_type)
        # Write the settings
        xls.append_to_notes("Input %ss file\t%s" % (feature_type,
                                                    feature_file))
        xls.append_to_notes("Input peaks file\t%s" % peak_file)
        xls.append_to_notes("Maximum cutoff distance (bp)\t%d" % max_distance)
        xls.append_to_notes("Maximum no. of hits to report\t%d" % max_closest)
        xls.append_to_notes("Promoter region (bp from TSS)\t-%d to %d" %
                            promoter_region)
        if tss_only:
            xls.append_to_notes("Distances calculated from\tTSS only")
        else:
            xls.append_to_notes("Distances calculated from\tTSS or TES")
        xls.append_to_notes("Only use differentially expressed %ss\t%s" %
                            (feature_type,
                             "Yes" if use_differentially_expressed else "No"))
        # Add features to peaks
        xls.write_features_to_peaks(peak_fields)
        xls.add_result_sheet('%ss' % feature_type.title(),
                             basename+"_features_per_peak.txt")
        if options.summary:
            xls.add_result_sheet('%ss (summary)' % feature_type.title(),
                                 basename+"_features_per_peak_summary.txt")
        # Add peaks to features
        xls.write_peaks_to_features(feature_fields)
        xls.add_result_sheet('Peaks',basename+"_peaks_per_feature.txt")
        if options.summary:
            xls.add_result_sheet('Peaks (summary)',
                                 basename+"_peaks_per_feature_summary.txt")
        xls.write(basename+'.xls')
        print "Wrote %s" % basename+'.xls'
        print

    # Finished
    print "Done"
