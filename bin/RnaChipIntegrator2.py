#!/usr/bin/env python
#
#     RnaChipIntegrator.py: analyse RNA-seq and ChIP-seq data
#     Copyright (C) University of Manchester 2011-14 Peter Briggs, Leo Zeef
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
import rnachipintegrator.analysis_redux as analysis
import rnachipintegrator.output as output

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
                              version="%prog "+__version__)
    p.add_option('--cutoff',action='store',dest='max_distance',
                 type='int',default=max_distance,
                 help="Maximum distance allowed between peaks and "
                 "features before being omitted from the analysis "
                 "(default %d bp; set to zero for no cutoff)" %
                 max_distance)
    p.add_option('--number',action='store',dest='max_closest',
                 type='int',default=max_closest,
                 help="Maximum number of hits to report from the analyses "
                 "(default %d; set to zero to report all hits)" %
                 max_closest)
    p.add_option('--promoter_region',action="store",dest="promoter_region",
                 default="%d,%d" % promoter_region,
                 help="Define promoter region with respect to feature TSS "
                 "in the form UPSTREAM,DOWNSTREAM (default -%d to %d bp of "
                 "TSS)" %  promoter_region)
    p.add_option('--compact',action='store_true',dest='compact',default=False,
                 help="Output minimal information in a compact format")
    p.add_option('--no-DE',action='store_true',
                 dest='no_differential_expression',default=False,
                 help="Ignore differential expression flags (even if "
                 "present in input)")
    p.add_option('--pad',action="store_true",dest="pad_output",
                 help="Where less than MAX_CLOSEST hits are found, pad "
                 "output with blanks to ensure that MAX_CLOSEST hits "
                 "are still reported")
    options,args = p.parse_args()

    # Input files
    if len(args) != 2:
        p.error("need to supply 2 files (features and peaks)")
    feature_file,peak_file = args

    # Report version and authors
    p.print_version()
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
    if options.pad_output:
        raise NotImplementedError("--pad not implemented")

    # Reporting formats
    if options.compact:
        mode = output.SINGLE_LINE
        peak_fields = ('chr','start','end','feature_list(feature.id)')
        feature_fields = ('feature.id','peak_list(chr,start,end,dist_closest)')
    else:
        mode = output.MULTI_LINE
        peak_fields = ('chr','start','end',
                       'feature.id',
                       'order','TSS','TES',
                       'dist_closest',
                       'dist_TSS','dist_TES',
                       'strand',
                       'in_the_gene',
                       'overlap_feature',
                       'overlap_promoter')
        feature_fields = ('feature.id','chr','start','end','order')

    # Report settings
    print "Input features file: %s" % feature_file
    print "Input peaks file   : %s" % peak_file
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
    use_differential_expression = False
    if features.isFlagged():
        print "\tFeature data include differential expression flag"
        print "\t%d features flagged as differentially expressed" % \
            len(features.filterByFlag(1))
        if not options.no_differential_expression:
            print
            print "*** Only differentially expressed features will used"
            print "*** Rerun with --no-DE to analyse all features "
            use_differential_expression = True
    print

    # Read in peak data
    peaks = PeakSet(peak_file)
    if not len(peaks):
        logging.error("No peak data read in")
        sys.exit(1)
    print "%d peak records read in" % len(peaks)
    if peaks.isSummit():
        print "\tPeak data are summits"
    else:
        print "\tPeak data are regions"
    print

    # Do the analyses
    print "**** Nearest features to peaks (TSS only) ****"
    reporter = output.AnalysisReporter(mode,peak_fields,
                                       promoter_region=promoter)
    for peak,nearest_features in analysis.find_nearest_features(
            peaks,features,tss_only=True,
            distance=max_distance,
            max_closest=max_closest,
            only_differentially_expressed=use_differential_expression):
        for line in reporter.report_nearest_features(peak,nearest_features):
            print line

    print "**** Nearest features to peaks (TSS/TES) ****"
    reporter = output.AnalysisReporter(mode,peak_fields,
                                       promoter_region=promoter)
    for peak,nearest_features in analysis.find_nearest_features(
            peaks,features,tss_only=False,
            distance=max_distance,
            max_closest=max_closest,
            only_differentially_expressed=use_differential_expression):
        for line in reporter.report_nearest_features(peak,nearest_features):
            print line

    print "**** Nearest peaks to features ****"
    reporter = output.AnalysisReporter(mode,feature_fields)
    for feature,nearest_peaks in analysis.find_nearest_peaks(
            features,peaks,
            distance=max_distance,
            max_closest=max_closest,
            only_differentially_expressed=use_differential_expression):
        for line in reporter.report_nearest_peaks(feature,nearest_peaks):
            print line
