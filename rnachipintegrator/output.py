#!/bin/env python
#
#     output.py: functions for outputing analysis results
#     Copyright (C) University of Manchester 2015 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
output.py

Functions for outputing analysis results

"""
import distances
from analysis_redux import distance_closest_edge
from analysis_redux import distance_tss,distance_tes

#######################################################################
# Constants
#######################################################################

MULTI_LINE=0
SINGLE_LINE=1

#######################################################################
# Classes
#######################################################################

class AnalysisReporter:
    """
    Class to handle reporting of analysis results

    Once initialised the reporter can be used to generate 'reports'
    of each peak along with the nearest features (using the
    'report_nearest_features' method) or for each feature along
    with the nearest peaks (using 'report_nearest_peaks').

    Output can be in either 'multi-line' (one line per result pair),
    or 'single-line' format (one line containing all results).

    For each method a list of the fields to be reported can be
    specified. Available fields are:

    - (peak.)chr: chromosome for the peak
    - (peak.)start: peak start position
    - (peak.)end: peak end position
    - (feature.)id: feature ID
    - feature.chr: chromosome for the feature
    - feature.start: feature start position
    - feature.end: feature end position
    - (feature.)TSS: feature TSS
    - (feature.)TES: feature TES
    - (feature.)strand: feature strand
    - dist_closest: closest distance between peak and feature
    - dist_TSS: distance between peak and feature TSS
    - dist_TES: distance between peak and feature TES
    - overlap_feature: 'YES' if peak overlaps the feaure, 'NO' if not
    - overlap_promoter: 1 if peak overlaps the promoter region, 0 if not
    - in_the_gene: synonym for 'overlap_feature'
    - 'differentially_expressed': flag value for feature

    (In the field names above, the parts in (...) are optional e.g.
    'chr' == 'peak.chr' etc.)

    For multi-line output these additional fields are available:

    - order: the 'order' of the feature/peak pair (e.g. '1 of 4')

    For single-line output these additional fields are available:

    - number_of_peaks
    - peak_list(...): output all peaks
    - number_of_features
    - feature_list(...): output all features

    For the '*_list' options, the paranthese should enclose a list
    of fields to output for each peak or feature in the list e.g.
    'peak_list(chr,start,dist_closest)' or 'feature_list(feature.id)'.

    The following fields have not been implemented:

    - 
    - features_inbetween

    """
    def __init__(self,mode,fields,promoter_region=None):
        """
        Create new AnalysisReporter instance

        Arguments:
          promoter_region (tuple): promoter region extent (optional)
          mode (int): either SINGLE_LINE or MULTI_LINE
          fields (list): list of fields to output

        """
        self._fields = fields
        self._mode = mode
        self._promoter_region = promoter_region
        self._context_peak = None
        self._context_feature = None

    def report_nearest_features(self,peak,features):
        """
        Return details of nearest features for a peak
        
        Arguments:
          peak (Peak): peak of interest
          features (FeatureSet): list of nearest features

        Yields:
          string: line(s) of text reporting the results

        """
        # Initialise and set the context
        self._context_peak = peak
        nresults = len(features)
        if self._mode == SINGLE_LINE:
            # Report everything on a single line
            line = []
            for field in self._fields:
                if field == 'number_of_features':
                    value = nresults
                elif field.startswith('feature_list('):
                    # Extract the subfields
                    subfields = field[:-1].split('(')[1].split(',')
                    # Report list of features
                    value = []
                    for feature in features:
                        self._context_feature = feature
                        for subfield in subfields:
                            value.append(self.value_for(subfield))
                    value = '\t'.join([str(x) for x in value])
                else:
                    # All other fields
                    value = self.value_for(field)
                line.append(str(value))
            # Return (yield) the line
            yield '\t'.join(line)
        elif self._mode == MULTI_LINE:
            # Report each result pair on a new line
            i = 0
            for feature in features:
                self._context_feature = feature
                i += 1
                line = []
                for field in self._fields:
                    if field == 'order':
                        value = '%d of %d' % (i,nresults)
                    else:
                        value = self.value_for(field)
                    line.append(str(value))
                # Return (yield) the line
                yield '\t'.join(line)
        # Reset the context
        self._context_peak = None
        self._context_feature = None

    def report_nearest_peaks(self,feature,peaks):
        """
        Return details of nearest peaks for a feature
        
        Arguments:
          feature (Feature): feature of interest
          peaks (PeakSet): list of nearest peaks
          mode (int): either SINGLE_LINE or MULTI_LINE
          fields (list): list of fields to output

        Returns:
          string: block of text reporting the results

        """
        # Initialise and set the context
        self._context_feature = feature
        nresults = len(peaks)
        if self._mode == SINGLE_LINE:
            # Report everything on a single line
            line = []
            for field in self._fields:
                if field == 'number_of_peaks':
                    value = nresults
                elif field.startswith('peak_list('):
                    # Extract the subfields
                    subfields = field[:-1].split('(')[1].split(',')
                    # Report peaks
                    value = []
                    for peak in peaks:
                        self._context_peak = peak
                        for subfield in subfields:
                            value.append(self.value_for(subfield))
                    value = '\t'.join([str(x) for x in value])
                else:
                    # All other fields
                    value = self.value_for(field)
                line.append(str(value))
            # Return (yield) the line
            yield '\t'.join(line)
        elif self._mode == MULTI_LINE:
            # Report each result pair on a new line
            i = 0
            for peak in peaks:
                line = []
                i += 1
                self._context_peak = peak
                # Assemble line from fields
                for field in self._fields:
                    if field == 'order':
                        value = '%d of %d' % (i,nresults)
                    else:
                        value = self.value_for(field)
                    line.append(str(value))
                # Return (yield) the line
                yield '\t'.join(line)
        # Reset the context
        self._context_peak = None
        self._context_feature = None

    def value_for(self,attr):
        """
        Return the value for the specified attribute

        Given the name of a field/attribute (see above for
        a list and definition of each), return the value
        for the current peak/feature pair (which should have
        been set by the calling method in the '_context_peak'
        and '_context_feature' properties).

        Arguments:
          attr (string): attribute name

        Returns:
          Value of the field for the current peak/feature
          pair
        
        """
        peak = self._context_peak
        feature = self._context_feature
        if attr == 'chr' or attr == 'peak.chr':
            return peak.chrom
        elif attr == 'peak.start' or attr == 'start':
            return peak.start
        elif attr == 'peak.end' or attr == 'end':
            return peak.end
        elif attr == 'id' or attr == 'feature.id':
            return feature.id
        elif attr == 'feature.chr':
            return feature.chrom
        elif attr == 'feature.start':
            return feature.start
        elif attr == 'feature.end':
            return feature.end
        elif attr == 'TSS':
            return feature.tss
        elif attr == 'TES':
            return feature.tes
        elif attr == 'strand' or attr == 'feature.strand':
            return feature.strand
        elif attr == 'differentially_expressed':
            return feature.flag
        elif attr == 'dist_closest':
            return distance_closest_edge(peak,feature)
        elif attr == 'dist_TSS':
            return distance_tss(peak,feature)
        elif attr == 'dist_TES':
            return distance_tes(peak,feature)
        elif attr == 'overlap_feature' or attr == 'in_the_gene':
            if distances.regions_overlap((peak.start,peak.end),
                                         (feature.tss,feature.tes)):
                overlap_feature = 1
            else:
                overlap_feature = 0
            if attr == 'in_the_gene':
                overlap_feature = ('YES' if overlap_feature == 1 else 'NO')
            return overlap_feature
        elif attr == 'overlap_promoter':
            if self._promoter_region is not None:
                promoter = feature.getPromoterRegion(*self._promoter_region)
                if distances.regions_overlap((peak.start,peak.end),
                                             promoter):
                    overlap_promoter = 1
                else:
                    overlap_promoter = 0
            else:
                raise Exception("'overlap_promoter' requested but no "
                                "promoter region has been defined")
            return overlap_promoter
        elif attr == 'features_inbetween':
            raise NotImplementedError("'features_inbetween' not implemented")
        else:
            raise KeyError("Unrecognised report field: '%s'" % attr)
