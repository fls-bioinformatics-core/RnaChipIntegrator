#!/bin/env python
#
#     analysis_redux.py: reimplement functions for performing analyses
#     Copyright (C) University of Manchester 2011-15 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
analysis_redux.py

Reimplementation of analysis functions.

"""
import os
import logging
import distances
from Features import FeatureSet
from Peaks import PeakSet

#######################################################################
# Analysis functions
#######################################################################

# NB there are a few potential optimisations, e.g.:
# - caching the results of filtering on chromosomes
# - improving the search over distances e.g. by sorting once on
#   distances from an arbitrary zero, and then calculating distances
#   between things with reference to that?

def find_nearest_features(peaks,features,distance=None,
                          max_closest=None,tss_only=False,
                          only_differentially_expressed=False):
    """
    Locate nearest features for each peak

    Arguments:
      features (FeatureList): list of features
      peaks (PeakList): list of peaks
      distance (int): optional cut-off distance to apply
      max_closest (int): optional maximum number of peaks
        to find per feature
      tss_only (bool): only consider distances from the
        feature TSS (default is to consider distances from
        both the TSS and TES)
      only_differentially_expressed (bool): only consider
        features that are flagged as differentially expressed

    Yields:
      tuple: Peak object and a FeatureSet object with the
        nearest features, for each peak in the input PeakSet

    """
    # Reimplemtation & generalisation of:
    # - AnalyseNearestTSSToSummits
    # - AnalyseNearestTranscriptsToPeakEdges
    for peak in peaks:
        # Only consider features on same chromosome
        feature_list = features.filterByChr(peak.chrom)
        # Differentially-expressed features only?
        if only_differentially_expressed:
            feature_list = feature_list.filterByFlag(1)
        if tss_only:
            feature_list.sortByClosestTSSTo(peak.start,peak.end)
        else:
            feature_list.sortByClosestEdgeTo(peak.start,peak.end)
        # Apply distance cut-off
        if distance is not None:
            closest = FeatureSet()
            for feature in feature_list:
                if tss_only:
                    if distance_tss(peak,feature) > distance:
                        break
                else:
                    if distance_closest_edge(peak,feature) > distance:
                        break
                closest.addFeature(feature)
            feature_list = closest
        # Reduce to maximum number of features
        if max_closest is not None:
            feature_list = feature_list[:max_closest]
        # Return result
        yield (peak,feature_list)

def find_nearest_peaks(features,peaks,distance=None,
                       max_closest=None,
                       only_differentially_expressed=False):
    """
    Locate nearest peaks for each feature

    Arguments:
      features (FeatureList): list of features
      peaks (PeakList): list of peaks
      distance (int): optional cut-off distance to apply
      max_closest (int): optional maximum number of peaks
        to find per feature
      only_differentially_expressed (bool): only consider
        features that are flagged as differentially expressed

    Yields:
      tuple: Feature object and a PeakSet object with the
        nearest peaks, for each feature in the input
        FeatureSet

    """
    # Reimplementation & generalisation of
    # AnalyseNearestPeaksToTranscripts function
    # Reduce to set of differentially expressed features
    if only_differentially_expressed:
        features = features.filterByFlag(1)
    # Find nearest peaks for each feature
    for feature in features:
        # Only consider peaks on same chromosome
        peak_list = peaks.filterByChr(feature.chrom)
        # Sort into distance order
        peak_list.sortByDistanceFrom(feature.tss)
        # Apply distance cut-off
        if distance is not None:
            closest = PeakSet()
            for peak in peak_list:
                #if min(abs(peak.start - feature.tss),
                #       abs(peak.end - feature.tss)) > distance:
                if distance_tss(peak,feature) > distance:
                    break
                closest.addPeak(peak)
            peak_list = closest
        # Reduce to maximum number of peaks
        if max_closest is not None:
            peak_list = peak_list[:max_closest]
        # Return results
        yield (feature,peak_list)

def distance_closest_edge(peak,feature):
    """
    Get distance from a peak to a feature

    Arguments:
      peak (Peak): peak
      feature (Feature): feature

    Returns:
      int: Smallest absolute distance from the peak to
        the feature.

    """
    return min(distances.closestDistanceToRegion(peak.start,
                                                 feature.tss,feature.tes,
                                                 zero_inside_region=True),
               distances.closestDistanceToRegion(peak.end,
                                                 feature.tss,feature.tes,
                                                 zero_inside_region=True))

def distance_tss(peak,feature):
    """
    Get distance from a peak to a feature TSS

    Arguments:
      peak (Peak): peak
      feature (Feature): feature

    Returns:
      int: Smallest absolute distance from the peak to
        the feature TSS.

    """
    return distances.closestDistanceToRegion(feature.tss,
                                             peak.start,peak.end,
                                             zero_inside_region=True)

def distance_tes(peak,feature):
    """
    Get distance from a peak to a feature TES

    Arguments:
      peak (Peak): peak
      feature (Feature): feature

    Returns:
      int: Smallest absolute distance from the peak to
        the feature TES.

    """
    return distances.closestDistanceToRegion(feature.tes,
                                             peak.start,peak.end,
                                             zero_inside_region=True)

