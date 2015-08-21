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

def find_nearest_features(peaks,features,distance=None,tss_only=False,
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
      pad (bool): add extra 'None' items to output
        FeatureSet so that it contains max_closest results

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
            sort_features_by_tss_distances(peak,feature_list)
        else:
            sort_features_by_edge_distances(peak,feature_list)
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
        # Return at least one (null) result
        if not feature_list:
            feature_list.addFeature(None)
        # Return result
        yield (peak,feature_list)

def find_nearest_peaks(features,peaks,distance=None,tss_only=False,
                       only_differentially_expressed=False):
    """
    Locate nearest peaks for each feature

    Arguments:
      features (FeatureList): list of features
      peaks (PeakList): list of peaks
      distance (int): optional cut-off distance to apply
      tss_only (bool): only consider distances from the
        feature TSS (default is to consider distances from
        both the TSS and TES)
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
        if tss_only:
            sort_peaks_by_tss_distances(feature,peak_list)
        else:
            sort_peaks_by_edge_distances(feature,peak_list)
        # Apply distance cut-off
        if distance is not None:
            closest = PeakSet()
            for peak in peak_list:
                if distance_tss(peak,feature) > distance:
                    break
                closest.addPeak(peak)
            peak_list = closest
        # Return at least one (null) result
        if not peak_list:
            peak_list.addPeak(None)
        # Return results
        yield (feature,peak_list)

def sort_features_by_edge_distances(peak,features):
    """
    Sort features by edge-to-edge distances to a peak

    Arguments:
      peak (Peak): peak instance
      features (FeatureSet): set of features that will
        be sorted into order according to the smallest
        distance of their edges from the edges of the
        peak. The sorting is done in place.

    """
    features.features = sorted(features.features,
                               key = lambda f:
                               edge_distances(peak,f))

def sort_features_by_tss_distances(peak,features):
    """
    Sort features by TSS-to-edge distances to a peak

    Arguments:
      peak (Peak): peak instance
      features (FeatureSet): set of features that will
        be sorted into order according to the smallest
        distance of their TSS positions to the edges of
        the peak. The sorting is done in place.

    """
    features.features = sorted(features.features,
                               key = lambda f:
                               tss_distances(peak,f))

def sort_peaks_by_edge_distances(feature,peaks):
    """
    Sort peaks by edge-to-edge distances to a feature

    Arguments:
      feature (Feature): feature instance
      peaks (PeakSet): set of peaks that will be
        sorted into order according to the smallest
        distance of their edges from the edges of the
        feature. The sorting is done in place.

    """
    peaks.peaks = sorted(peaks.peaks,
                         key = lambda p:
                         edge_distances(p,feature))

def sort_peaks_by_tss_distances(feature,peaks):
    """
    Sort peaks by edge-to-TSS distances to a feature

    Arguments:
      feature (Feature): feature instance
      peaks (PeakSet): set of peaks that will be
        sorted into order according to the smallest
        distance of their edges from the TSS position
        of the feature. The sorting is done in place.

    """
    peaks.peaks = sorted(peaks.peaks,
                         key = lambda p:
                         tss_distances(p,feature))

def edge_distances(peak,feature):
    """
    Return distances between peak and feature edges

    Arguments:
      peak (Peak): peak instance
      feature (Feature): feature instance

    Returns:
      tuple: a pair of distances, the first being the
        smallest distance between an edge of the peak
        region and an edge of the feature region and
        the second being the distance from the same peak
        edge to the opposite feature edge.
        Both distances are zero if the peak is entirely
        contained within the feature (or vice versa).
        If there is partial overlap then the smallest
        distance is returned as zero, the second is
        the amount of the peak that lies outside the
        feature.

    """
    if ((peak.start >= feature.start and
         peak.end <= feature.end) or
        (feature.start >= peak.start and
         feature.end <= peak.end)):
        # Peak entirely contained in feature (or vice versa)
        return (0,0)
    # All possible distances
    ds_fs = abs(peak.start - feature.start)
    ds_fe = abs(peak.start - feature.end)
    de_fs = abs(peak.end - feature.start)
    de_fe = abs(peak.end - feature.end)
    if peak.start >= feature.start and peak.start <= feature.end:
        # Peak start inside feature:
        #
        # |              PPPPPPP     |
        # |--------FFFFFFFFF---------|
        #
        # Return distances from peak end to feature end
        return (0,de_fe)
    if peak.end >= feature.start and peak.end <= feature.end:
        # Peak end inside feature
        #
        # |    PPPPPPP               |
        # |--------FFFFFFFFF---------|
        #
        # Return distances from peak start to feature start
        return (0,ds_fs)
    if ds_fe < de_fs:
        # No overlap: peak start to feature end is smallest
        #
        # |                PPPPPPP    |
        # |----FFFFFFFFF--------------|
        #
        return (ds_fe,ds_fs)
    else:
        # No overlap: peak end to feature start is smallest
        #
        # |    PPPPPPP                |
        # |-------------FFFFFFFFF-----|
        #
        return (de_fs,de_fe)

def tss_distances(peak,feature):
    """
    Return distances between peak edges and feature TSS

    Arguments:
      peak (Peak): peak instance
      feature (Feature): feature instance

    Returns:
      tuple: a pair of distances, the first being the
        smallest distance between an edge of the peak
        region to the TSS of the feature region, and
        the second being the distance from the other peak
        edge to the TSS.
        If the TSS position is entirely contained within
        the peak then the smallest distance is returned as
        zero, the second is the average of the distances
        to the two peak edges from the TSS.

    """
    # TSS is considered as a point
    if (feature.tss >= peak.start and
        feature.tss <= peak.end):
        # TSS entirely contained in the peak
        # Rank using the average distance of the peak
        # edges from the TSS
        return (0,(abs(feature.tss - peak.start) +
                   abs(feature.tss - peak.end))/2)
    # Possible distances to TSS
    d_tss = [abs(feature.tss - peak.start),
             abs(feature.tss - peak.end)]
    d_tss.sort()
    return tuple(d_tss)

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

