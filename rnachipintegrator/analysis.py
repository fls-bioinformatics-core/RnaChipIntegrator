#!/bin/env python
#
#     analysis.py: analyses of peaks vs features and vice versa
#     Copyright (C) University of Manchester 2011-15 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
analysis.py

Functions for analysing peaks against features, and vice versa:

- find_nearest_features
- find_nearest_peaks

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
    # Find nearest features to each peak
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
                    if distances.distance_tss(peak,feature) > distance:
                        break
                else:
                    if distances.distance_closest_edge(peak,feature) > distance:
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
                if distances.distance_tss(peak,feature) > distance:
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
                               distances.edge_distances(peak,f))

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
                               distances.tss_distances(peak,f))

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
                         distances.edge_distances(p,feature))

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
                         distances.tss_distances(p,feature))
