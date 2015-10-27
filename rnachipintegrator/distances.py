#!/bin/env python
#
#     distances.py: functions for determining distances and overlaps
#     Copyright (C) University of Manchester 2011-15 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
distances.py

Functions for determining distances and overlaps:

- closestDistanceToRegion:    distance from a reference position to a
                              coordinate or region
- regions_overlap:            determine if two regions have any overlap
- GetNearestTranscriptToPeak: compare two transcripts against the same
                              region
- direction:                  indicate whether two regions are upstream,
                              downstream, or overlapping

Also defines constants:

- UPSTREAM
- DOWNSTREAM
- OVERLAP

"""
# Module constants (used for 'direction' function)
UPSTREAM = -1
DOWNSTREAM = 1
OVERLAP = 0

def closestDistanceToRegion(reference,position1,position2=None,
                            zero_inside_region=False):
    """Return distance from a reference position to a coordinate or region

    For a single specified position, return the absolute distance between
    that position and the reference point.
        
    If a second position is given (specifying a region) then return the
    smallest absolute distance of (reference,position1) and
    (reference,position2).
        
    By default there is no special treatment when the reference lies inside
    the region specified by two positions; to return zero distance in
    these cases, set the 'zero_inside_region' argument to True.
    """
    # Only one position given
    if position2 is None:
        return abs(reference - position1)
    # Is point inside specified region?
    if zero_inside_region:
        if position1 < position2:
            upper,lower = position2,position1
        else:
            upper,lower = position1,position2
        if (lower < reference and reference < upper):
            return 0
    # Outside specified region - return absolute distance
    return min(abs(reference - position1),abs(reference - position2))

def regions_overlap(region1,region2):
    """Determine if two regions have any overlap

    Returns True if any part of the first region overlaps any part of
    the second region (including one region being completely inside the
    other).

    Arguments:
      region1: tuple of numbers representing the limits of the first
        region, can be in any order e.g. (0,10) or (10,0)
      region2: tuple of numbers representing the limits of the second
        region.

    Returns:
      True if there is overlap, False otherwise.
    """
    # Widest region
    if (abs(region1[0] - region1[1]) > abs(region2[0] - region2[1])):
        wide,narrow = region1,region2
    else:
        wide,narrow = region2,region1
    # Determine upper/lower limits of region #1
    if wide[0] < wide[1]:
        lower,upper = wide[0],wide[1]
    else:
        lower,upper = wide[1],wide[0]
    # Regions overlap if either start or end of region #2 lies
    # within region #1 (or vice versa)
    return ((lower <= narrow[0] and narrow[0] <= upper) or
            (lower <= narrow[1] and narrow[1] <= upper))

def GetNearestTranscriptToPeak(rna_data1,rna_data2,chip_peak):
    """Compare two RNA-seq transcripts against the same ChIP-seq peak.

    Given two RNASeqDataLine objects (one for each of the transcripts)
    and a ChIPSeqDataLine object (for the ChIP peak), return the
    transcript which lies closest (regardless of direction) to the ChIP
    peak.

    Note: if the two transcripts are equidistant from the ChIP peak
    then return the first one.
    """

    # Check that one or both are defined
    if not rna_data1:
        return rna_data2
    elif not rna_data2:
        return rna_data1

    # Get the distances from the peak
    distance1 = rna_data1.getTSS() - chip_peak.start
    distance2 = rna_data2.getTSS() - chip_peak.start

    # See which is closest
    dist_diff = abs(distance2) - abs(distance1)
    if dist_diff < 0:
        # Transcript 2 is nearest
        return rna_data2
    else:
        # Transcript 1 is nearest
        return rna_data1

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
    return min(closestDistanceToRegion(peak.start,
                                       feature.tss,feature.tes,
                                       zero_inside_region=True),
               closestDistanceToRegion(peak.end,
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
    return closestDistanceToRegion(feature.tss,
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
    return closestDistanceToRegion(feature.tes,
                                   peak.start,peak.end,
                                   zero_inside_region=True)

def direction(region1,region2,strand=None):
    """
    Direction (up- or downstream) of one region relative to another

    Given two regions (either Peak or Feature instances),
    determine the direction of the second region relative to
    the first.

    For positive strand, 'upstream' means that the test region
    occurs "before" the reference (i.e. has smaller coordinate
    positions), 'downstream' indicates it occurs after. For the
    negative strand the meaning is reversed.

    If region1 partially overlaps region2 then the upstream or
    downstream position will still be determined. However if
    one region wholly overlaps the other then 'overlap' will
    be returned instead.

    Arguments:
      region1 (Peak/Feature): reference region
      region2 (Peak/Feature): region to get direction of
        relative to region1
      strand (str): (optional) explicit strand direction
        ('+' or '-').

    Returns:
      int: UPSTREAM, DOWNSTREAM or OVERLAP (module constants)
        indicating direction of test region relative to
        reference region.

    """
    # Get strandedness info
    if strand is None:
        try:
            strand = region1.strand
        except AttributeError:
            try:
                strand = region2.strand
            except AttributeError:
                strand = '+'
    # Check for overlapping regions
    if (region1.start >= region2.start and region1.end <= region2.end) or \
       (region2.start >= region1.start and region2.end <= region1.end):
        return OVERLAP
    if region2.start >= region1.start and region2.start <= region1.end:
        # region2 start inside feature:
        #
        # |              2222222     |
        # |--------111111111---------|
        #
        # Return distances from peak end to feature end
        return (DOWNSTREAM if strand == '+' else UPSTREAM)
    if region2.end >= region1.start and region2.end <= region1.end:
        # region2 end inside region1
        #
        # |    2222222               |
        # |--------111111111---------|
        #
        return (UPSTREAM if strand == '+' else DOWNSTREAM)
    if region1.end < region2.start:
        # No overlap: peak start to feature end is smallest
        #
        # |                2222222    |
        # |----111111111--------------|
        #
        return (DOWNSTREAM if strand == '+' else UPSTREAM)
    else:
        # No overlap: region2 is downstream of region1
        #
        # |    2222222                |
        # |-------------111111111-----|
        #
        return (UPSTREAM if strand == '+' else DOWNSTREAM)
