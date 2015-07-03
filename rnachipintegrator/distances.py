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

"""

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
