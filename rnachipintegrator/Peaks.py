#!/bin/env python
#
#     Peaks.py: classes for handling peak data
#     Copyright (C) University of Manchester 2011-15 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
Peaks.py

Classes for handling peak data

"""

import logging
from utils import make_errline

class PeakSet:
    """Class for storing a set of peaks

    ChIP-seq data consists of ChIP peak information, each of which are
    stored individually in Peak objects. This class is a container for
    a collection of Peak objects and provides methods to operate on the
    collection, by creating subsets by filtering, and sorting the data
    based on various criteria.

    """
    def __init__(self,peaks_file=None,peaks_list=None):
        """Create a new PeakSet instance
        
        Arguments:
          peaks_file (str): (optional) the name of an input file
            to read the peak data from
          peaks_list (list): (optional) list of Peak objects to
            populate the PeakSet with

        """
        self.peaks = []
        if peaks_file:
            self.loadPeaksFromFile(peaks_file)
        elif peaks_list:
            for peak in peaks_list:
                self.addPeak(peak)

    def loadPeaksFromFile(self,peaks_file):
        """Read peaks data from a file and populate the object

        Arguments:
          peaks_file: the name of the input file to read peaks
          data from.

        """
        fp = open(peaks_file,'rU')
        for line in fp:
            # Skip lines that start with a # symbol
            if line.startswith('#'):
                logging.debug("Peaks file: skipped line: %s" % line.strip())
                continue
            # Lines are tab-delimited and have at least 3 columns:
            # chr  start  end
            items = line.strip().split('\t')
            if len(items) < 3:
                logging.warning("Peaks file: skipped line: %s" % line.strip())
                logging.warning("Insufficient number of fields (%d)" % \
                                    len(items))
                continue
            # Check that items in 2nd and 3rd columns are digits
            if not items[1].isdigit() or not items[2].isdigit():
                logging.warning("Peaks file: skipped line: %s" % line.strip())
                # Indicate problem field(s)
                errline = []
                for i in range(len(items)):
                    if i == 1 or i == 2:
                        if not items[i].isdigit():
                            errline.append("^"*len(items[i]))
                        else:
                            errline.append(" "*len(items[i]))
                    else:
                        errline.append(" "*len(items[i]))
                logging.warning("                         %s" % \
                                    ('\t'.join(errline)))
                continue
            # Store in a new Peak object
            peak = Peak(items[0],
                        items[1],
                        items[2])
            self.peaks.append(peak)
        fp.close()
        # Return a reference to this object
        return self

    def addPeak(self,peak):
        """Append a Peak to the PeakSet object

        Arguments:
          peak: a Peak instance.

        """
        self.peaks.append(peak)

    def isSummit(self):
        """Check whether peak set consists of summits only

        Checks the difference between start and end positions for
        stored Peak objects - if all differences are equal to 1 then
        returns True, indicating that the peaks are described as
        summits; otherwise returns False.

        """
        for peak in self.peaks:
            if (peak.end - peak.start) > 1:
                return False
        return True

    def filterByChr(self,matchChr):
        """Return a subset of data filtered by specified chromosome name

        Returns a new PeakSet object containing only the data from
        the current object which matches the specified criteria.

        """
        # Make a new (empty) PeakSet object
        peaks_subset = PeakSet()
        # Populate with only the matching data lines
        for peak in self.peaks:
            if peak.chrom == matchChr:
                peaks_subset.addPeak(peak)
        return peaks_subset

    def filterByPosition(self,limit1,limit2,exclude_limits=False):
        """Return a subset of peaks filtered by position

        Returns a new PeakSet object containing only the data from
        the current object where the start positions fall within a
        region defined by upper and lower limits.

        limits can be supplied in either order (i.e. highest/lowest
        or lowest/highest).

        If exclude_limits is False (the default) then start positions
        that fall exactly on one of the boundaries are counted as
        being within the region; if it is True then these start
        positions will not be considered to lie inside the region.

        """
        # Sort out upper and lower limits
        if limit1 > limit2:
            upper,lower = limit1,limit2
        else:
            upper,lower = limit2,limit1
        # Make a new (empty) PeakSet object
        peaks_subset = PeakSet()
        # Populate with only the matching data lines
        for peak in self.peaks:
            position = peak.start
            if exclude_limits:
                if lower < position and position < upper:
                    peaks_subset.addPeak(peak)
            else:
                if lower <= position and position <= upper:
                    peaks_subset.addPeak(peak)
        return peaks_subset

    def sortByDistanceFrom(self,position):
        """Sort the peaks into order based on distance from a position
    
        Sorts the Peak objects into order of absolute distance of
        their start to the specified position (closest first).
        
        Note that this operates on the current object.

        """
        self.peaks = sorted(self.peaks,
                            key=lambda record: min(abs(record.start - position),
                                                   abs(record.end - position)))
        return self

    def __iter__(self):
        return iter(self.peaks)

    def __getitem__(self,key):
        try:
            start = key.start
            stop = key.stop
            step = key.step
            slice_ = PeakSet()
            for peak in self.peaks[start:stop:step]:
                slice_.addPeak(peak)
            return slice_
        except AttributeError:
            return self.peaks[key]

    def __len__(self):
        return len(self.peaks)

    def __eq__(self,other):
        if len(self) != len(other):
            return False
        for p1,p2 in zip(self,other):
            if p1 != p2:
                return False
        return True

    def __ne__(self,other):
        if len(self) != len(other):
            return True
        for p1,p2 in zip(self,other):
            if p1 != p2:
                return True
        return False

class Peak:
    """Class for storing a peak

    Access the data from the line using the object's properties:

      chrom
      start
      end

    There are also convenience methods (e.g. insideRegion).

    """
    def __init__(self,chrom,start,end):
        self.chrom = chrom.strip('"\'')
        self.start = int(start)
        self.end = int(end)

    def __repr__(self):
        return "%s\t%s\t%s" % (self.chrom,
                               self.start,
                               self.end)

    def __eq__(self,other):
        return \
            (self.chrom == other.chrom) and \
            (self.start == other.start) and \
            (self.end == other.end)

    def __ne__(self,other):
        return \
            (self.chrom != other.chrom) or \
            (self.start != other.start) or \
            (self.end != other.end)

    def insideRegion(self,limit1,limit2,exclude_limits=False):
        """Check if the peak is contained within a defined region

        The region is defined by upper and lower limits: if the start
        position lies within these limits then return True, otherwise
        return false.

        limits can be supplied in either order (i.e. highest/lowest
        or lowest/highest).

        By default if the start position lies exactly on one of the
        limits then it is also counted as being within the region;
        if exclude_limits is set to True then the limits are not
        considered part of the region.
        
        """
        if limit1 > limit2:
            upper,lower = limit1,limit2
        else:
            upper,lower = limit2,limit1
        if not exclude_limits:
            return (lower <= self.start and self.start <= upper)
        else:
            return (lower < self.start and self.start < upper)
