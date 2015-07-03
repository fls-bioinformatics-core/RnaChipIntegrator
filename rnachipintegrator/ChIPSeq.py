#!/bin/env python
#
#     ChIPSeq.py: classes for handling ChIP-seq peak data
#     Copyright (C) University of Manchester 2011-15 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
ChIPSeq.py

Classes for handling ChIP-seq peak data

"""

import logging

class ChIPSeqData:
    """Class for storing ChIP-seq data

    ChIP-seq data consists of ChIP peak information, each of which are
    stored individually in ChIPSeqDataLine objects. This class is a
    container for a collection of ChIPSeqDataLine objects and provides
    methods to operate on the collection, by creating subsets by
    filtering, and sorting the data based on various criteria.
    """
    def __init__(self,chipseq_file=None):
        """Create a new ChIPSeqData instance
        
        Arguments:
          chipseq_file: (optional) the name of an input file to read
            the ChIP-seq data from.
        """
        self.data = []
        if chipseq_file:
            self.loadDataFromFile(chipseq_file)

    def loadDataFromFile(self,chipseq_file):
        """Read data from a file and populate the object

        Arguments:
          chipseq_file: the name of the input file to read ChIP-seq data from.
        """
        fp = open(chipseq_file,'rU')
        for line in fp:
            # Skip lines that start with a # symbol
            if line.startswith('#'):
                logging.debug("ChIP file: skipped line: %s" % line.strip())
                continue
            # Lines are tab-delimited and have at least 3 columns:
            # chr  start  end
            items = line.strip().split('\t')
            if len(items) < 3:
                logging.warning("ChIP file: skipped line: %s" % line.strip())
                logging.warning("Insufficient number of fields (%d)" % \
                                    len(items))
                continue
            # Check that items in 2nd and 3rd columns are digits
            if not items[1].isdigit() or not items[2].isdigit():
                logging.warning("ChIP file: skipped line: %s" % line.strip())
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
            # Store in a new ChIPSeqDataLine object
            dataline = ChIPSeqDataLine(items[0],
                                       items[1],
                                       items[2])
            self.data.append(dataline)
        fp.close()
        # Return a reference to this object
        return self

    def addDataLine(self,dataline):
        """Append a line of data to the ChIPSeqData object

        Arguments:
          dataline: a ChIPSeqDataLine instance.
        """
        self.data.append(dataline)

    def isSummit(self):
        """Check whether ChIP data consists of summits only

        Checks the difference between start and end positions for
        ChIPSeqDataLines - if all differences are equal to 1 then
        returns True, indicating that the peaks are described as
        summits; otherwise returns False.
        """
        for data in self.data:
            if (data.end - data.start) > 1:
                return False
        return True

    def filterByChr(self,matchChr):
        """Return a subset of data filtered by specified chromosome name

        Returns a new ChIPSeqData object containing only the data from
        the current object which matches the specified criteria.
        """
        # Make a new (empty) ChIPSeqData object
        chip_data_subset = ChIPSeqData()
        # Populate with only the matching data lines
        for dataline in self.data:
            if dataline.chr == matchChr:
                chip_data_subset.addDataLine(dataline)
        return chip_data_subset

    def filterByPosition(self,limit1,limit2,exclude_limits=False):
        """Return a subset of data filtered by peak position

        Returns a new ChIPSeqData object containing only the data from
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
        # Make a new (empty) ChipSeqData object
        chip_data_subset = ChIPSeqData()
        # Populate with only the matching data lines
        for dataline in self.data:
            position = dataline.start
            if exclude_limits:
                if lower < position and position < upper:
                    rna_data_subset.addDataLine(dataline)
            else:
                if lower <= position and position <= upper:
                    chip_data_subset.addDataLine(dataline)
        return chip_data_subset

    def sortByDistanceFrom(self,position):
        """Sort the data into order based on distance from a position
    
        Sorts the ChIPSeqDataLines into order of absolute distance of
        their start to the specified position (closest first).
        
        Note that this operates on the current object.
        """
        self.data = sorted(self.data,
                           key=lambda record: abs(record.start - position))
        return self

    def __getitem__(self,key):
        return self.data[key]

    def __len__(self):
        return len(self.data)

class ChIPSeqDataLine:
    """Class for storing a line of ChIP-seq data (ChIP peak)

    Access the data from the line using the object's properties:

      chr
      start
      end

    NB THIS CLASS MAKES THE ASSUMPTION THAT THE CHIP PEAK DATA IS
    THE SUMMIT..

    There are also convenience methods (e.g. insideRegion).
    """
    def __init__(self,chip_chr,start,end):
        self.chr = chip_chr.strip('"\'')
        self.start = int(start)
        self.end = int(end)

    def __repr__(self):
        return "%s\t%s\t%s" % (self.chr,
                               self.start,
                               self.end)

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
        considered part of the region."""
        if limit1 > limit2:
            upper,lower = limit1,limit2
        else:
            upper,lower = limit2,limit1
        if not exclude_limits:
            return (lower <= self.start and self.start <= upper)
        else:
            return (lower < self.start and self.start < upper)
