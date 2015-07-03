#!/bin/env python
#
#     RNASeq.py: classes for handling RNA-seq data
#     Copyright (C) University of Manchester 2011-15 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
RNASeq.py

Classes for handling RNA-seq transcript/feature data.

"""

import logging
from distances import closestDistanceToRegion

class RNASeqData:
    """Class for storing RNA-seq data

    RNA-seq data consists of genes/transcripts/isomers, which are
    stored individually in RNASeqDataLine objects. This class is a
    container for a collection of RNASeqDataLine objects and provides
    methods to operate on the collection, by creating subsets by
    filtering, and sorting the data based on various criteria.
    """
    def __init__(self,rnaseq_file=None):
        """Create a new RNASeqData instance

        Raises an exception if there are errors in the input file data
        (non-numeric fields for start/end positions, end positions
        occurring before start positions, or illegal strand values).
        
        Arguments:
          rnaseq_file: (optional) the name of an input file to read
            the RNA-seq data from.
        """
        self.data = []
        if rnaseq_file:
            self.loadDataFromFile(rnaseq_file)

    def loadDataFromFile(self,rnaseq_file):
        """Read data from a file and populate the object

        Arguments:
          rnaseq_file: the name of the input file to read RNA-seq data from.
        """
        # Local flags etc
        line_index = 0
        critical_error = False
        # Read in data from file
        fp = open(rnaseq_file,'rU')
        for line in fp:
            # Increment index
            line_index += 1
            # Skip lines starting with #
            if line.startswith('#'):
                logging.debug("RNA file: skipped line: %s" % line.strip())
                continue
            # Lines are tab-delimited and have at least 5 columns:
            # ID  chr  start  end  strand
            items = line.strip().split('\t')
            if len(items) < 5:
                logging.warning("RNA file: skipped line: %s" % line.strip())
                logging.warning("Insufficient number of fields (%d)" % \
                                    len(items))
                continue
            # Check line is valid i.e. start and stop should be
            # numbers, strand should be + or -
            problem_fields = []
            if not items[2].isdigit(): problem_fields.append(2)
            if not items[3].isdigit(): problem_fields.append(3)
            if not (items[4] == '+' or  items[4] == '-'): problem_fields.append(4)
            if problem_fields:
                # If this is the first line then assume it's a header and ignore
                if line_index == 1:
                    logging.warning("%s: first line ignored as header: %s" % 
                                    (rnaseq_file,line.strip()))
                else:
                    # Indicate problem field(s)
                    logging.error("%s: critical error line %d: bad values:" %
                                  (rnaseq_file,line_index))
                    logging.error("%s" % line.strip())
                    logging.error("%s" % make_errline(line.strip(),problem_fields))
                    # This is a critical error: update flag
                    critical_error = True
                # Continue to next line
                continue
            elif int(items[2]) >= int(items[3]):
                # Start position is same or higher than end
                logging.error("%s: critical error line %d: 'end' comes before 'start':" %
                              (rnaseq_file,line_index))
                logging.error("%s" % line.strip())
                logging.error("%s" % make_errline(line.strip(),(2,3)))
                # This is a critical error: update flag but continue reading
                critical_error = True
                continue
            # Store in a new RNASeqDataLine object
            dataline = RNASeqDataLine(items[0],
                                      items[1],
                                      items[2],
                                      items[3],
                                      items[4])
            # Additional flag
            if len(items) >= 6:
                # Is column 6 a flag?
                try:
                    flag_value = int(items[5])
                    if flag_value != 0 and flag_value != 1:
                        flag_value = None
                except ValueError:
                    flag_value = None
                # Store value
                dataline.flag = flag_value

            # Store data
            self.data.append(dataline)
        fp.close()
        # Deal with postponed critical errors
        if critical_error:
            raise Exception, "critical error(s) in '%s'" % rnaseq_file
        # Return a reference to this object
        return self

    def addDataLine(self,dataline):
        """Append a line of data to the RNASeqData object

        Arguments:
          dataline: an RNASeqDataLine instance.
        """
        self.data.append(dataline)

    def filterByChr(self,matchChr):
        """Return a subset of data filtered by specified chromosome name

        Returns a new RNASeqData object containing only the data from
        the current object which matches the specified criteria.
        """
        # Make a new (empty) RNASeqData object
        rna_data_subset = RNASeqData()
        # Populate with only the matching data lines
        for dataline in self.data:
            if dataline.chr == matchChr:
                rna_data_subset.addDataLine(dataline)
        return rna_data_subset

    def filterByStrand(self,matchStrand):
        """Return a subset of data filtered by specified strand

        Returns a new RNASeqData object containing only the data from
        the current object which matches the specified criteria.
        """
        # Make a new (empty) RNASeqData object
        rna_data_subset = RNASeqData()
        # Populate with only the matching data lines
        for dataline in self.data:
            if dataline.strand == matchStrand:
                rna_data_subset.addDataLine(dataline)
        return rna_data_subset

    def filterByFlag(self,matchFlag):
        """Return a subset of data filtered by flag value

        Returns a new RNASeqData object containing only the data from
        the current object which matches the specified criteria.

        Note that if there is no flag (the "isFlagged()" function returns
        False) then an empty set will be returned.
        """
        # Make a new (empty) RNASeqData object
        rna_data_subset = RNASeqData()
        # Populate with only the matching data lines
        for dataline in self.data:
            if dataline.flag == matchFlag:
                rna_data_subset.addDataLine(dataline)
        return rna_data_subset

    def filterByTSS(self,limit1,limit2,exclude_limits=False):
        """Return a subset of data filtered by TSS position

        Returns a new RNASeqData object containing only the data from
        the current object where the TSS positions fall within a
        region defined by upper and lower limits.

        limits can be supplied in either order (i.e. highest/lowest
        or lowest/highest).

        If exclude_limits is False (the default) then TSS positions
        that fall exactly on one of the boundaries are counted as
        being within the region; if it is True then these TSS
        positions will not be considered to lie inside the region.
        """
        # Sort out upper and lower limits
        if limit1 > limit2:
            upper,lower = limit1,limit2
        else:
            upper,lower = limit2,limit1
        # Make a new (empty) RNASeqData object
        rna_data_subset = RNASeqData()
        # Populate with only the matching data lines
        for dataline in self.data:
            TSS = dataline.getTSS()
            if exclude_limits:
                if lower < TSS and TSS < upper:
                    rna_data_subset.addDataLine(dataline)
            else:
                if lower <= TSS and TSS <= upper:
                    rna_data_subset.addDataLine(dataline)
        return rna_data_subset

    def sortByDistanceFrom(self,position):
        """Sort the data into order based on distance from a position
    
        Sorts the RNASeqDataLines into order of absolute distance of
        their TSS to the specified position (closest first).
        
        Note that this operates on the current object.
        """
        self.data = sorted(self.data,
                           key=lambda record: abs(record.getTSS()-position))
        return self

    def sortByClosestEdgeTo(self,position1,position2=None):
        """Sort the data into order based on closest edge (TSS or TES)

        Sorts the RNASeqDataLines into order of smallest absolute distance
        to the specified position (closest first), considering both TSS
        and TES, using the getClosestEdgeDistanceTo method of the
        RNASeqDataLine class.
        
        Note that this operates on the current object."""
        self.data = sorted(self.data,
                           key=lambda record: 
                           record.getClosestEdgeDistanceTo(position1,
                                                           position2))
        return self

    def sortByClosestTSSTo(self,position1,position2=None):
        """Sort the data into order based on closest edge to TSS

        Sorts the RNASeqDataLines into order of smallest absolute distance
        to the specified position (closest first) to the TSS position,
        using the getClosestTSSDistanceTo method of the RNASeqDataLine class.
        
        Note that this operates on the current object."""
        self.data = sorted(self.data,
                           key=lambda record:
                               record.getClosestTSSDistanceTo(position1,
                                                              position2))
        return self

    def isFlagged(self):
        """Check whether RNA-seq data includes flags

        Checks whether all the RNA-seq records also have a valid flag
        associated with them - if yes then returns True (indicating the
        dataset as a whole is flagged), otherwise returns False.
        """
        # Check all data and look for any None flags
        for data in self.data:
            if data.flag is None:
                return False
        # All flags valid
        return True

    def __getitem__(self,key):
        return self.data[key]

    def __len__(self):
        return len(self.data)

class RNASeqDataLine:
    """Class for storing a line of RNA-seq data (gene/transcript/isomer)

    Access the data from the line using the object's properties:

      id
      chr
      start
      end
      strand

    There are also convenience methods (getTSS, getTES, getPromoterRegion)
    and methods for calculating various distances.
    """
    def __init__(self,rna_id,rna_chr,start,end,strand):
        self.id = rna_id
        self.chr = rna_chr
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.flag = None

    def __repr__(self):
        items = [self.id,
                 self.chr,
                 str(self.start),
                 str(self.end),
                 self.strand]
        if self.flag != None:
            items.append(str(self.flag))
        return '\t'.join(items)

    def getTSS(self):
        """Return the TSS coordinate

        TTS (transcription start site) is the start position for a +ve 
        strand, or end for a -ve strand.
        """
        if self.strand == '+':
            return self.start
        elif self.strand == '-':
            return self.end
        # FIXME Should raise an exception
        return None

    def getTES(self):
        """Return the TES coordinate

        TES (transcription end site) is the start position for a +ve
        strand, or end for a -ve strand.
        """
        if self.strand == '+':
            return self.end
        elif self.strand == '-':
            return self.start
        # FIXME Should raise an exception
        return None

    def containsPosition(self,coordinate):
        """Check whether a coordinate is within the gene coordinates

        Returns True if coordinate lies within start and end, False
        otherwise.
        """
        return (self.start <= coordinate and coordinate <= self.end)

    def getClosestTSSDistanceTo(self,position1,position2=None,
                                zero_inside_region=False):
        """Return distance from TSS to a coordinate or region

        For a single specified position, return the absolute distance
        between the position and the TSS.
        
        If a second position is given (specifying a region) then return
        smallest absolute distance of (TSS,position1) and (TSS,position2).

        By default there is no special treatment when the TSS lies inside
        the region specified by two positions; to return zero distance in
        these cases, set the 'zero_inside_region' argument to True.
        """
        return closestDistanceToRegion(self.getTSS(),
                                       position1,position2,
                                       zero_inside_region)

    def getClosestTESDistanceTo(self,position1,position2=None,
                                zero_inside_region=False):
        """Return distance from TES to a coordinate or region

        For a single specified position, return the absolute distance
        between the position and the TES.
        
        If a second position is given (specifying a region) then return
        smallest absolute distance of (TES,position1) and (TES,position2).

        By default there is no special treatment when the TES lies inside
        the region specified by two positions; to return zero distance in
        these cases, set the 'zero_inside_region' argument to True.
        """
        return closestDistanceToRegion(self.getTES(),
                                       position1,position2,
                                       zero_inside_region)

    def getClosestEdgeDistanceTo(self,position1,position2=None,
                                 zero_inside_region=False):
        """Return closest edge distance to a coordinate or region

        For a single specified position, the closest edge is whichever
        of the TSS or TES is nearest (smallest absolute distance) from
        that position i.e. the smallest distance of (TSS,position) and
        (TES,position).
        
        If a second position is given (specifying a region) then
        the closest edge is whichever of the TSS/TES is closest to
        either position1 or position2 i.e. the smallest distance of
        (TSS,position1), (TES,position1), (TSS,position2) and
        (TES,position2).

        By default there is no special treatment when either the TSS
        or TES lie inside the region specified by two positions; to
        set this to zero, set the 'zero_inside_region' argument to
        True.
        """
        return min(self.getClosestTSSDistanceTo(position1,
                                                position2,
                                                zero_inside_region),
                   self.getClosestTESDistanceTo(position1,
                                                position2,
                                                zero_inside_region))

    def getPromoterRegion(self,to_TSS,from_TSS):
        """Return the coordinates of the promoter region

        The promoter region is a region of coordinates around the
        TSS of a gene, defined by the supplied distances 'to_TSS'
        (the distance downstream from the TSS) and 'from_TSS' (the
        distance upstream from the TSS).

        Returns a tuple containing the start and end coordinates
        defining the promoter region.
        """
        if self.strand == '+':
            return (self.getTSS() - to_TSS,
                    self.getTSS() + from_TSS)
        else:
            return (self.getTSS() + to_TSS,
                    self.getTSS() - from_TSS)