#!/bin/env python
#
#     RnaChipIntegrator.py: analyse RNA-seq and ChIP-seq data
#     Copyright (C) University of Manchester 2011 Peter Briggs
#
########################################################################
#
# RnaChipIntegrator.py
#
#########################################################################

"""RnaChipIntegrator.py

Provides classes and functions for integrating/analysing RNA-seq and
ChIP-seq data. __main__ provides a program for running these analyses
from the command line.

Usage: RnaChipIntegrator.py [OPTIONS] <rna-data> <chip-data>

RNA-seq data has at least 5 columns of data:

# ID  chr  start  end  strand

Optionally there might be a 6th column (flag)

ChIP-seq data has 3 columns of data:
# chr  start  end"""

#######################################################################
# Module metadata
#######################################################################

__version__ = "0.2.0"

#######################################################################
# Import modules that this module depends on
#######################################################################

import sys
import os
import optparse

# Set default logging level and output
import logging
logging.getLogger().setLevel(logging.ERROR)
logging.basicConfig(format='%(levelname)s: %(message)s')

# Get the Spreadsheet module
try:
    import Spreadsheet
except ImportError:
    logging.error("Failed to import the Spreadsheet module")
    logging.error("Set your PYTHONPATH to include the directory with this module, or get the")
    logging.error("latest version from github via:")
    logging.error("https://github.com/fls-bioinformatics-core/genomics/blob/master/share/Spreadsheet.py")
    logging.error("and ensure that the underlying xlwt, xlrd and xlutils libraries are installed")
    sys.exit(1)

#######################################################################
# Class definitions
#######################################################################

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

class AnalysisResult:
    """Class to hold results from RNA-seq/ChIP-seq analysis

    This class is a general container for results from an analysis.

    Each entry is a Python dictionary which associates a ChIPSeqDataLine
    and an RNASeqDataLine pair along with any additional arbitrary data,
    via the addResult method, e.g.

    >>> result = AnalysisResult()
    >>> result.addResult(chip,rna,distance_to_TSS=12345,in_the_gene='YES')

    'distance_to_TSS' and 'in_the_gene' are arbitrary stored data.

    The results can be written to a file using the 'output' method,
    specifying a file name and a header line which names the arbitrary
    stored data items to write in each column, e.g.

    >>> result.output('analysis.out',('distance_to_TSS','in_the_gene'))
    """
    def __init__(self):
        """Create a new AnalysisResult instance.

        Use the 'addResult' method to populate."""
        # List of results
        self.results = []
        # List of fields
        self.fields = []

    def addResult(self,chip_data,rna_data,**args):
        """Add a result to the AnalysisResult instance.

        A result associates a ChIP peak with a transcript,

        Arguments:
          chip_data: a ChIPSeqDataLine
          rna_data: an RNASeqDataLine
          
        In addition this method will accept arbitrary named arguments
        which can be supplementary data associated with the result.
        """
        # Create a new dictionary
        result = {'chip_seq': chip_data,
                  'rna_seq': rna_data }
        # Add the arbitrary items
        for arg in args:
            if arg not in self.fields:
                self.fields.append(arg)
            result[arg] = args[arg]
        self.results.append(result)

    def output(self,filename,header,pad=''):
        """Write results out to a file

        The results are written to a tab-delimited file named 'filename',
        with the items in 'header' defining the content of each line.
        These header items should be a subset of the custom items
        stored by calls to addResult.

        Arguments:
          filename: name of the file to write to
          header: list of items from each result line to be written
          pad: (optional) string to use for padding empty fields
        """
        fp = open(filename,'w')
        # Write header line
        fp.write('#'+'\t'.join(header)+'\n')
        # Write results
        for result in self.results:
            items = []
            for field in header:
                try:
                    items.append(str(result[field]))
                except KeyError:
                    items.append(pad)
            fp.write('\t'.join(items)+'\n')
        fp.close()

    def output_xls(self,xls,title,header,pad=''):
        """Write results to an XLS spreadsheet
        """
        ws = xls.addSheet(title)
        ws.addText('#'+'\t'.join(header))
        # Keep track of line number and sheet number
        row_index = 1
        nsheet = 1
        # Loop over all results
        for result in self.results:
            items = []
            for field in header:
                try:
                    item = str(result[field])
                except KeyError:
                    item = pad
                # Split too-long items into multiple items/columns
                # This is because the Spreadsheet writer truncates cells
                # that exceed the spreadsheet cell character limit
                char_limit = Spreadsheet.MAX_LEN_WORKSHEET_CELL_VALUE
                while len(item) > char_limit:
                    logging.warning("Split value across multiple cells in sheet '%s' row %d"
                                    % (ws.title,row_index))
                    # Split on ';'
                    try:
                        # Locate nearest semicolon to the character limit
                        i = item[:char_limit].rindex(';')
                        items.append(item[:i])
                        item = item[i:].strip(';')
                    except ValueError:
                        # Unable to locate semicolon so split on the
                        # character limit
                        items.append(item[:char_limit])
                        item = item[char_limit:]
                items.append(item)
            ws.addText('\t'.join(items))
            # Update line number and check if the maximum number of lines has been exceeded
            row_index += 1
            if row_index == Spreadsheet.MAX_NUMBER_ROWS_PER_WORKSHEET:
                logging.warning("Maximum number of rows in XLS sheet '%s' exceeded (%d)" %
                                (ws.title,Spreadsheet.MAX_NUMBER_ROWS_PER_WORKSHEET))
                # Make a new spreadsheet for the excess rows
                nsheet += 1
                ws = xls.addSheet("%s(%d)" % (title,nsheet))
                logging.warning("Created new sheet '%s' to store additional results" %
                                ws.title)
                # Add the header
                ws.addText('#'+'\t'.join(header))
                # Reset the line index counter
                row_index = 1

    def __getitem__(self,key):
        return self.results[key]

    def __len__(self):
        return len(self.results)

#######################################################################
# Utility Functions
#######################################################################

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
        lower, upper = wide[0],wide[1]
    else:
        upper,lower = wide[1],wide[0]
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

def make_errline(line,bad_fields=[]):
    """Return an 'error line' indicating problem fields in a string

    Given a tab-delimited line and a list of integer indices
    indicating which fields in the line have problems, this function
    returns a tab-delimited string where the original fields are
    replaced by either spaces or '^' characters.

    When printed beneath the original line, the '^'s indicate which
    fields are 'bad' according to the supplied indices, e.g.

    Input line: 'good    good    bad    bad    good'
    Error line: '                ^^^    ^^^        '

    Arguments:
      line: string where tabs delimit fields
      bad_fields: list of integer indices corresponding to 'bad'
        values in 'line'

    Returns:
      Tab-delimited 'error line' to be printed beneath the original
      line, to indicate which fields are 'bad'.
    """
    # Indicate problem field(s)
    errline = []
    items = line.split('\t')
    for i in range(len(items)):
        if i in bad_fields:
            errline.append("^"*len(items[i]))
        else:
            errline.append(" "*len(items[i]))
    return '\t'.join(errline)

#######################################################################
# Analysis Functions
#######################################################################

def AnalyseClosestTranscriptsToPeaksInEachDirection(chip_seq,rna_seq):
    """Get closest up- and down-stream transcripts for + and - strands to ChIP peak"
    """
    print "Analysis #1a: ChIP-seq perspective"
    print "Get closest upstream and downstream transcripts to each ChIP peak"
    print "on both + and - strands"
    print ""
    results = []
    for chip_data in chip_seq:
        print "\tMatching to ChIP peak: %s" % chip_data
        # Initialise the hits
        print "\t%d total matches in RNA-seq data" % \
            len(rna_seq.filterByChr(chip_data.chr))
        # Process the subset of RNA-seq data with matching chromosome
        # on the + strand
        closest_upstream_on_plus_strand = None
        closest_downstream_on_plus_strand = None
        print "\t%d matches on + strand" % \
            len(rna_seq.filterByChr(chip_data.chr).filterByStrand('+'))
        for rna_data in rna_seq.filterByChr(chip_data.chr).filterByStrand('+'):
            # Get the distance from the ChIP peak
            distance = rna_data.getTSS() - chip_data.start
            # Check if this is closer than anything seen previously
            # for each case
            if distance >= 0:
                # Upstream (i.e. +ve distance)
                closest_upstream_on_plus_strand = \
                    GetNearestTranscriptToPeak(closest_upstream_on_plus_strand,\
                                                   rna_data,chip_data)
            else:
                # Downstream (i.e. -ve distance)
                closest_downstream_on_plus_strand = \
                    GetNearestTranscriptToPeak(closest_downstream_on_plus_strand,\
                                                   rna_data,chip_data)
        # Process the subset of RNA-seq data with matching chromosome
        # on the - strand
        closest_upstream_on_minus_strand = None
        closest_downstream_on_minus_strand = None
        print "\t%d matches on - strand" % \
            len(rna_seq.filterByChr(chip_data.chr).filterByStrand('-'))
        for rna_data in rna_seq.filterByChr(chip_data.chr).filterByStrand('-'):
            # Get the distance from the ChIP peak
            distance = rna_data.getTSS() - chip_data.start
            # Check if this is closer than anything seen previously
            # for each case
            if distance >= 0:
                # Upstream (i.e. +ve distance)
                closest_upstream_on_minus_strand =  \
                    GetNearestTranscriptToPeak(closest_upstream_on_minus_strand,
                                               rna_data,chip_data)
            else:
                # Downstream (i.e. -ve distance)
                closest_downstream_on_minus_strand = \
                    GetNearestTranscriptToPeak(closest_downstream_on_minus_strand,
                                              rna_data,chip_data)
        # Store the data about the hits to output later
        nearest_positions = {
            "upstream +": closest_upstream_on_plus_strand,
            "downstream +": closest_downstream_on_plus_strand,
            "upstream -": closest_upstream_on_minus_strand,
            "downstream +": closest_downstream_on_minus_strand,
            }
        for hit in nearest_positions.keys():
            rna_data = nearest_positions[hit]
            if not rna_data:
                # No data stored for this position
                print "*** No hit for %s ***" % hit
                continue
            # Derive the distances
            distance_to_TSS = rna_data.getTSS() - chip_data.start
            distance_to_TES = rna_data.getTES() - chip_data.start
            # ChIP peak is in the gene?
            if rna_data.containsPosition(chip_data.start):
                in_the_gene = "YES"
            else:
                in_the_gene = "NO"
            # Add to results list
            results.append([chip_data.chr,
                            chip_data.start,
                            rna_data.id,
                            distance_to_TSS,
                            distance_to_TES,
                            rna_data.strand,
                            in_the_gene])

    # Write the results to file
    results.insert(0,["#chr","start","ID","distance_to_TSS","distance_to_TES",
                      "strand","in_the_gene"])
    output_results("RNA-seq_to_ChIP-seq.txt",results)

def AnalyseNearestTSSToSummits(chip_seq,rna_seq,max_distance,
                               max_closest=4,
                               filename=None,
                               xls=None):
    """Find nearest RNA transcripts to a set of ChIP peaks

    Given a set of ChIP peaks in a ChIPSeqData object, and a set of
    gene transcripts in an RNASeqData object, this analysis aims to
    match genes to peaks.

    The procedure is:

    - For each ChIP peak:
    -- Consider only genes that are flagged as "significant" and
       with the same chromosome
    -- Report the nearest genes (by TSS) within a cutoff distance
       of the peak position (up to 'max_closest' hits)
    -- For each hit:
    --- Report all the genes (by TSS) from the full list of genes
        which lie between the ChIP peak and the TSS of the hit

    Arguments:
      chip_seq: ChIP peaks in ChIPSeqData object
      rna_seq: RNA-seq gene transcripts in RNASeqData object
      max_distance: cutoff distance for reporting nearest genes
        to ChIP peaks (in units of bases)
      max_closest: (optional) maximum number of peaks to report
        (defaults to 4)
      filename: (optional) if not None then specifies the file
        name to write the results to as tab-delimited data
      xls: (optional) if not None then specifies the XLS file name
        to add the results to as a new sheet
    """
    # Loop over ChIP peaks and sort RNA transcripts into order
    # for each based on the absolute distance of their TSS from the peak
    # After sorting the top results will be the closest
    logging.debug("Starting AnalyseNearestTranscriptsToPeaks:")

    # Check that we have summit data for ChIP peaks
    if not chip_seq.isSummit():
        logging.warning("The supplied ChIP data only defines regions")
        logging.warning("This analysis is intended to work with summit data")

    # Create subset of "significant" i.e. flagged RNA-seq data
    if rna_seq.isFlagged():
        significant_rna_seq = rna_seq.filterByFlag(matchFlag=1)
        logging.debug("%d RNA-seq records after filtering on flag" % \
                          len(significant_rna_seq))
    else:
        significant_rna_seq = rna_seq

    # Run the analysis on the significant RNA-seq transcripts
    results = AnalysisResult()
    for chip_data in chip_seq:
        logging.debug("\t%s" % chip_data)
        rna_chr = significant_rna_seq.\
            filterByChr(chip_data.chr).\
            filterByTSS(chip_data.start-max_distance,\
                            chip_data.start+max_distance).\
                            sortByDistanceFrom(chip_data.start)
        rna_chr_full = rna_seq.filterByChr(chip_data.chr)
        closest = rna_chr[:max_closest]
        i = 0
        for rna_data in closest:
            logging.debug("\t\t%s\t%d" % (rna_data,
                                          (rna_data.getTSS()-chip_data.start)))
            # Get the transcripts from the full list which
            # lie between the transcript TSS and the ChIP peak
            # These can be on either strand
            transcripts = rna_chr_full.filterByTSS(chip_data.start,
                                                   rna_data.getTSS(),
                                                   exclude_limits=True)
            # Put them into distance order from the peak
            transcripts.sortByDistanceFrom(chip_data.start)
            # Make a list of just the ids
            transcript_ids = []
            for transcript in transcripts:
                transcript_ids.append(transcript.id)
                logging.debug("\t\t\t%s" % transcript.id)
            # Derive the distances
            distance_to_TSS = rna_data.getTSS() - chip_data.start
            distance_to_TES = rna_data.getTES() - chip_data.start
            # ChIP peak is in the gene?
            if rna_data.containsPosition(chip_data.start):
                in_the_gene = "YES"
            else:
                in_the_gene = "NO"
            # "Nearest" column reports "1 of 4", "2 of 4" etc
            i += 1
            nearest = "%d of %d" % (i,len(closest))
            # Add to the analysis result
            results.addResult(chip_data,rna_data,
                              chr=chip_data.chr,
                              start=chip_data.start,
                              ID=rna_data.id,
                              nearest=nearest,
                              TSS=rna_data.getTSS(),
                              distance_to_TSS=distance_to_TSS,
                              distance_to_TES=distance_to_TES,
                              strand=rna_data.strand,
                              in_the_gene=in_the_gene,
                              transcripts_inbetween=len(transcripts),
                              transcript_ids_inbetween=\
                                  ';'.join(transcript_ids))
        # Report peaks with no significant genes in the cut-off region
        if len(closest) == 0:
            results.addResult(chip_data,None,
                              chr=chip_data.chr,
                              start=chip_data.start)
            logging.debug("\t\tNo transcripts found")
        logging.debug("")
    # Write the results to file
    if filename:
        results.output(filename,('chr','start','ID','nearest','TSS',
                                 'distance_to_TSS',
                                 'distance_to_TES',
                                 'strand','in_the_gene',
                                 'transcripts_inbetween',
                                 'transcript_ids_inbetween'))
    # Write the results to a spreadsheet
    if xls: results.output_xls(xls,'TSSToSummits',
                               ('chr','start','ID','nearest','TSS',
                                'distance_to_TSS',
                                'distance_to_TES',
                                'strand','in_the_gene',
                                'transcripts_inbetween',
                                'transcript_ids_inbetween'))
    # Return results object
    return results

def AnalyseNearestPeaksToTranscripts(rna_seq,chip_seq,window_width,
                                     filename=None,xls=None):
    """Find nearest ChIP peaks to RNA-seq gene transcripts

    Given a set of ChIP peaks in a ChIPSeqData object, and a set of
    gene transcripts in an RNASeqData object, this analysis aims to
    match peaks to genes.

    The procedure is:

    - For each RNA-seq gene transcript:
    -- If the gene is flagged as "significant" then:
    --- Report the ChIP peaks on the same chromosome which lie
        within +/- window_width of the gene TSS

    Note that all genes are reported in the output file.

    Arguments:
      rna_seq: RNA-seq gene transcripts in RNASeqData object
      chip_seq: ChIP peaks in ChIPSeqData object
      window_width: cutoff distance for reporting nearest ChIP peaks
        to RNA-seq gene TSS positions (in units of bases)
      filename: (optional) if not None then specifies the file
        name to write the results to as tab-delimited data
      xls: (optional) if not None then specifies the XLS file name
        to add the results to as a new sheet
    """
    logging.debug("Starting AnalyseNearestPeaksToTranscripts:")
    # Check that we have summit data for ChIP peaks
    if not chip_seq.isSummit():
        logging.warning("The supplied ChIP data only defines regions")
        logging.warning("This analysis is intended to work with summit data")
    # Get flag status of RNA-seq data
    rna_seq_is_flagged = rna_seq.isFlagged()
    # Do the analysis
    results = AnalysisResult()
    max_peaks = 0
    for rna_data in rna_seq:
        logging.debug("\t%s" % rna_data)
        # Define a "window" either side of the TSS
        window = (rna_data.getTSS()-window_width,
                  rna_data.getTSS()+window_width)
        logging.debug("\t\tWindow: %d - %d" % (window[0],window[1]))
        # Get ChIP peaks on the matching chromosome that lie in this window
        # Only do this for the significant (i.e. flagged) genes, unless
        # the data is not flagged
        if not rna_seq_is_flagged or rna_data.flag == 1:
            chip_peaks = chip_seq.\
                filterByChr(rna_data.chr).filterByPosition(window[0],window[1])
            chip_peaks.sortByDistanceFrom(rna_data.getTSS())
        else:
            # No hits - empty result object
            chip_peaks = ChIPSeqData()
            logging.debug("\t\tNo peaks found")
        # Store results
        result = { 'ID': rna_data.id,
                   'chr_RNA': rna_data.chr,
                   'start': rna_data.start,
                   'end': rna_data.end,
                   'strand': rna_data.strand,
                   'differentially_expressed': rna_data.flag,
                   'number_of_peaks': len(chip_peaks) }
        for i in range(len(chip_peaks)):
            chip_data = chip_peaks[i]
            # Derive the distance
            distance = rna_data.getTSS() - chip_data.start
            # Add to results for output
            result['chr_ChIP_'+str(i+1)] = chip_data.chr
            result['summit_'+str(i+1)] = chip_data.start
            result['distance_'+str(i+1)] = distance
            logging.debug("\t\t%s (%d)" % (chip_data,distance))
        logging.debug("")
        # Add line to results to be written to the output file
        results.addResult(chip_peaks,rna_data,**result)
        # Keep track of the maximum number of peaks for output
        max_peaks = max(len(chip_peaks),max_peaks)
    # Construct header line for output
    if filename or xls:
        header = ["ID","chr_RNA","start","end","strand",
                  "differentially_expressed","number_of_peaks"]
        for i in range(max_peaks):
            header.extend(["chr_ChIP_"+str(i+1),
                           "summit_"+str(i+1),
                           "distance_"+str(i+1)])
    # Write to output file
    if filename: results.output(filename,header,pad='---')
    # Write the results to a spreadsheet
    if xls: results.output_xls(xls,'PeaksToTranscripts',
                               header,pad='---')
    # Finished
    return results

def AnalyseNearestTranscriptsToPeakEdges(chip_seq,rna_seq,
                                         promoter_region=(10000,2500),
                                         max_closest=4,
                                         max_distance=0,
                                         TSS_only=False,
                                         filename=None,
                                         xls=None):
    """Find nearest RNA transcripts to a set of ChIP peaks based on edges

    Given a set of ChIP peaks in a ChIPSeqData object, and a set of
    gene transcripts in an RNASeqData object, this analysis aims to
    match genes to peaks.

    The procedure is:

    - For each ChIP peak:
    -- Consider only genes on the same chromosome
    -- Report the nearest genes by finding either the nearest TSS or TES, or
       just the nearest TSS, to the peak edges (start/stop) (up to four hits)
       (Both TSS and TES are used if the 'TSS_only' flag is set to False;
       otherwise only TSS is used).
    -- For each hit:
    --- Report whether the peak overlaps the transcript region
    --- Report whether the peak overlaps the promoter region

    Arguments:
      chip_seq: ChIP peaks in ChIPSeqData object
      rna_seq: RNA-seq gene transcripts in RNASeqData object
      promoter_region: (optional) a tuple (leading,trailing) which
        specifies the promoter region by the "leading" (upstream) and
        "trailing" (downstream) edges of the region with respect to the
        gene TSS. Defaults to (10000,2500).
      max_closest: (optional) maximum number of peaks to report
        (defaults to 4).
      max_distance: (optional) set cutoff gene-to-edge distance, beyond
        which genes are not reported. Defaults to zero (don't apply a
        cutoff and report all genes regardless of distance).
      TSS_only: if set to True then determine closest distance from TSS
        only, otherwise use both TSS and TES (the default).
      filename: (optional) if specified then write output to file with
        with this name; also write a 'summary' file with just the top
        hit for each ChIP peak.
      xls: (optional) if not None then specifies the XLS file name
        to add the results to as a new sheet
    """
    logging.debug("Starting AnalyseNearestTranscriptsToPeakEdges:")
    # Check that we have region data for ChIP peaks
    if chip_seq.isSummit():
        logging.warning("The supplied ChIP data only defines summits")
        logging.warning("This analysis is intended to work with regions")
    # Check mode of operation
    if not TSS_only:
        logging.debug("Use both TSS and TES in analysis (TSS_only = %s)" %
                      TSS_only)
    else:
        logging.debug("Use only TSS in analysis (TSS_only = %s)" % TSS_only)
    results = AnalysisResult()
    # Loop over all peaks
    for chip_peak in chip_seq:
        logging.debug("\t%s" % chip_peak)
        # Sort transcripts into order of the nearest TSS or TES to
        # the edge of the peak/binding region
        # Note that this method does not discriminate between situations
        # where the transcript lies partially or wholly within the
        # binding region
        rna_chr = rna_seq.filterByChr(chip_peak.chr)
        if not TSS_only:
            rna_chr = rna_chr.sortByClosestEdgeTo(chip_peak.start,
                                                  chip_peak.end)
        else:
            rna_chr = rna_chr.sortByClosestTSSTo(chip_peak.start,
                                                 chip_peak.end)
        # Get the closest peaks (filter by distance and max number)
        closest = []
        for rna_data in rna_chr:
            # Get closest edge distances
            distance = rna_data.getClosestEdgeDistanceTo(chip_peak.start,
                                                         chip_peak.end)
            # Apply distance cutoff
            if max_distance > 0 and distance > max_distance:
                # Break out of the loop
                logging.debug("Exceeded maximum distance, stopping")
                break
            distance_to_TSS = rna_data.getClosestTSSDistanceTo(\
                chip_peak.start,chip_peak.end)
            distance_to_TES = rna_data.getClosestTESDistanceTo(\
                chip_peak.start,chip_peak.end)
            # Determine if transcript and peak overlap
            binding_region = (chip_peak.start,chip_peak.end)
            if regions_overlap(binding_region,(rna_data.getTSS(),
                                               rna_data.getTES())):
                overlap_transcript = 1
            else:
                overlap_transcript = 0
            # Determine if promoter region and peak overlap?
            promoter = rna_data.getPromoterRegion(*promoter_region)
            if regions_overlap(binding_region,promoter):
                overlap_promoter = 1
            else:
                overlap_promoter = 0
            # Add to results
            results.addResult(chip_peak,rna_data,
                              chr=chip_peak.chr,
                              start=chip_peak.start,
                              end=chip_peak.end,
                              ID=rna_data.id,
                              TSS=rna_data.getTSS(),
                              TES=rna_data.getTES(),
                              strand=rna_data.strand,
                              dist_closest_edge=distance,
                              dist_TSS=distance_to_TSS,
                              dist_TES=distance_to_TES,
                              overlap_transcript=overlap_transcript,
                              overlap_promoter=overlap_promoter)
            # Report
            logging.debug("\t\t%s (%d,%s,%s)" % \
                             (rna_data,
                              distance,
                              overlap_transcript,
                              overlap_promoter))
            # Check whether we've reached the limit to report
            closest.append(rna_data)
            if len(closest) == max_closest:
                logging.debug("Found %d closest" % max_closest)
                break
        # Finished loop
        if len(closest) == 0:
            logging.debug("\t\tNo transcripts found")
        logging.debug("")
    # Construct header line and summary results for output
    if filename or xls:
        # Header line
        header = ["chr","start","end","ID","strand","TSS","TES",
                  "dist_closest_edge",
                  "dist_TSS","dist_TES",
                  "overlap_transcript",
                  "overlap_promoter"]
        # Summary results: one gene per ChIP peak (i.e. top hit for
        # each peak)
        summary_results = AnalysisResult()
        last_chip_peak = None
        for result in results:
            if result['chip_seq'] != last_chip_peak:
                summary_results.addResult(result['chip_seq'],
                                          result['rna_seq'],
                                          **result)
                last_chip_peak = result['chip_seq']
            else:
                pass
    # Write to output file
    if filename:
        # Full results
        results.output(filename,header)
        # Summary
        summary_filename = os.path.splitext(filename)[0]+"_summary"+\
            os.path.splitext(filename)[1]
        summary_results.output(summary_filename,header)
    # Write the results to a spreadsheet
    if xls:
        # Full results
        if not TSS_only:
            xls_title = 'TranscriptsToPeakEdges'
        else:
            xls_title = 'TSSToPeakEdges'
        results.output_xls(xls,xls_title,header)
        # Summary results
        xls_title += '(summary)'
        summary_results.output_xls(xls,xls_title,header)
            
    # Return (full) analysis results
    return results

#######################################################################
# Descriptions of each method for XLS notes
#######################################################################

xls_notes_for_nearest_TSS_to_summits = \
"""<style font=bold bgcolor=gray25>Nearest TSS to Peak Summits</style>

<style font=bold>Input parameters:</style>
Cutoff distance from peaks\t%d bp

<style font=bold>Description of output fields:</style>
chr\tchromosome
start\tstart position for the peak (assumed to be the summit)
ID\tID for a closest gene/transcript
nearest\ta string of the form "1 of 4", "2 of 3" etc, indicating how many transcripts are listed for the peak, and which one of these the current transcript is.
TSS\tthe TSS position for the gene/transcript
distance_to_TSS\tdistance from the peak start (= summit) to the gene TSS
distance_to_TES\tdistance from the peak start (= summit) to the gene TES
strand\tthe strand direction
in_the_gene\tindicates whether the peak start position (= summit) lies within the gene coordinates (either `YES` or `NO`)
transcripts_inbetween\tnumber of genes lying between the peak and the current gene
transcript_ids_inbetween\tlist of gene names lying between the peak and the current gene
"""

xls_notes_for_nearest_transcripts_to_peak_edges = \
"""<style font=bold bgcolor=gray25>Nearest Transcripts to Peak Edges/Nearest TSS to Peak Edges</style>

<style font=bold>Input parameters:</style>
Promoter region:\t%s bp
Maximum number of transcripts to report\t%d
Cutoff distance from peaks\t%d bp

<style font=bold>Description of output fields:</style>
chr\tchromosome
start\tpeak start position
end\tpeak end position
ID\tID for a closest gene/transcript
strand\tthe strand direction
TSS\tgene TSS position
TES\tgene TES position
dist_closest_edge\tclosest distance between the edges of the peak and gene regions.
dist_TSS\tdistance from the closest edge to the gene TSS.
dist_TES\tdistance from the closest edge to the gene TES.
overlap_transcript\tindicates whether the gene region overlaps the the peak region at any point (1 indicates an overlap, 0 no overlap).
overlap_promoter\tindicates whether the gene promoter region overlaps the peak region at any point (1 indicates an overlap, 0 no overlap).

<style font=bold>NB "summary" pages list only the top hits for each peak or transcript.</style>
"""

xls_notes_for_nearest_peaks_to_transcripts = \
"""<style font=bold bgcolor=gray25>Nearest Peaks to Transcripts</style>

<style font=bold>Input parameters:</style>
Window width:\t%d bp

<style font=bold>Description of output fields:</style>
ID\tgene/transcript ID
chr_RNA\tchromosome
start\tgene start position
end\tgene end position
strand\tthe strand direction
differentially_expressed\tthe value of the "significance flag" supplied on input
number_of_peaks\tthe number of ChIP peaks found within the "window" distance of the gene TSS

Then for each peak (closest first) there are three columns:
chr_ChIP_#\tchromosome (same as `chr_RNA` above)
summit_#\tpeak start (=summit)
distance_#\tdistance from the peak summit to the gene TSS
"""

#######################################################################
# Non-core Functions
#######################################################################

def count_unique_TSS(rna_seq):
    """Count the number of unique TSS positions

    Given an RNASeqData object with a list of transcripts, returns
    the number of unique TSS positions found in the list."""
    unique = []
    for rna_data in rna_seq:
        if not rna_data.getTSS() in unique:
            unique.append(rna_data.getTSS())
    return len(unique)

def output_results(filen,results):
    """Write list of results to tab-delimited file

    Given a (Python) list of results, each item is written to the
    specified file as a set of tab-delimited fields.

    If the file already exists then it will be overwritten, and
    there are no checks on the integrity of the data.
    """
    fp = open(filen,'w')
    for result in results:
        items = [str(x) for x in result]
        fp.write('\t'.join(items)+'\n')
    fp.close()

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":

    # Initialisations
    do_chip_analyses = False
    do_rna_analyses = False
    xls_out = None
    max_distance = 130000
    max_edge_distance = 0
    window_width = 20000
    promoter_region = (10000,2500)
    max_closest = 4

    # Set default logging level
    logging.getLogger().setLevel(logging.INFO)

    p = optparse.OptionParser(usage="%prog [options] RNA-seq_data ChIP-seq_data",
                              version="%prog "+__version__,
                              description=
                              "Implements various methods for reporting the nearest ChIP peaks "+
                              "to genes/gene transcripts from RNA-seq data, and vice versa. "
                              "ChIP-centric analyses report the nearest transcripts to each "+
                              "ChIP peak; RNA-seq-centric analyses report the nearest peaks to "+
                              "each transcript. 'RNA-seq_data' file must contain tab-delimited "+
                              "columns 'ID,chr,start,end,strand[,flag,...]'. 'ChIP-seq_data' "+
                              "file must contain tab-delimited columns 'chr,start,stop' and "+
                              "defines either summits (start/stop differ by 1 bp) or regions "+
                              "(start/stop extend over several bps)."+
                              "The outputs are: one tab-delimited file for each analysis "+
                              "performed (named after the appropriate input file unless "+
                              "overriden by the --project option), and an XLS spreadsheet with "+
                              "one worksheet per analysis.")

    # General options
    p.add_option('--chip',action="store_true",dest="do_chip_analyses",
                     help="Do ChIP-seq-centric analyses")
    p.add_option('--rna',action="store_true",dest="do_rna_analyses",
                     help="Do RNA-seq-centric analyses")
    p.add_option('--project',action="store",dest="basename",
                     help="Set basename for output files; output from each "+
                     "analysis method will use this, with the method name appended"+
                     " (defaults to the input file names)")
    p.add_option('--no-xls',action="store_false",dest="write_xls_out",default=True,
                 help="Don't write an XLS file")
    p.add_option('--debug',action="store_true",dest="debug",
                     help="Verbose output for debugging")

    # Options for NearestTSSToSummits
    group = optparse.OptionGroup(p,"NearestTSSToSummits (ChIP-seq)",
                                 description="For each ChIP peak summit, reports the "+
                                 "transcripts with TSS positions that lie within the specified "+
                                 "cut-off distance of the peak summit.")
    group.add_option('--cutoff',action="store",dest="max_distance",
                     default=max_distance,type='int',
                     help="Maximum distance a transcript TSS can be from each peak "+
                     "before being omitted from the analysis "+
                     "(default %d bp)" % max_distance)
    p.add_option_group(group)

    # Options for NearestPeaksToTranscripts
    group = optparse.OptionGroup(p,"NearestPeaksToTranscripts (RNA-seq)",
                                 description="For each transcript, reports the peaks with summit "+
                                 "positions that lie within the specified 'window' distance of "+
                                 "the transcript TSS.")
    group.add_option('--window',action="store",dest="window_width",
                     default=window_width,type='int',
                     help="Maximum distance a peak can be from each transcript TSS "+
                     "before being omitted from analysis "+
                     "(default %d bp)" % window_width)
    p.add_option_group(group)

    # Options for NearestTranscriptsToPeakEdge/NearestTSSToPeakEdge
    group = optparse.OptionGroup(p,"NearestTranscriptsToPeakEdge/NearestTSSToPeakEdge (ChIP-seq)",
                                 description="For each ChIP peak, reports the transcripts that "+
                                 "lie closest to either 'edge' of the peak region, by "+
                                 "considering the TSS alone (NearestTSSToPeakEdge) or by "+
                                 "considering both the TSS and TES positions "+
                                 "(NearestTranscriptsToPeakEdge).")
    group.add_option('--edge-cutoff',action="store",dest="max_edge_distance",
                     default=max_edge_distance,type='int',
                     help="Maximum distance a transcript TSS can be from the peak "+
                     "edge before being omitted from analysis. Set to "+
                     "zero to indicate no cut-off (default %d bp)" % max_edge_distance)
    group.add_option('--number',action="store",dest="max_closest",
                     default=max_closest,type='int',
                     help="Maximum number of transcripts to report from "+
                     "from the analysis (default %d)" % max_closest)
    group.add_option('--promoter_region',action="store",dest="promoter_region",
                     default="%d,%d" % promoter_region,
                     help="Define promoter region with respect to gene TSS "+
                     "(default -%d to %d bp of TSS)" %  promoter_region)
    p.add_option_group(group)

    # Process the command line
    options,arguments = p.parse_args()

    # Input files
    if len(arguments) != 2:
        p.error("incorrect number of arguments")
    else:
        rnaseq_file = arguments[0]
        chipseq_file = arguments[1]

    # Report version
    p.print_version()

    # Sort out analysis settings
    do_chip_analyses = options.do_chip_analyses
    do_rna_analyses = options.do_rna_analyses
    if not (do_chip_analyses or do_rna_analyses):
        # Neither explicitly requested - do both
        do_chip_analyses = True
        do_rna_analyses = True

    # Handle options
    max_distance = options.max_distance
    window_width = options.window_width
    max_edge_distance = options.max_edge_distance
    max_closest = options.max_closest
    write_xls_out = options.write_xls_out

    # Promoter region
    promoter_region = (abs(int(options.promoter_region.split(',')[0])),
                       abs(int(int(options.promoter_region.split(',')[1]))))

    # Output basename
    if options.basename:
        chip_basename = options.basename + "_peaks"
        rna_basename = options.basename + "_transcripts"
        if write_xls_out:
            xls_out = options.basename + ".xls"
    else:
        rna_basename = os.path.basename(os.path.splitext(rnaseq_file)[0])
        chip_basename = os.path.basename(os.path.splitext(chipseq_file)[0])
        if write_xls_out:
            xls_out = chip_basename + "_summary.xls"

    # Debugging output
    if options.debug: logging.getLogger().setLevel(logging.DEBUG)

    # Report settings
    print "Input transcripts file (RNA-seq) : %s" % rnaseq_file
    print "Input peaks file (ChIP-seq)      : %s" % chipseq_file
    if do_chip_analyses:
        print ""
        print "ChIP analyses:"
        print "\tMaximum cutoff distance  : %d (bp)" % max_distance
        print "\tMaximum edge distance    : %d (bp)" % max_edge_distance
        print "\tMax no. of closest genes : %d" % max_closest
        print "\tPromoter region          : -%d to %d (bp from TSS)" % \
            promoter_region
        print "\tBasename for output files: %s" % chip_basename
    if do_rna_analyses:
        print ""
        print "RNA-seq analyses:"
        print "\tWindow width             : %d (bp)" % window_width
        print "\tBasename for output files: %s" % rna_basename
    if xls_out:
        print ""
        print "Outputting results to XLS file   : %s" % xls_out

    # Initialise the data objects
    try:
        rna_seq = RNASeqData(rnaseq_file)
    except Exception, ex:
        logging.critical("Failed to read in RNA-seq data: %s" % ex)
        print "Please fix errors in input file before running again"
        sys.exit(1)
    chip_seq = ChIPSeqData(chipseq_file)

    # Check we have data
    if not len(rna_seq):
        print "ERROR No RNA-seq data read in, stopping"
        sys.exit(1)
    else:
        print "%d RNA-seq records read in" % len(rna_seq)
        if rna_seq.isFlagged():
            print "RNA-seq data is flagged"
    if not len(chip_seq):
        print "ERROR No ChIP-seq data read in, stopping"
        sys.exit(1)
    else:
        print "%d ChIP-seq records read in" % len(chip_seq)
        print
        if chip_seq.isSummit():
            print "ChIP data appears to be peak summits, the following analyses will be run:"
            print "\tNearestTSSToSummits"
            print "\tNearestPeaksToTranscripts"
        else:
            print "ChIP data appears to be regions, the following analyses will be run:"
            print "\tNearestTranscriptsToPeakEdges"
            print "\tNearestTranscriptsToPeakEdges (TSS only)"
        print

    if xls_out:
        # Create initial XLS document
        xls = Spreadsheet.Workbook()
        xls_notes = xls.addSheet('Notes')
    else:
        xls = None

    # Analysis #1a: ChIP-seq perspective
    # NB this analysis disabled for now
    if False:
        AnalyseClosestTranscriptsToPeaksInEachDirection(chip_seq,rna_seq)

    # ChIP-seq-based analyses
    if do_chip_analyses:
        if chip_seq.isSummit():
            # "Nearest TSS to summit" analysis
            outfile = chip_basename+"_NearestTSSToSummits.txt"
            print "Running AnalyseNearestTSSToSummits"
            print "\tWriting output to %s" % outfile
            AnalyseNearestTSSToSummits(chip_seq,rna_seq,max_distance,
                                       max_closest,outfile,
                                       xls=xls)
            print "\tDone"
            if xls: xls_notes.addText(xls_notes_for_nearest_TSS_to_summits %
                                      max_distance)
        else:
            # "Nearest edge to peak region" analysis
            outfile = chip_basename+"_NearestTranscriptsToPeakEdges.txt"
            print "Running AnalyseNearestTranscriptsToPeakEdges (TSS/TES)"
            print "\tWriting output to %s" % outfile
            AnalyseNearestTranscriptsToPeakEdges(chip_seq,rna_seq,
                                                 promoter_region,
                                                 max_closest,
                                                 max_edge_distance,
                                                 TSS_only=False,
                                                 filename=outfile,
                                                 xls=xls)
            print "\tDone"
            if xls: xls_notes.addText(
                xls_notes_for_nearest_transcripts_to_peak_edges %
                (promoter_region,max_closest,max_edge_distance))

            # "Nearest TSS to peak region" analysis
            outfile = chip_basename+"_NearestTSSToPeakEdges.txt"
            print "Running AnalyseNearestTranscriptsToPeakEdges (TSS only)"
            print "\tWriting output to %s" % outfile
            AnalyseNearestTranscriptsToPeakEdges(chip_seq,rna_seq,
                                                 promoter_region,
                                                 max_closest,
                                                 max_edge_distance,
                                                 TSS_only=True,
                                                 filename=outfile,
                                                 xls=xls)
            print "\tDone"

    # RNA-seq-based analysese
    if do_rna_analyses:
        if chip_seq.isSummit():
            # "Nearest peak summits to TSS" analysis
            outfile = rna_basename+"_NearestPeaksToTranscripts.txt"
            print "Running AnalyseNearestPeaksToTranscripts"
            print "\tWriting output to %s" % outfile
            AnalyseNearestPeaksToTranscripts(rna_seq,chip_seq,window_width,
                                             filename=outfile,
                                             xls=xls)
            print "\tDone"
            if xls: xls_notes.addText(
                xls_notes_for_nearest_peaks_to_transcripts % window_width)
        else:
            # No analyses available for ChIP peak regions
            pass

    # Finish off spreadsheet output
    if xls:
        # Add the program version information to the spreadsheet
        xls_notes.addText("Produced by %s" % p.get_version())
        # Write the XLS file to disk
        xls.save(xls_out)

    # Finished
    print "Done"
    sys.exit()
                                      
