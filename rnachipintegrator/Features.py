#!/bin/env python
#
#     Features.py: classes for handling feature data
#     Copyright (C) University of Manchester 2011-15 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
Features.py

Classes for handling feature data.

"""

import logging
from distances import closestDistanceToRegion
from utils import make_errline

class FeatureSet:
    """Class for storing a set of features

    RNA-seq features consists of genes/transcripts/isomers, which
    are stored individually in Feature objects. This class is a
    container for a collection of Feature objects and provides
    methods to operate on the collection, by creating subsets by
    filtering, and sorting the features based on various criteria.

    """
    def __init__(self,features_file=None,features_list=None):
        """Create a new FeatureSet instance

        Raises an exception if there are errors in the input file data
        (non-numeric fields for start/end positions, end positions
        occurring before start positions, or illegal strand values).
        
        Arguments:
          features_file (str): (optional) the name of an input
            file to read the feature data from
          features_list (list): (optional) list of Feature objects
            to populate the FeatureSet with

        """
        self.features = []
        if features_file:
            self.loadFeaturesFromFile(features_file)
        elif features_list:
            for feature in features_list:
                self.addFeature(feature)

    def loadFeaturesFromFile(self,features_file):
        """Read features from a file and populate the object

        Arguments:
          features_file: the name of the input file to read features from.

        """
        # Local flags etc
        line_index = 0
        critical_error = False
        # Read in data from file
        fp = open(features_file,'rU')
        for line in fp:
            # Increment index
            line_index += 1
            # Skip lines starting with #
            if line.startswith('#'):
                logging.debug("Feature file: skipped line: %s" % line.strip())
                continue
            # Lines are tab-delimited and have at least 5 columns:
            # ID  chr  start  end  strand
            items = line.strip().split('\t')
            if len(items) < 5:
                logging.warning("Feature file: skipped line: %s" % line.strip())
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
                                    (features_file,line.strip()))
                else:
                    # Indicate problem field(s)
                    logging.error("%s: critical error line %d: bad values:" %
                                  (features_file,line_index))
                    logging.error("%s" % line.strip())
                    logging.error("%s" % make_errline(line.strip(),problem_fields))
                    # This is a critical error: update flag
                    critical_error = True
                # Continue to next line
                continue
            elif int(items[2]) >= int(items[3]):
                # Start position is same or higher than end
                logging.error("%s: critical error line %d: 'end' comes before 'start':" %
                              (features_file,line_index))
                logging.error("%s" % line.strip())
                logging.error("%s" % make_errline(line.strip(),(2,3)))
                # This is a critical error: update flag but continue reading
                critical_error = True
                continue
            # Store in a new Feature object
            feature = Feature(items[0],
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
                feature.flag = flag_value

            # Store data
            self.features.append(feature)
        fp.close()
        # Deal with postponed critical errors
        if critical_error:
            raise Exception, "critical error(s) in '%s'" % features_file
        # Return a reference to this object
        return self

    def addFeature(self,feature):
        """Append a feature to the FeatureSet object

        Arguments:
          feature: a Feature instance.

        """
        self.features.append(feature)

    def filterByChr(self,matchChr):
        """Return a subset of features filtered by specified chromosome name

        Returns a new FeatureSet object containing only the data from
        the current object which matches the specified criteria.

        """
        # Make a new (empty) FeatureSet object
        feature_subset = FeatureSet()
        # Populate with only the matching features
        for feature in self.features:
            if feature.chrom == matchChr:
                feature_subset.addFeature(feature)
        return feature_subset

    def filterByStrand(self,matchStrand):
        """Return a subset of features filtered by specified strand

        Returns a new FeatureSet object containing only the data from
        the current object which matches the specified criteria.

        """
        # Make a new (empty) FeatureSet object
        feature_subset = FeatureSet()
        # Populate with only the matching features
        for feature in self.features:
            if feature.strand == matchStrand:
                feature_subset.addFeature(feature)
        return feature_subset

    def filterByFlag(self,matchFlag):
        """Return a subset of features filtered by flag value

        Returns a new FeatureSet object containing only the features from
        the current object which matches the specified criteria.

        Note that if there is no flag (the "isFlagged()" function returns
        False) then an empty set will be returned.

        """
        # Make a new (empty) RNASeqData object
        feature_subset = FeatureSet()
        # Populate with only the matching features
        for feature in self.features:
            if feature.flag == matchFlag:
                feature_subset.addFeature(feature)
        return feature_subset

    def filterByTSS(self,limit1,limit2,exclude_limits=False):
        """Return a subset of features filtered by TSS position

        Returns a new FeatureSet object containing only the features
        from the current object where the TSS positions fall within a
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
        # Make a new (empty) FeatureSet object
        feature_subset = FeatureSet()
        # Populate with only the matching features
        for feature in self.features:
            TSS = feature.getTSS()
            if exclude_limits:
                if lower < TSS and TSS < upper:
                    feature_subset.addFeature(feature)
            else:
                if lower <= TSS and TSS <= upper:
                    feature_subset.addFeature(feature)
        return feature_subset

    def sortByDistanceFrom(self,position):
        """Sort the features into order based on distance from a position
    
        Sorts the features into order of absolute distance of
        their TSS to the specified position (closest first).
        
        Note that this operates on the current object.

        """
        self.features = sorted(self.features,
                               key=lambda record:
                               abs(record.getTSS()-position))
        return self

    def sortByClosestEdgeTo(self,position1,position2=None):
        """Sort the features into order based on closest edge (TSS or TES)

        Sorts the features into order of smallest absolute distance
        to the specified position (closest first), considering both TSS
        and TES, using the getClosestEdgeDistanceTo method of the
        Feature class.
        
        Note that this operates on the current object.

        """
        self.features = sorted(self.features,
                               key=lambda record: 
                               record.getClosestEdgeDistanceTo(position1,
                                                               position2))
        return self

    def sortByClosestTSSTo(self,position1,position2=None):
        """Sort the features into order based on closest edge to TSS

        Sorts the features into order of smallest absolute distance
        to the specified position (closest first) to the TSS position,
        using the getClosestTSSDistanceTo method of the Feature class.
        
        Note that this operates on the current object.

        """
        self.features = sorted(self.features,
                               key=lambda record:
                               record.getClosestTSSDistanceTo(position1,
                                                              position2))
        return self

    def isFlagged(self):
        """Check whether feature data includes flags

        Checks whether all the Feature records also have a valid flag
        associated with them - if yes then returns True (indicating the
        dataset as a whole is flagged), otherwise returns False.

        """
        # Check all data and look for any None flags
        for feature in self.features:
            if feature.flag is None:
                return False
        # All flags valid
        return True

    def __iter__(self):
        return iter(self.features)

    def __getitem__(self,key):
        try:
            start = key.start
            stop = key.stop
            step = key.step
            slice_ = FeatureSet()
            for feature in self.features[start:stop:step]:
                slice_.addFeature(feature)
            return slice_
        except AttributeError:
            return self.features[key]

    def __len__(self):
        return len(self.features)

    def __eq__(self,other):
        if len(self) != len(other):
            return False
        for f1,f2 in zip(self,other):
            if f1 != f2:
                return False
        return True

    def __ne__(self,other):
        if len(self) != len(other):
            return True
        for f1,f2 in zip(self,other):
            if f1 != f2:
                return True
        return False

class Feature:
    """Class for storing an 'feature' (gene/transcript/isomer)

    Access the data for the feature using the object's properties:

      id
      chrom
      start
      end
      strand
      tss
      tes

    There are also convenience methods (getTSS, getTES, getPromoterRegion)
    and methods for calculating various distances.

    """
    def __init__(self,feature_id,chrom,start,end,strand):
        self.id = feature_id
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.flag = None
        # Set the TSS and TES
        if self.strand == '+':
            self.tss = self.start
            self.tes = self.end
        elif self.strand == '-':
            self.tss = self.end
            self.tes = self.start
        else:
            raise Exception("Bad strand: '%s'" % self.strand)

    def __repr__(self):
        items = [self.id,
                 self.chrom,
                 str(self.start),
                 str(self.end),
                 self.strand]
        if self.flag != None:
            items.append(str(self.flag))
        return '\t'.join(items)

    def __eq__(self,other):
        return \
            (self.id == other.id) and \
            (self.strand == other.strand) and \
            (self.start == other.start) and \
            (self.end == other.end)

    def __ne__(self,other):
        return \
            (self.id != other.id) or \
            (self.strand != other.strand) or \
            (self.start != other.start) or \
            (self.end != other.end)

    def getTSS(self):
        """Return the TSS coordinate

        TTS (transcription start site) is the start position for a +ve 
        strand, or end for a -ve strand.

        This is a wrapper for accessing the 'tss' property.

        """
        return self.tss

    def getTES(self):
        """Return the TES coordinate

        TES (transcription end site) is the start position for a +ve
        strand, or end for a -ve strand.

        This is a wrapper for accessing the 'tes' property.

        """
        return self.tes

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
