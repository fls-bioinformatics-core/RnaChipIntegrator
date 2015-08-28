#!/bin/env python
#
#     output.py: functions for outputting analysis results
#     Copyright (C) University of Manchester 2015 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
output.py

Functions for outputing analysis results

"""
import distances
from Peaks import Peak

#######################################################################
# Constants
#######################################################################

MULTI_LINE=0
SINGLE_LINE=1
FIELDS = {
    'chr': "chromosome",
    'start': "peak start position",
    'end': "peak end position",
    'id': "feature ID",
    'strand': "feature strand direction",
    'TSS': "feature TSS position",
    'TES': "feature TES position",
    'peak.chr': "chromosome of the peak",
    'peak.start': "peak start position",
    'peak.end': "peak end position",
    'feature.id': "feature ID",
    'feature.chr': "chromosome of the feature",
    'feature.start': "feature start position",
    'feature.end': "feature end position",
    'feature.TSS': "feature TSS position",
    'feature.TES': "feature TES position",
    'feature.strand': "feature strand direction",
    'dist_closest': "closest distance between peak and feature",
    'dist_TSS': "distance between peak and feature TSS",
    'dist_TES': "distance between peak and feature TES",
    'overlap_feature': "1 if peak overlaps the feature, 0 if not",
    'overlap_promoter': "1 if peak overlaps the promoter region, 0 if not",
    'in_the_gene': "'YES' if peak overlaps the feaure, 'NO' if not",
    'differentially_expressed': "1 if feature is differentially expressed, 0 if not",
    'order': "the 'order' of the feature/peak pair (e.g. '1 of 4')",
    'number_of_results': "number of hits being reported",
}

#######################################################################
# Classes
#######################################################################

class AnalysisReporter:
    """
    Class to handle reporting of analysis results

    Once initialised the reporter can be used to generate 'reports'
    of each peak along with the nearest features (using the
    'report_nearest_features' method) or for each feature along
    with the nearest peaks (using 'report_nearest_peaks').

    Output can be in either 'multi-line' (one line per result pair),
    or 'single-line' format (one line containing all results).

    For each method a list of the fields to be reported can be
    specified. Available fields are:

    - (peak.)chr: chromosome for the peak
    - (peak.)start: peak start position
    - (peak.)end: peak end position
    - (feature.)id: feature ID
    - feature.chr: chromosome for the feature
    - feature.start: feature start position
    - feature.end: feature end position
    - (feature.)TSS: feature TSS
    - (feature.)TES: feature TES
    - (feature.)strand: feature strand
    - dist_closest: closest distance between peak and feature
    - dist_TSS: distance between peak and feature TSS
    - dist_TES: distance between peak and feature TES
    - overlap_feature: 'YES' if peak overlaps the feaure, 'NO' if not
    - overlap_promoter: 1 if peak overlaps the promoter region, 0 if not
    - in_the_gene: synonym for 'overlap_feature'
    - 'differentially_expressed': flag value for feature

    (In the field names above, the parts in (...) are optional e.g.
    'chr' == 'peak.chr' etc.)

    For multi-line output these additional fields are available:

    - order: the 'order' of the feature/peak pair (e.g. '1 of 4')

    For single-line output these additional fields are available:

    - number_of_results
    - list(...): output all results (peaks or features, as
      appropriate)

    For the 'list' options, the paranthese should enclose a list
    of fields to output for each peak or feature in the list e.g.
    'list(chr,start,dist_closest)' or 'list(feature.id)'.

    The following fields have not been implemented:

    - features_inbetween

    """
    def __init__(self,mode,fields,promoter_region=None,
                 max_hits=None,pad=False,
                 null_placeholder='.'):
        """
        Create new AnalysisReporter instance

        Arguments:
          mode (int): either SINGLE_LINE or MULTI_LINE
          fields (list): list of fields to output
          promoter_region (tuple): promoter region extent (optional)
          max_hits (int): optional maximum number of hits to
            report for each set of results
          null_placeholder (str): placeholder to use in output for
            fields which evaluate to 'null'
          pad (bool): add extra 'None' items to output hits to
            pad out to max_closest results

        """
        self._fields = fields
        self._mode = mode
        self._promoter_region = promoter_region
        self._placeholder = null_placeholder
        self._max_hits = max_hits
        self._pad = pad
        self._context_peak = None
        self._context_feature = None

    def report_nearest(self,reference,results):
        """
        Return details of nearest objects to a reference

        This is a generic reporting method which can handle
        either nearest features to a reference peak (in which
        case ``reference`` should be a Peak and ``results``
        the corresponding FeatureSet), or nearest peaks to a
        reference Feature (when ``reference`` is a Feature and
        ``results`` is a PeakSet).
        
        Arguments:
          reference (Object): reference object (i.e.
            Peak or Feature of interest)
          results (Object): list of corresponding results
            i.e. FeatureSet (for reference Peak) or
            PeakSet (reference Feature)

        Yields:
          string: line(s) of text reporting the results

        """
        # Initialise and set the context
        if isinstance(reference,Peak):
            self._context_peak = reference
            is_features = True
        else:
            self._context_feature = reference
            is_features = False
        # Reduce to maximum number of hits
        if self._max_hits is not None:
            results = results[:self._max_hits]
        else:
            results = results[:]
        nresults = len(results)
        # Pad with null results
        if self._pad:
            while len(results) < self._max_hits:
                if is_features:
                    results.addFeature(None)
                else:
                    results.addPeak(None)
        # Write the results
        if self._mode == SINGLE_LINE:
            # Report everything on a single line
            line = []
            for field in self._fields:
                if field == 'number_of_results':
                    value = nresults
                elif field.startswith('list('):
                    # Extract the subfields
                    subfields = field[:-1].split('(')[1].split(',')
                    # Report list of features
                    value = []
                    for result in results:
                        if is_features:
                            self._context_feature = result
                        else:
                             self._context_peak = result
                        for subfield in subfields:
                            value.append(self.value_for(subfield))
                    value = '\t'.join([str(x) for x in value])
                else:
                    # All other fields
                    value = self.value_for(field)
                line.append(str(value))
            # Return (yield) the line
            yield '\t'.join(line)
        elif self._mode == MULTI_LINE:
            # Report each result pair on a new line
            i = 0
            for result in results:
                if is_features:
                    self._context_feature = result
                else:
                    self._context_peak = result
                i += 1
                line = []
                for field in self._fields:
                    if field == 'order' and result is not None:
                        value = '%d of %d' % (i,nresults)
                    else:
                        value = self.value_for(field)
                    line.append(str(value))
                # Return (yield) the line
                yield '\t'.join(line)
        # Reset the context
        self._context_peak = None
        self._context_feature = None

    def report_nearest_features(self,peak,features):
        """
        Return details of nearest features for a peak

        This is a wrapper for ``report_nearest``.

        Arguments:
          peak (Peak): peak of interest
          features (FeatureSet): list of nearest features

        Yields:
          string: line(s) of text reporting the results

        """
        for line in self.report_nearest(peak,features):
            yield line

    def report_nearest_peaks(self,feature,peaks):
        """
        Return details of nearest peaks for a feature

        This is a wrapper for ``report_nearest``.
        
        Arguments:
          feature (Feature): feature of interest
          peaks (PeakSet): list of nearest peaks

        Returns:
          string: block of text reporting the results

        """
        for line in self.report_nearest(feature,peaks):
            yield line

    def value_for(self,attr):
        """
        Return the value for the specified attribute

        Wraps '_value_for' method, and returns the null
        placeholder value in the event of an AttributeError
        being raised.

        Arguments:
          attr (string): attribute name

        Returns:
          Value of the field for the current peak/feature
          pair

        """
        try:
            return self._value_for(attr)
        except (AttributeError,KeyError):
            return self._placeholder

    def _value_for(self,attr):
        """
        Return the value for the specified attribute

        Given the name of a field/attribute (see above for
        a list and definition of each), return the value
        for the current peak/feature pair (which should have
        been set by the calling method in the '_context_peak'
        and '_context_feature' properties).

        Arguments:
          attr (string): attribute name

        Returns:
          Value of the field for the current peak/feature
          pair

        Raises:
          AttributeError: if valid ``attr`` cannot be derived
          KeyError: if ``attr`` is not a recognised attribute
            name
        
        """
        peak = self._context_peak
        feature = self._context_feature
        if attr == 'chr' or attr == 'peak.chr':
            return peak.chrom
        elif attr == 'peak.start' or attr == 'start':
            return peak.start
        elif attr == 'peak.end' or attr == 'end':
            return peak.end
        elif attr == 'id' or attr == 'feature.id':
            return feature.id
        elif attr == 'feature.chr':
            return feature.chrom
        elif attr == 'feature.start':
            return feature.start
        elif attr == 'feature.end':
            return feature.end
        elif attr == 'TSS':
            return feature.tss
        elif attr == 'TES':
            return feature.tes
        elif attr == 'strand' or attr == 'feature.strand':
            return feature.strand
        elif attr == 'differentially_expressed':
            return feature.flag
        elif attr == 'dist_closest':
            return distances.distance_closest_edge(peak,feature)
        elif attr == 'dist_TSS':
            return distances.distance_tss(peak,feature)
        elif attr == 'dist_TES':
            return distances.distance_tes(peak,feature)
        elif attr == 'overlap_feature' or attr == 'in_the_gene':
            if distances.regions_overlap((peak.start,peak.end),
                                         (feature.tss,feature.tes)):
                overlap_feature = 1
            else:
                overlap_feature = 0
            if attr == 'in_the_gene':
                overlap_feature = ('YES' if overlap_feature == 1 else 'NO')
            return overlap_feature
        elif attr == 'overlap_promoter':
            if self._promoter_region is not None:
                promoter = feature.getPromoterRegion(*self._promoter_region)
                if distances.regions_overlap((peak.start,peak.end),
                                             promoter):
                    overlap_promoter = 1
                else:
                    overlap_promoter = 0
            else:
                raise Exception("'overlap_promoter' requested but no "
                                "promoter region has been defined")
            return overlap_promoter
        elif attr == 'features_inbetween':
            raise NotImplementedError("'features_inbetween' not implemented")
        else:
            raise KeyError("Unrecognised report field: '%s'" % attr)

    def make_header(self):
        """
        Create a 'header' line for output

        Builds a header line which can be incorporated into
        an output file, based on the fields that are being
        reported.

        Returns:
          str: tab-delimited header line

        """
        if self._mode == MULTI_LINE:
            return '\t'.join(self._fields)
        elif self._mode == SINGLE_LINE:
            header_fields = []
            for f in self._fields:
                try:
                    # Handle fields in a list(...)
                    subfields = f[:-1].split('(')[1].split(',')
                    if self._max_hits is None:
                        continue
                    for i in range(1,self._max_hits+1):
                        for s in subfields:
                            header_fields.append("%s_%d" % (s,i))
                except IndexError:
                    # Not a list
                    header_fields.append(f)
            return '\t'.join(header_fields)

class AnalysisReportWriter(AnalysisReporter):
    """
    Write analysis results to file

    Wrapper for AnalysisReporter that writes the results to
    a file; optionally it can also write 'summary' files
    (only top result is reported).

    """
    def __init__(self,mode,fields,promoter_region=None,
                 null_placeholder='.',max_hits=None,pad=None,
                 outfile=None,summary=None):
        """
        Create new AnalysisReportWriter instance

        Arguments:
          mode (int): either SINGLE_LINE or MULTI_LINE
          fields (list): list of fields to output
          promoter_region (tuple): promoter region extent (optional)
          null_placeholder (str): placeholder to use in output for
            fields which evaluate to 'null'
          max_hits (int): optional maximum number of hits to
            report for each set of results
          pad (bool): add extra 'None' items to output hits to
            pad out to max_closest results
          outfile (str): name of output file to write results to
          summary (str): optional, name of file to write summary
            results to

        """
        AnalysisReporter.__init__(self,mode,fields,
                                  promoter_region=promoter_region,
                                  null_placeholder=null_placeholder,
                                  pad=pad,max_hits=max_hits)
        if outfile is not None:
            # Open output file and write header
            self._fp = open(outfile,'w')
            self._fp.write("#%s\n" % self.make_header())
        else:
            self._fp = None
        if summary is not None:
            # Open summary file and write header
            self._summary = open(summary,'w')
            self._summary.write("#%s\n" % self.make_header())
        else:
            self._summary = None

    def write_nearest_features(self,peak,features):
        """
        Write a set of features to the output file(s)

        Arguments:
          peak (Peak): peak of interest
          features (FeatureSet): list of nearest features

        """
        lines = list(self.report_nearest_features(peak,features))
        if self._fp is not None:
            self._fp.write("%s\n" % '\n'.join(lines))
        if self._summary is not None:
            self._summary.write("%s\n" % lines[0])

    def write_nearest_peaks(self,feature,peaks):
        """
        Write a set of peaks to the output file(s)

        Arguments:
          feature (Feature): feature of interest
          peaks (PeakSet): list of nearest peaks

        """
        lines = list(self.report_nearest_peaks(feature,peaks))
        if self._fp is not None:
            self._fp.write("%s\n" % '\n'.join(lines))
        if self._summary is not None:
            self._summary.write("%s\n" % lines[0])

    def close(self):
        """
        Close the files associated with the writer

        """
        try:
            self._fp.close()
        except AttributeError:
            pass
        try:
            self._summary.close()
        except AttributeError:
            pass

def describe_fields(fields):
    """
    Return list of field descriptions

    Creates a list consisting of (FIELD,DESC) tuples
    where FIELD is the name of the field and DESC is
    its corresponding description text.

    For example if the supplied fields were:

    "chr,start,id,dist_closest"

    then the resulting desciption list would look like:

    [('chr','Chromosome'),
     ('start','Peak start position'),
     ('id','Feature ID'),
     ('dist_closest','Closest distance between peak and feature')]

    Arguments:
      fields (list): list of fields

    Returns:
      list: list of (field,description) tuples.

    """
    descriptions = []
    for attr in fields:
        try:
            descriptions.append(("%s" % attr,
                                 "%s" % FIELDS[attr]))
        except KeyError:
            if attr.startswith('list('):
                descriptions.append(('For each hit:',))
                sub_attrs = attr[:-1].split('(')[1].split(',')
                for sub_attr in sub_attrs:
                    descriptions.append(("%s_#" % sub_attr,
                                         "%s" % FIELDS[sub_attr]))
    return descriptions
