#!/bin/env python
#
#     output.py: functions for outputing analysis results
#     Copyright (C) University of Manchester 2015 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
output.py

Functions for outputing analysis results

"""
import distances
from Peaks import Peak
from analysis_redux import distance_closest_edge
from analysis_redux import distance_tss,distance_tes
import Spreadsheet

#######################################################################
# Constants
#######################################################################

MULTI_LINE=0
SINGLE_LINE=1

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
                 null_placeholder='.'):
        """
        Create new AnalysisReporter instance

        Arguments:
          promoter_region (tuple): promoter region extent (optional)
          mode (int): either SINGLE_LINE or MULTI_LINE
          fields (list): list of fields to output
          null_placeholder (str): placeholder to use in output for
            fields which evaluate to 'null'

        """
        self._fields = fields
        self._mode = mode
        self._promoter_region = promoter_region
        self._placeholder = null_placeholder
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
        nresults = len(results)
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
            return distance_closest_edge(peak,feature)
        elif attr == 'dist_TSS':
            return distance_tss(peak,feature)
        elif attr == 'dist_TES':
            return distance_tes(peak,feature)
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

    def make_header(self,max_hits):
        """
        Create a 'header' line for output

        Builds a header line which can be incorporated into
        an output file, based on the fields that are being
        reported.
        
        Arguments:
          max_hits (int): number of repeats to include for
            '*_list' based fields.

        Returns:
          str: tab-delimited header line

        """
        if self._mode == MULTI_LINE:
            return '\t'.join(self._fields)
        elif self._mode == SINGLE_LINE:
            header_fields = []
            for f in self._fields:
                try:
                    subfields = f[:-1].split('(')[1].split(',')
                    for i in range(1,max_hits+1):
                        for s in subfields:
                            header_fields.append("%s_%d" % (s,i))
                except IndexError:
                    header_fields.append(f)
            return '\t'.join(header_fields)


class XLS:
    """
    Class to assemble XLS output file

    Utility class to help build an XLS file from existing
    output TSV files.

    Example usage:

    >>> xls = XLS()
    >>> xls.add_result_sheet('results','results.tsv')
    >>> xls.write('results.xls')

    """
    def __init__(self):
        """
        Create a new XLS instance

        """
        self._xls = Spreadsheet.Workbook()
        self._char_limit = Spreadsheet.MAX_LEN_WORKSHEET_CELL_VALUE
        self._line_limit = Spreadsheet.MAX_NUMBER_ROWS_PER_WORKSHEET

    def add_result_sheet(self,title,tsv_file):
        """
        Add a sheet populated from a file

        Creates a new sheet in the spreadsheet with the
        supplied title and populates using the contents
        of a tab-delimited file.

        Arguments:
          title (str): a title for the sheet
          tsv_file (str): path to a tab-delimited file

        """
        ws = self._xls.addSheet(title)
        with open(tsv_file,'r') as fp:
            for i,line in enumerate(fp):
                if i == self._line_limit:
                    raise Exception("Too many rows for spreadsheet")
                ws.addText(line.rstrip('\n'))

    def write(self,xls_file):
        """
        Write XLS to a file

        Arguments:
          xls_file (str): name or path of output file

        """
        self._xls.save(xls_file)
