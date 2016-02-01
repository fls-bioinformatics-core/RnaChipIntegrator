#!/bin/env python
#
#     xls_output.py: functions for writing analysis results to Excel files
#     Copyright (C) University of Manchester 2015 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
xls_output.py

Functions for outputting analysis results to XLS spreadsheet

"""
import datetime
import Spreadsheet
import output
import utils

NOTES = dict()
NOTES['preamble'] = """<style font=bold bgcolor=gray25>%s</style>

Find nearest peaks to %ss (and vice versa)

Bioinformatics Core Facility, Faculty of Life Sciences, University of Manchester
http://fls-bioinformatics-core.github.com/RnaChipIntegrator/
Run on %s

<style font=bold bgcolor=gray25>Settings</style>"""
NOTES['peak_centric'] = """
<style font=bold bgcolor=gray25>'Peak-centric': nearest %ss to each peak</style>
Column\tDescription"""
NOTES['feature_centric'] = """
<style font=bold bgcolor=gray25>'%s-centric': nearest peaks to each %s</style>
Column\tDescription"""

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
    def __init__(self,program_version,feature_type=None):
        """
        Create a new XLS instance

        Arguments:
          program_version (str): name and version of the program
            that is writing the spreadsheet
          feature_type (str): if not None then replace 'feature'
            with 'feature_type' (e.g. 'gene', 'transcript' etc) in
            the output

        """
        self._xls = Spreadsheet.Workbook()
        self._title_limit = Spreadsheet.MAX_LEN_WORKSHEET_TITLE
        self._char_limit = Spreadsheet.MAX_LEN_WORKSHEET_CELL_VALUE
        self._line_limit = Spreadsheet.MAX_NUMBER_ROWS_PER_WORKSHEET
        self._notes = self._xls.addSheet("Notes")
        self._feature_type = ('gene' if feature_type is None
                              else feature_type)
        self._notes.addText(NOTES['preamble'] % (program_version,
                                                 self._feature_type,
                                                 datetime.date.today()))

    def append_to_notes(self,text):
        """
        Append arbitrary text to the 'notes' page

        Arguments:
          text (str): text that will be added to the
            end of the notes.

        """
        self._notes.addText(text)

    def write_peak_centric(self,fields):
        """
        Write details of the 'peak-centric' results to XLS notes

        Arguments:
          fields (list): list of fields in the output

        """
        self.append_to_notes(NOTES['peak_centric'] %
                             self._feature_type)
        self.append_to_notes(self._field_descriptions(fields))

    def write_feature_centric(self,fields):
        """
        Write details of the 'feature-centric' results to XLS notes

        Arguments:
          fields (list): list of fields in the output

        """
        self.append_to_notes(NOTES['feature_centric'] %
                             (self._feature_type.title(),
                              self._feature_type))
        self.append_to_notes(self._field_descriptions(fields))

    def _field_descriptions(self,fields):
        """
        Generate field (column) descriptions for XLS notes

        Arguments:
          fields (list): list of fields to describe

        Returns:
          string: text with one field name/description pair
            (separated by a tab) per line

        """
        return '\n'.join(['\t'.join(x) for x in
                          output.describe_fields(fields)]).\
                             replace('feature',self._feature_type).\
                             replace('Feature',self._feature_type.title())

    def add_result_sheet(self,title,tsv_file):
        """
        Add a sheet populated from a file

        Creates a new sheet in the spreadsheet with the
        supplied title and populates using the contents
        of a tab-delimited file.

        If there are more lines than can be written to a
        single worksheet then creates additional sheets
        as required.

        Arguments:
          title (str): a title for the sheet
          tsv_file (str): path to a tab-delimited file

        """
        ws = self._xls.addSheet(title)
        n = 0
        # Get header line
        with open(tsv_file,'r') as fp:
            header = fp.readline().rstrip('\n')
        with open(tsv_file,'r') as fp:
            i = 0
            for line in fp:
                if not i%self._line_limit and i:
                    # Start a new sheet
                    n += 1
                    new_title = utils.truncate_text("%s(%d)" % (title,n),
                                                    self._title_limit)
                    print "Wrapping onto new sheet: %s" % new_title
                    ws = self._xls.addSheet(new_title)
                    ws.addText(header)
                    i += 1
                ws.addText(line.rstrip('\n'))
                i += 1

    def write(self,xls_file):
        """
        Write XLS to a file

        Arguments:
          xls_file (str): name or path of output file

        """
        self._xls.save(xls_file)
