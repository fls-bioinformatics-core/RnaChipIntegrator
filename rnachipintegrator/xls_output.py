#!/bin/env python
#
#     xls_output.py: functions for writing analysis results to Excel files
#     Copyright (C) University of Manchester 2015 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
xls_output.py

Functions for outputting analysis results

"""
import datetime
import Spreadsheet
import xls_notes

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
    def __init__(self,program_version):
        """
        Create a new XLS instance

        Arguments:
          program_version (str): name and version of the program
            that is writing the spreadsheet

        """
        self._xls = Spreadsheet.Workbook()
        self._char_limit = Spreadsheet.MAX_LEN_WORKSHEET_CELL_VALUE
        self._line_limit = Spreadsheet.MAX_NUMBER_ROWS_PER_WORKSHEET
        self._notes = self._xls.addSheet("Notes")
        self._notes.addText(xls_notes.preamble % (program_version,
                                                  datetime.date.today()))

    def append_to_notes(self,text):
        """
        Append arbitrary text to the 'notes' page

        Arguments:
          text (str): text that will be added to the
            end of the notes.

        """
        self._notes.addText(text)

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
