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
import Spreadsheet

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
