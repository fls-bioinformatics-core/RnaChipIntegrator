#!/bin/env python
#
#     xls_output.py: functions for writing analysis results to Excel files
#     Copyright (C) University of Manchester 2015-16 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
xls_output.py

Functions for outputting analysis results to XLSX spreadsheet

"""
import datetime
import xlsxwriter
import re
import output
import utils

# Regular expressions for styling tags
RE_STYLE = re.compile(r"^<style +([^>]*)>(.*)</style>$")

# Notes text
NOTES = dict()
NOTES['preamble'] = """<style font=bold bgcolor=gray>%s</style>

Find nearest peaks to %ss (and vice versa)

Bioinformatics Core Facility, Faculty of Life Sciences, University of Manchester
http://fls-bioinformatics-core.github.com/RnaChipIntegrator/
Run on %s

<style font=bold bgcolor=gray>Settings</style>"""
NOTES['peak_centric'] = """
<style font=bold bgcolor=gray>'Peak-centric': nearest %ss to each peak</style>
Column\tDescription"""
NOTES['feature_centric'] = """
<style font=bold bgcolor=gray>'%s-centric': nearest peaks to each %s</style>
Column\tDescription"""

class XLSX:
    """
    Class to assemble XLSX output file

    Utility class to help build an XLSX file from existing
    output TSV files.

    Example usage:

    >>> xlsx = XLS('results.xlsx')
    >>> xlsx.add_result_sheet('results','results.tsv')
    >>> xlsx.write()

    """
    def __init__(self,xlsx_file,program_version,feature_type=None):
        """
        Create a new XLSX instance

        Arguments:
          xlsx_file (str): name or path of output file
          program_version (str): name and version of the program
            that is writing the spreadsheet
          feature_type (str): if not None then replace 'feature'
            with 'feature_type' (e.g. 'gene', 'transcript' etc) in
            the output

        """
        self._xlsx = xlsxwriter.Workbook(xlsx_file)
        self._sheets = {}
        self._rows = {}
        self._widths = {}
        self._styles = {}
        self._feature_type = ('gene' if feature_type is None
                              else feature_type)
        self.add_sheet("Notes")
        self.append_to_notes(NOTES['preamble'] % (program_version,
                                                  self._feature_type,
                                                  datetime.date.today()))

    def add_sheet(self,name):
        """
        Create a new worksheet in the XLSX file

        Arguments:
          name (str): title for the new sheet (must be
            unique across the XLSX file)

        Returns:
          WorkSheet: new worksheet.

        """
        if name in self._sheets:
            raise KeyError("'%s': worksheet already exists")
        ws = self._xlsx.add_worksheet(name)
        self._sheets[name] = ws
        self._rows[name] = 0
        self._widths[name] = []
        return ws

    def get_format(self,*args):
        """
        Return a cell format object matching arguments

        Returns a Format object matching the supplied
        arguments, which should be strings of the form

        'KEY=VALUE'

        Formats are cached so there will be one Format
        per unique set of key/value pairs.

        """
        # Create a name for this style
        name = list(args)[:]
        name.sort()
        name = "_".join(name)
        # See if it's already defined
        if name not in self._styles:
            # Create a new style (cell_format)
            fmt = self._xlsx.add_format()
            for style in args:
                if style == "font=bold":
                    fmt.set_bold(True)
                elif style.startswith("bgcolor="):
                    color = style.split('=')[1]
                    fmt.set_bg_color(color)
                else:
                    raise NotImplementedError("%s: not implemented" %
                                              style)
            self._styles[name] = fmt
        # Return the cell format for this style
        return self._styles[name]

    def add_text(self,name,text):
        """
        Add (append) arbitrary text to a worksheet

        Arguments:
          name (str): name of the worksheet
          text (str): text that will be added to the
            end of the worksheet

        """
        ws = self._sheets[name]
        i = self._rows[name]
        for line in text.split('\n'):
            j = 0
            for item in line.split('\t'):
                # Check for styles
                style_match = RE_STYLE.match(item)
                if style_match:
                    item = style_match.group(2)
                    style = style_match.group(1).split()
                    fmt = self.get_format(*style)
                else:
                    fmt = None
                # Write the item
                ws.write(i,j,item,fmt)
                # Update the widths
                try:
                    self._widths[name][j] = max(self._widths[name][j],
                                                len(item))
                except IndexError:
                    self._widths[name].append(len(item))
                # Increment column counter
                j += 1
            # Increment row counter
            i += 1
        self._rows[name] = i

    def append_to_notes(self,text):
        """
        Append arbitrary text to the 'notes' page

        Arguments:
          text (str): text that will be added to the
            end of the notes.

        """
        self.add_text("Notes",text)

    def write_peak_centric(self,fields):
        """
        Write details of the 'peak-centric' results to XLSX notes

        Arguments:
          fields (list): list of fields in the output

        """
        self.append_to_notes(NOTES['peak_centric'] %
                             self._feature_type)
        self.append_to_notes(self._field_descriptions(fields))

    def write_feature_centric(self,fields):
        """
        Write details of the 'feature-centric' results to XLSX notes

        Arguments:
          fields (list): list of fields in the output

        """
        self.append_to_notes(NOTES['feature_centric'] %
                             (self._feature_type.title(),
                              self._feature_type))
        self.append_to_notes(self._field_descriptions(fields))

    def _field_descriptions(self,fields):
        """
        Generate field (column) descriptions for XLSX notes

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
        ws = self.add_sheet(title)
        # Get header line
        with open(tsv_file,'r') as fp:
            i = self._rows[title]
            for line in fp:
                j = 0
                for value in line.rstrip('\n').split('\t'):
                    ws.write(i,j,value)
                    try:
                        self._widths[title][j] = max(self._widths[title][j],
                                                     len(value))
                    except IndexError:
                        self._widths[title].append(len(value))
                    j += 1
                i += 1
            self._rows[title] = i
        # Freeze the header
        ws.freeze_panes(1,0)

    def write(self):
        """
        Write XLSX to file

        """
        # Set the column widths
        for name in self._sheets:
            ws = self._sheets[name]
            for j,w in enumerate(self._widths[name]):
                ws.set_column(j,j,w*1.2)
        # Close to write to file
        self._xlsx.close()
