#!/bin/env python
#
#     rearrange_columns.py: rearrange columns in tabular data files
#     Copyright (C) University of Manchester 2011 Peter Briggs
#
########################################################################
#
# rearrange_columns.py
#
#########################################################################

"""rearrange_columns.py

Utility program to allow manipulation of tab delimited data files.

Usage: python rearrange_columns.py <column_list> <input_file>

<column_list> is a comma-separated list of column indices (starting
from zero) specifying the order that columns from <input_file> should
be written to the output file. For example:

> python rearrange_columns.py 5,0,1,2,4 transcripts.bed

The output is written to stdout; redirect to a file name to create a
new file."""

#######################################################################
# Import modules that this module depends on
#######################################################################
import sys
import os

#######################################################################
# Main program
#######################################################################

if __name__ == "__main__":
    # Get the command line
    if len(sys.argv) != 3:
        print "Usage: %s <column_list> <input_file>" % \
            os.path.basename(sys.argv[0])
        sys.exit(1)
    # Deal with arguments
    #
    # Column list: comma-separated list of integers
    # specifying columns to copy from the <input_file>
    col_list = sys.argv[1].strip().split(',')
    #
    # Input file
    input_file = sys.argv[2]

    # Process input file line-by-line
    fp = open(input_file,'r')
    for linein in fp:
        lineout = []
        items = linein.strip().split('\t')
        for col in col_list:
            lineout.append(items[int(col)])
        print '\t'.join(lineout)
    # Reached end of file
    fp.close()
    sys.exit()
            
    
