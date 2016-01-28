#!/bin/env python
#
#     utils.py: utility functions for RnaChipIntegrator
#     Copyright (C) University of Manchester 2011-15 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
utils.py

Utility functions for RnaChipIntegrator:

- make_errline: highlight problem fields in a string
- truncate_text: truncate a text string to a specified length

"""

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
    items = line.rstrip().split('\t')
    for i in range(len(items)):
        if i in bad_fields:
            errline.append("^"*len(items[i]))
        else:
            errline.append(" "*len(items[i]))
    return '\t'.join(errline)

def truncate_text(text,max_len):
    """Truncate a text string

    Given a title and an optional extension, remove characters
    and replace with ellipsis (i.e. ...) so that it fit into
    the maxium number of characters (max_len).

    """
    len_text = len(text)
    if len_text <= max_len:
        return text
    text = text[len_text-max_len:]
    return '...' + text[3:]
