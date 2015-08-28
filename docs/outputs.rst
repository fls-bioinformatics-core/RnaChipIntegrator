Output files
============

The output of the program consists of a pair of tab-delimited files:

* ``BASENAME_features_to_peak.txt``
* ``BASENAME_peaks_to_features.txt``

By default the output file ``BASENAME`` is taken from the name of the
input feature file; use the ``--name`` option to set a custom basename.

Additional files may be produced depending on the options that have
been specified on the command line:

* ``BASENAME_features_to_peak_summary.txt`` and 
  ``BASENAME_peaks_to_features_summary.txt``: tab-delimited files which
  contain only the top (i.e. closest) feature/peak pairs (if
  ``--summary`` is specified).
* ``BASENAME.xls``: Excel spreadsheet containing the results from all
  the tab-delimited files, plus a 'notes' sheet with additional
  information about the results from each analysis.

Compact output
**************

By default each feature/peak pair will be output on a separate line, for
example::

    #chr	start	end	feature.id	strand	TSS	TES	dist_closest	dist_TSS	dist_TES	overlap_feature	overlap_promoter
    chr2R	4959563	4959564	CG8084-RA	+	4956606	4965060	0	2957	5496	1	0
    chr2R	4959563	4959564	CG8193-RA	-	4932214	4929765	27349	2734	929798	0	0
    chr3R	12882217	12882218	CG3937-RB	-	12921260	12917257	35039	39042	35039	0	0
    ...

Specifying the ``--compact`` option changes the ouput so that all the
features closest to each peak (and vice versa) are written on a
single line, for example::

    #chr	start	end	feature.id_1	feature.id_2	feature.id_3	feature.id_4
    chr2R	4959563	4959564	CG8084-RA	CG8193-RA
    chr3R	12882217	12882218	CG3937-RB

.. warning::

   ``--compact`` is not compatible with ``--summary``.

.. note::

   In this mode only a "minimal" set of fields are reported.

Output padding
**************

If the ``--pad`` option is specified then where fewer than the
maximum number of pairs would be reported, additional 'blank'
lines are inserted to make up the number of lines to the maximum.

For example::

    #chr	start	end	feature.id	strand	TSS	TES	dist_closest	dist_TSS	dist_TES	overlap_feature	overlap_promoter
    chr2R	4959563	4959564	CG8084-RA	+	4956606	4965060	0	2957	5496	1	0
    chr2R	4959563	4959564	CG8193-RA	-	4932214	4929765	27349	27349	29798	0	0
    chr2R	4959563	4959564	---	---	---	---	---	---	---	---	---
    chr2R	4959563	4959564	---	---	---	---	---	---	---	---	---

