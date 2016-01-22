.. _outputs:

Output files
============

Overview
--------

The default output of the program consists of a pair of tab-delimited
files:

* | **<BASENAME>_features_per_peak.txt**
  | shows all of the features associated with each peak

* | **<BASENAME>_peaks_per_feature.txt**
  | shows all of the peaks associated with each feature

By default the output file ``BASENAME`` is taken from the name of the
input feature file; use the ``--name`` option to set a custom basename.

Additional files may be produced depending on the options that have
been specified on the command line.

The format and content of each file is described in the following sections.

Features associated with each peak
----------------------------------

By default the 'features to peaks' file has one line for each
peak/feature pair that is being reported, with the following
columns of data for each:

================ ================================================
Name             Description
================ ================================================
peak.chr	 chromosome of the peak
peak.start	 peak start position
peak.end	 peak end position
feature.id	 feature ID
strand	         feature strand direction
TSS	         feature TSS position
TES	         feature TES position
dist_closest	 closest distance between peak and feature considering
                 all edges (zero if there is overlap)
dist_TSS	 distance between peak and feature TSS
dist_TES	 distance between peak and feature TES
direction        'U' if hit is upstream, 'D' if downstream, '.' if
                 overlapped
overlap_feature	 1 if peak overlaps the feature, 0 if not
overlap_promoter 1 if peak overlaps the promoter region, 0 if not
================ ================================================

Use the ``--compact`` option to output all the features for each peak
on a single line (:ref:`compact_output`).

Peaks associated with each feature
----------------------------------

By default the 'peaks to features' file has one line for each
feature/peak pair that is being reported, with the following
columns of data for each:

=============== ====================================================
Name            Description
=============== ====================================================
feature.id	feature ID
feature.chr	chromosome of the feature
feature.start	feature start position
feature.end	feature end position
feature.strand	feature strand direction
peak.chr	chromosome of the peak
peak.start	peak start position
peak.end	peak end position
dist_closest	closest distance between peak and feature considering
                all edges (zero if there is overlap)
dist_TSS	distance between peak and feature TSS
direction       'U' if hit is upstream, 'D' if downstream, '.' if
                overlapped
in_the_feature  'YES' if peak overlaps the feaure, 'NO' if not
=============== ====================================================

Use the ``--compact`` option to output all the peaks for each feature
on a single line (:ref:`compact_output`).

.. _summary_files:

Summary files (``--summary``)
-----------------------------

Using the ``--summary`` option outputs an additional pair of
tab-delimited files:

* ``BASENAME_features_per_peak_summary.txt``
* ``BASENAME_peaks_per_feature_summary.txt``

These will only contain the 'top' (i.e. closest) feature/peak pairs,
with the same columns of data as the 'full' versions of the files.

.. _xls_file:

Excel spreadsheet (``--xls``)
-----------------------------

Using the ``--xls`` option outputs an additional Excel spreadsheet
file ``BASENAME.xls``, which contains the results from all the
tab-delimited files (including the summaries, if ``--summary`` was
also specified), plus a 'notes' sheet with additional information
about the results from each analysis.

.. note::

   The XLS-writing library has a limit on the number of rows that
   can be written to a sheet; if the number of results exceeds this
   limit then the results will be broken into multiple sheets in
   the output XLS file.

.. _compact_output:

Compact output format (``--compact``)
-------------------------------------

By default each feature/peak pair will be output on a separate line, for
example::

    #chr   start    end      feature.id  strand  TSS      TES      dist_closest dist_TSS dist_TES  overlap_feature  overlap_promoter
    chr2R  4959563  4959564  CG8084-RA   +       4956606  4965060  0            2957     5496      1                0
    chr2R  4959563  4959564  CG8193-RA   -       4932214  4929765  27349        27349    29798     0                0
    chr3R  12882217 12882218 CG3937-RB   -       12921260 12917257 35039        39042    35039     0                0
    ...

Specifying the ``--compact`` option changes the ouput so that all the
features closest to each peak (and vice versa) are written on a
single line, for example::

    #chr   start    end      feature.id_1 feature.id_2 feature.id_3 feature.id_4
    chr2R  4959563  4959564  CG8084-RA    CG8193-RA
    chr3R  12882217 12882218 CG3937-RB

.. warning::

   ``--compact`` is not compatible with ``--summary``.

.. _output_padding:

Output padding (``--pad``)
--------------------------

If the ``--pad`` option is specified then where fewer than the
maximum number of pairs would be reported, additional 'blank'
lines are inserted to make up the number of lines to the maximum.

For example::

    #chr   start    end      feature.id  strand  TSS      TES      dist_closest dist_TSS dist_TES  overlap_feature  overlap_promoter
    chr2R  4959563  4959564  CG8084-RA   +       4956606  4965060  0            2957     5496      1                0
    chr2R  4959563  4959564  CG8193-RA   -       4932214  4929765  27349        27349    29798     0                0
    chr2R  4959563  4959564  ---         ---     ---      ---      ---          ---      ---       ---              ---
    chr2R  4959563  4959564  ---         ---     ---      ---      ---          ---      ---       ---              ---

.. _feature_type:

Feature type (``--feature``)
----------------------------

By default the program uses the generic term "feature" in its outputs
to describe the genomic features being examined.

A specific feature type can be specified using the ``--feature``
option, for example::

    --feature=gene

in which case the work "feature" will be replaced by "gene" in output
headers and so on.

.. note::

   The feature type is purely cosmetic and has no effect on the
   distance calculations.

