.. _outputs:

Output files
============

Overview
--------

The default output of the program consists of a pair of tab-delimited
files:

* | **<BASENAME>_peak_centric.txt**
  | shows all of the genes associated with each peak ('peak-centric' analysis)

* | **<BASENAME>_gene_centric.txt**
  | shows all of the peaks associated with each gene ('gene-centric' analysis)

By default the output file ``BASENAME`` is taken from the name of the
input 'genes' file; use the ``--name`` option to set a custom basename.

Additional files may be produced depending on the options that have
been specified on the command line.

The format and content of each file is described in the following sections.

Genes associated with each peak ('peak-centric' output)
-------------------------------------------------------

By default the 'peak-centric' output file has one line for each
peak/gene pair that is being reported, with the following
columns of data for each:

================ ================================================
Name             Description
================ ================================================
peak.chr	 chromosome of the peak
peak.start	 peak start position
peak.end	 peak end position
gene.id	         gene ID
strand	         gene strand direction
TSS	         gene TSS position
TES	         gene TES position
dist_closest	 closest distance between peak and gene considering
                 all edges (zero if there is overlap)
dist_TSS	 distance between peak and feature TSS
dist_TES	 distance between peak and feature TES
direction        'U' if hit is upstream, 'D' if downstream, '.' if
                 overlapped
overlap_gene	 1 if peak overlaps the gene, 0 if not
overlap_promoter 1 if peak overlaps the promoter region, 0 if not
================ ================================================

Use the ``--compact`` option to output all the genes for each peak
on a single line (:ref:`compact_output`).

Peaks associated with each feature ('gene-centric' output)
----------------------------------------------------------

By default the 'gene-centric' file has one line for each
gene/peak pair that is being reported, with the following
columns of data for each:

=============== ====================================================
Name            Description
=============== ====================================================
gene.id	        gene ID
gene.chr	chromosome of the gene
gene.start	gene start position
gene.end	gene end position
gene.strand	gene strand direction
peak.chr	chromosome of the peak
peak.start	peak start position
peak.end	peak end position
dist_closest	closest distance between peak and gene considering
                all edges (zero if there is overlap)
dist_TSS	distance between peak and feature TSS
direction       'U' if hit is upstream, 'D' if downstream, '.' if
                overlapped
in_the_gene     'YES' if peak overlaps the feaure, 'NO' if not
=============== ====================================================

Use the ``--compact`` option to output all the peaks for each genes
on a single line (see :ref:`compact_output`).

.. _number:

Number of results to report (``--number``)
------------------------------------------

By default, all gene/peak pairs that are located within the
specified cut-off distance (see :ref:`distance_cutoff`) will be
reported in the output files.

To restrict the maximum number of pairs that are reported per gene
or peak use the ``--number`` to specify a limit. Even if more pairs
are found, only this number of pairs will be output.

.. warning::

   Be aware that if used, this number limit is applied rigidly.
   For example, even if the fourth and fifth gene/peak pairs both
   have the same distance separation then using ``--number=4``
   will only include the first of these and reject the second.

.. _summary_files:

Summary files (``--summary``)
-----------------------------

Using the ``--summary`` option outputs an additional pair of
tab-delimited files:

* ``BASENAME_peak_centric_summary.txt``
* ``BASENAME_gene_centric_summary.txt``

These will only contain the 'top' (i.e. closest) gene/peak pairs,
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

By default each gene/peak pair will be output on a separate line, for
example::

    #chr   start    end      gene.id     strand  TSS      TES      dist_closest dist_TSS dist_TES  overlap_gene  overlap_promoter
    chr2R  4959563  4959564  CG8084-RA   +       4956606  4965060  0            2957     5496      1             0
    chr2R  4959563  4959564  CG8193-RA   -       4932214  4929765  27349        27349    29798     0             0
    chr3R  12882217 12882218 CG3937-RB   -       12921260 12917257 35039        39042    35039     0             0
    ...

Specifying the ``--compact`` option changes the ouput so that all the
genes closest to each peak (and vice versa) are written on a single
line, for example::

    #chr   start    end      gene.id_1  gene.id_2  gene.id_3  gene.id_4
    chr2R  4959563  4959564  CG8084-RA  CG8193-RA
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

    #chr   start    end      gene.id     strand  TSS      TES      dist_closest dist_TSS dist_TES  overlap_gene     overlap_promoter
    chr2R  4959563  4959564  CG8084-RA   +       4956606  4965060  0            2957     5496      1                0
    chr2R  4959563  4959564  CG8193-RA   -       4932214  4929765  27349        27349    29798     0                0
    chr2R  4959563  4959564  ---         ---     ---      ---      ---          ---      ---       ---              ---
    chr2R  4959563  4959564  ---         ---     ---      ---      ---          ---      ---       ---              ---

.. _feature_type:

Specifying feature type other than 'gene' etc (``--feature``)
-------------------------------------------------------------

By default the program uses the term 'gene' in its outputs
regardless of the nature of the genomic features being examined.
This term can be changed to refer to a different feature type
by using the ``--feature`` option.

For example::

    --feature=transcript

in which case the word 'gene' will be replaced by 'transcript' in
output headers and so on.

.. note::

   The feature type is purely cosmetic and has no effect on the
   input or output file formats, or the analyses performed.
