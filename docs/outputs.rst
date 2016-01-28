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
dist_TSS	 distance between peak and gene TSS
dist_TES	 distance between peak and gene TES
direction        'U' if hit is upstream, 'D' if downstream, '.' if
                 overlapped
overlap_gene	 1 if peak overlaps the gene, 0 if not
overlap_promoter 1 if peak overlaps the promoter region, 0 if not
================ ================================================

Each peak will appear as many times as there are nearest genes being
reported for that peak. For example::

    ...
    chr1  9619046  9619167  NM_178399_3110035E14Rik  +  9591248  9617222  ...
    chr1  9619046  9619167  NM_008651_Mybl1          -  9690280  9635825  ...
    chr1  9619046  9619167  NM_175236_Adhfe1         +  9538049  9570746  ...
    chr1  9619046  9619167  NM_021511_Rrs1           +  9535513	 9537532  ...
    ...

If there are no closest genes for a peak (based on the distance cutoff)
then the peak will still be reported but the remainder of the fields will
be filled with placeholders::

    chr17  23681171  23681172  ---  ---  ---  --- ...

Use the ``--compact`` option to output all the genes for each peak
on a single line (:ref:`compact_output`).

Peaks associated with each gene ('gene-centric' output)
-------------------------------------------------------

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
dist_TSS	distance between peak and gene TSS
direction       'U' if hit is upstream, 'D' if downstream, '.' if
                overlapped
in_the_gene     'YES' if peak overlaps the gene, 'NO' if not
=============== ====================================================

Each gene will appear as many times as there are nearest peaks being
reported for that gene. For example::

    ...
    BC021773_Glb1l  chr1  75193364  75207353  -  chr1  75481920  75482054  ...
    BC021773_Glb1l  chr1  75193364  75207353  -  chr1  75481920  75482054  ...
    ...

If there are no closest peaks to a gene (based on the distance cutoff)
then the gene will still be reported but the remainder of the fields
will be filled with placeholders::

    BC028767_3110009E18Rik  chr1  122017764  122114603  +  ---  ---  --- ...

Use the ``--compact`` option to output all the peaks for each genes
on a single line (see :ref:`compact_output`).

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
