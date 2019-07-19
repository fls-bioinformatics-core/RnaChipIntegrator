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
peak.id          (Optional) peak ID (if ``--peak_id`` option was
                 used)
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
direction        'U' if peak is upstream (5') of gene; 'D' if peak
                 is downstream (3') of gene; '.' if overlapping
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
peak.id         (Optional) peak ID (if ``--peak_id`` option was
                used)
peak.chr	chromosome of the peak
peak.start	peak start position
peak.end	peak end position
dist_closest	closest distance between peak and gene considering
                all edges (zero if there is overlap)
dist_TSS	distance between peak and gene TSS
direction       'U' if gene is upstream (5') of peak; 'D' if gene
                is downstream (3') of peak; '.' if overlapping
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

.. _xlsx_file:

Excel spreadsheet (``--xlsx``)
------------------------------

Using the ``--xlsx`` option outputs an additional Excel spreadsheet
file ``BASENAME.xlsx``, which contains the results from all the
tab-delimited files (including the summaries, if ``--summary`` was
also specified), plus a 'notes' sheet with additional information
about the results from each analysis.

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

.. _peak_id:

Specifying an ID for input peaks (``--peak_id``)
------------------------------------------------

If the ``--peak_id`` option is specified on the command line
then this indicates a column in the input peaks file which
should be used as names for each of the peaks in that file.

For example, if the input peaks file looks like::

    #Chrom	Start	End	Name
    chr1	9619046	9619167	P0001
    chr1	9619175	9619382	P0002
    chr1	10617233	10617437	P0003
    ...

then using::

    --peak_id=4

will associate the names ``P0001``, ``P0002``, ``P0003``...
with the corresponding peaks.

When specified, this ID is carried through to the output file
as an additional field, for example::

    ...
    P0001	chr1	9619046	9619167	NM_021511_Rrs1	+	9535513	9537532	81514	83533	81514	D	0	0
    P0002	chr1	9619175	9619382	NM_178399_3110035E14Rik	+	9591248	96172221953	27927	1953	D	0	0
    ...

.. _split_outputs:

Writing results to separate files in batch mode (``--split-outputs``)
---------------------------------------------------------------------

By default in 'batch' mode (i.e. when multiple cutoff distances and/or
multiple peak or genes files are supplied) all results for the
gene-centric analyses will be written to a single file (and similarly
for the peak-centric analyses).

To force ``RnaChipIntegrator`` to write the results of each batch to
a separate file, use the ``--split-outputs`` option. When this option
is specified a set of files will be generated for each peak, gene
and cutoff with appropriate names to indicate which files and cutoff
were used.

.. _additional_fields_for_batch_operation:

Additional fields for batch operation
-------------------------------------

When ``RnaChipIntegrator`` is run in 'batch' mode (that is, any mode
where multiple cutoffs have been supplied via the ``--cutoffs`` option,
and/or multiple input files have been supplied the ``--peaks`` and
``genes`` options), extra fields will be added for each reported
peak/gene pair to distinguish which analysis the result came from:

========== =====================================================
Name       Description
========== =====================================================
peak_file  Source file for the peak (if more than one peaks file
           was supplied via ``--peaks``)
gene_file  Source file for the gene (if more than one genes file
           was supplied via ``--genes``)
cutoff     maximum cutoff distance (if more than one cutoff
           distance was supplied via ``--cutoffs``)
========== =====================================================

For example::

    #peak_file       gene_file        cutoff  gene.id         gene.chr gene.start gene.end gene.strand ...
    /data/peaks1.txt /data/genes1.txt 50000   AF064749_Col6a3 chr1     92566771   92800755 -           ...

See the :ref:`multiple_distance_cutoffs` and :ref:`multiple_input_files`
sections for more information on these options.

.. note::

   Each of the additional fields will only appear if it is required
   in order to distinguish between the different analyses. For
   example, ``cutoff`` will only appear if more than one maxium
   cutoff distance was supplied.

.. _upstream_and_downstream:

Interpreting 'upstream' and 'downstream'
----------------------------------------

One of the attributes reported for each peak/gene pair found in the
analyses is the 'directionality' (in the ``direction`` column),
which can be either 'upstream' (``U``), 'downstream' (``D``) or
overlapped.

The intepretation of 'upstream' and 'downstream' for a given pairing
depends on the 'centricity' of the analysis and the strand direction.

For peak-centric analyses, the direction is from the point of view
of the peak::

                         ---Downstream-->    <---Upstream---

    + strand:  5' |----Gene1-------------Peak-------------Gene2----> 3'

In the example above, the peak is downstream of ``Gene1`` and upstream
of ``Gene2``.

(An analogy is that of a river which flows from the 5' to the 3' end;
the 'downstream' direction is the direction of flow from start to end,
while the 'upstream' direction is the opposite, from end to start.)

For the - strand this is reversed::

    - strand:  3' <----Gene3-------------Peak-------------Gene4----| 5'

                          ---Upstream--->    <---Downstream---

i.e. the peak is upstream of ``Gene3`` and downstream of ``Gene4``.

For gene-centric analyses, the direction is from the point of view
of the gene i.e. for the + strand::

                         ---Downstream-->    <---Upstream---

    + strand:  5' |----Peak1-------------Gene-------------Peak2----> 3'

(Here the gene is downstream of ``Peak1`` and upstream of ``Peak2``).

For the - strand::

    - strand:  3' <----Peak3-------------Gene-------------Peak4----| 5'

                          ---Upstream--->    <---Downstream---

(The gene is upstream of ``Peak3`` and downstream of ``Peak 4``).
