.. _usage:

Usage
=====

Simple usage
------------

The easiest form of usage is::

    RnaChipIntegrator GENES PEAKS

where ``GENES`` and ``PEAKS`` are tab-delimited files containing
the gene and peak data respectively (see :ref:`inputs` for details
of these files).

This will produce two output files:

 - ``GENES_peak_centric.txt``: reports the nearest genes
   for each peak ('peak-centric' analysis)
 - ``GENES_gene_centric.txt``: reports the nearest peaks
   for each gene ('gene-centric' analysis)

In both cases the files will contain one peak/gene pair per line
(see :ref:`outputs` for details of these files).

The program has various options that can be applied to control the
analyses that are performed and the outputs from each run, as outlined
in the following sections.

.. _distance_cutoff:

Specifying distance cutoff (``--cutoff``)
------------------------------------------

The ``--cutoff`` option specifies a maximum distance in bp that a
gene/peak pair can be apart and still be included in the analyses;
gene/peak pairs which are further apart than this distance will
not be reported.

For example::

    RnaChipIntegrator --cutoff=130000 GENES PEAKS

.. note::

   If a maximum cutoff distance is not explicitly specified then
   the default is 1000000 bp. Set the distance to 0 to turn off
   the cutoff limit and include all pairs regardless of distance.

Specifying how distances are measured between peaks and genes (``--edge``)
--------------------------------------------------------------------------

By default the distance between a peak and a gene is calculated
as the distance from the gene TSS to the nearest peak edge, for
example:

.. image:: nearest_tss.png
   :align: center

Alternatively distances can be calculated as the shortest distance
between either of the peak edges to either the TSS or the TES of
the gene, by specifying the ``--edge=both`` option::

    RnaChipIntegrator --edge=both GENES PEAKS

For example for the same arrangement as above this would generate a
much smaller closest distance:

.. image:: nearest_edge.png
   :align: center

.. note::

   Using ``--edge=both`` essentially makes the analyses
   'strand-agnostic'.

.. _using_differential_expression_data:

Only using differentially expressed genes (``--only-DE``)
---------------------------------------------------------

If the input genes data contains a differential expression flag
(see :ref:`genes_data_file`) then this can be used in the analysis
by turning on the ``--only-DE`` option::

    RnaChipIntegrator --only-DE GENES PEAKS

which will only included the flagged genes in the analyses.

.. note::

   Without the ``--only-DE`` option, all genes will be used
   regardless of the presence of a differential expression
   flag.

.. _number:

Limiting the number of results to report (``--number``)
-------------------------------------------------------

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

.. _promoter_region:

Specifying the promoter region (``--promoter_region``)
------------------------------------------------------

As part of its peak-centric analyses, for each peak/gene pair
``RnaChipIntegrator`` reports whether the peak overlaps the
promoter region of the gene.

By default, within the program the promoter region of a gene is
defined as starting 1000 bp upstream of the gene TSS and ending
100 bp downstream of the TSS.

The ``--promoter_region`` option can be used to define a different
set of limits for this region, using the general format::

    --promoter_region=UPSTREAM,DOWNSTREAM

For example::

    --promoter_region=1500,200

would define a promoter region starting 1500 bp upstream of the
TSS and ending 200 bp downstream.

.. _analyses_option:

Running either peak-centric or gene-centric analysis only (``--analyses``)
--------------------------------------------------------------------------

By default ``RnaChipIntegrator`` runs both peak-centric and
gene-centric analyses.

However it is possible to restrict the program to just one or
other of these, by using the ``--analyses`` option.

For example to run only the peak-centric analyses::

    --analyses=peak_centric

Or, to run only the gene-centric analyses::

    --analyses=gene_centric

The advantage of restricting the analyses is that it reduces the
program run time, and limits the outputs to only those specifically
requested.

Changing the output files and formats
-------------------------------------

There are a number of options to produce additional output files, and
to modify the format depending on requirements:

 * :ref:`xlsx_file`
 * :ref:`summary_files`
 * :ref:`compact_output`
 * :ref:`output_padding`
 * :ref:`feature_type`

Using RnaChipIntegrator in Galaxy
---------------------------------

In addition to the command-line version, we have also provided a tool
which allows ``RnaChipIntegrator`` to be run within the popular
`Galaxy <https://galaxyproject.org/>`_ bioinformatics platform:

 * https://toolshed.g2.bx.psu.edu/view/pjbriggs/rnachipintegrator/

The tool can be installed into a local instance of Galaxy directly from
the Galaxy Toolshed

See the documentation at http://getgalaxy.org/ on how to get a local
Galaxy up and running, and how to install tools from the Toolshed.
