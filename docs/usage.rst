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
 - ``GENES_feature_centric.txt``: reports the nearest peaks
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
   the default is 1000000 bp.

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

Changing the output files and formats
-------------------------------------

There are a number of options to produce additional output files, and
to modify the format depending on requirements:

 * :ref:`number`
 * :ref:`xls_file`
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
