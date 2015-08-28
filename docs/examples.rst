Examples
========

The ``examples`` subdirectory in the source distribution includes
the following example input files:

Genomic feature data:

* ``ExpressionData.txt``: sample expression data including differential
  expression flags

Peaks:

* ``ChIP_summits.txt``: sample ChIP peak 'summit' data (i.e. ``start``
  and ``end`` separated by 1bp)
* ``ChIP_regions.txt``: sample ChIP peak 'region' data (i.e. ``start``
  and ``end`` wider than 1bp)

The simplest example run::

    RnaChipIntegrator.py --cutoff=130000 ExpressionData.txt ChIP_regions.txt
    RnaChipIntegrator.py --cutoff=130000 ExpressionData.txt ChIP_summits.txt

which will find the nearest peaks and features using the feature TSS
position only.

To analyse the region data considering both the TSS and TES positions::

    RnaChipIntegrator.py --cutoff=130000 --edge=both ExpressionData.txt ChIP_regions.txt
