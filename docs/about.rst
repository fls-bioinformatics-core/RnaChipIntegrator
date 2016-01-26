.. _about:

What it does
============

``RnaChipIntegrator`` was designed to integrate genes/transcripts from
expression analysis (RNA-seq, microarrays) with ChIP-seq binding regions,
however it is flexible enough to allow the comparison of any genome
coordinate based data sets.

``RnaChipIntegrator`` answers the questions:

 * "Which genes are close to each of my ChIP-seq regions?", and
 * "Which ChIP-seq regions are close to each of my genes?".

The first data set, called **'genes'**, is strand specific and the genome
coordinates correspond to the transcription start site (TSS) and
transcription end site (TES), depending on the strand.

.. note::

   For strand and genome coordinates:

   * The start coordinate of a gene on the forward or '+' strand relates
     to the TSS;
   * The start coordinate of a gene on the reverse or '-' strand relates
     to the TES.

This is primarily gene or transcript annotation for the whole genome.
However, other non-gene features, such as CpG islands, can be used.

The second data set we call **'peaks'** are strand non-specific, including
only the start and end coordinate. This is primarily the coordinates of
ChIP-seq binding regions (a.k.a. peaks).

See the :ref:`inputs` section for more information about the input file
formats.

Example use cases ('gene' versus 'peak') include:

 * RNA-seq expressed genes versus ChIP-seq binding regions
 * Microarray expressed genes versus ChIP-chip binding regions
 * Total gene annotation versus ChIP-seq binding regions
 * Gene promoters versus CpG island annotation
