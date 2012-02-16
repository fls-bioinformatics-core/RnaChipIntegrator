RnaChipIntegrator: integration of gene expression and ChIP data
============================================================

RnaChipIntegrator is a utility that performs integrated analyses of gene expression
and ChIP data, identifying the nearest ChIP peaks to each transcript, and vice versa. 
Output of RnaChipIntegrator differs according to the criteria that  are used to 
calculate the distances between ChIP peaks and transcripts.

The program was originally written specifically for ChIP-Seq and RNA-Seq data but works
equally well for ChIP-chip and microarray expression data, and can also be used to integrate
any set of genomic features (e.g. canonical genes, CpG islands) with expression data.

Basic usage
-----------

The most basic form of usage is:

    RnaChipIntegrator <rna-data-file> <chip-data-file>

The RNA-seq data file must be a tab-delimited file with 5 required columns of data for each
genomic feature (one per line):

    ID  chr  start  end  strand

`chr` is the chromosome name, `start` and `end` define the limits of the feature,
and strand must be either `+` or `-`. `ID` is a name which is used to identify the
genomic feature in the output.

Optionally there can be a 6th column, indicating whether the transcript was
differentially expressed (= 1) or not (= 0).

The ChIP-seq data file must be a tab-delimited file with 3 columns of data for each
ChIP peak (one per line):

    chr  start  end

`chr` is the chromosome name (must match those in the RNA-seq file), and `start`
and `end` define the ChIP peaks - these can either be the summit of the binding region 
(in which case 'end' - 'start' = 1), or the entire binding region (with 'start' and 
'end' indicating the extent).

Note that different analyses will be performed depending on whether the ChIP peaks
are defined as summits or regions.

Outputs
-------

The output of the program consists of one tab-delimited file for each analysis that
was performed, plus an XLS spreadsheet which has all the results plus a "notes" page
that explains the data from each analysis.

For peaks defined as regions these are:

 * TranscriptsToPeakEdges: reports the nearest transcripts (up to 4) with the smallest
   distance from either their TSS or TES to the nearest peak edge.
   (A separate "summary" reports only the top hit for each peak edge.)

 * TSSToPeakEdges: reports the nearest transcripts (up to 4) with the smallest distance
   from their TSS to the nearest peak edge.
   (A separate "summary" reports only the top hit for each peak edge.)

For peaks defined as summits:

 * TSSToSummits: reports the nearest transcripts (up to 4) with the smallest distance
   from the TSS to the nearest peak summit.

 * PeaksToTranscripts: reports the nearest peak summits (up to 4) with the smallest
   distance to either the TSS or TES of each transcript.

For full descriptions of each analysis see the `MANUAL` document in the `doc`
subdirectory.

Installation
------------

See the `INSTALL` document.

Examples
--------

Example data files can be found in the `examples` subdirectory, and which can be used
as input to the program for test or demonstration purposes; see the `README` file in
the same directory for more information.

Licensing
---------

This software is licensed under the Artistic License 2.0; see the `LICENSE` document.

More information
----------------

See the `README` document in the `doc` subdirectory for more detailed descriptions of
each analysis, the available options, and the outputs from the program.

--