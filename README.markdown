RnaChipIntegrator: analyses of RNA-Seq and ChIP-Seq data
========================================================

RnaChipIntegrator is a utility that performs integrated analyses of RNA-Seq
and ChIP-Seq data, identifying the nearest ChIP peaks to each transcript, and
vice versa. The individual analyses differ from each other according to the
criteria that are used to calculate the distances between peaks and transcripts.

The program was originally written specifically for ChIP-Seq and RNA-Seq data but
can be used with any set of genomic features e.g. canonical genes, CpG islands etc,
or expression data e.g. microarrays.

Basic usage
-----------

The most basic form of usage is:

    RnaChipIntegrator <rna-data-file> <chip-data-file>

The RNA-seq data file must be a tab-delimited file with 5 columns of data for each
genomic feature (one per line):

    # ID  chr  start  end  strand

`chr` is the chromosome name, `start` and `end` define the limits of the feature,
and strand must be either `+` or `-`. `ID` is a name which is used to identify the
genomic feature in the output.

Optionally there can be a 6th column, indicating whether the gene was
differentially expressed (= 1) or not (= 0).

The ChIP-seq data file must be a tab-delimited file with 3 columns of data for each
ChIP peak (one per line):

    # chr  start  end

`chr` is the chromosome name (must match those in the RNA-seq file), and `start`
and `end` define the ChIP peaks - these can either be summits (in which case
'end' - 'start' = 1), or regions (with 'start' and 'end' indicating the extent).

Note that different analyses will be selected depending on whether the ChIP peaks
are defined as summits or regions.

Outputs
-------

The output of the program consists of one tab-delimited file for each analysis that
was performed, plus an XLS spreadsheet which has all the results plus a "notes" page
that explains the data from each analysis.

For full descriptions of each analysis see the `README` document in the `doc`
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