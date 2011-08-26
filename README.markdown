RnaChipIntegrator.py: analyse RNA-seq and ChIP-seq data
=======================================================

Introduction
------------

This program implements a number of different analyses of RNA-seq data combined
with ChIP-seq data (and vice versa):

* **ChIP-seq data** consists of a number of peaks, which are described by the
  chromosome name plus either the summit (essentially a single point) or by
  the binding region (essentially the width of the peak as defined by upper
  and lower limits).

* **RNA-seq data** consists of genes, transcripts or isomers. These are
  described by the chromosome name, start and end coordinates and direction
  (\+ indicating the forward strand, \- indicating the reverse strand).

The analyses take one of two perspectives on the data:

1. **ChIP-seq perspective:** for each ChIP peak, locate the "nearest" genes
   in the RNA-seq data;

2. **RNA-seq perspective:** for each gene, locate the "nearest" peaks in the
   ChIP-seq data.

The analyses differ in how "nearest" is defined, the exact nature of the
ChIP peak input data (i.e. summit versus peak width), and the information that
is reported.

### Nearest genes to ChIP peak ###

The following analysis methods are available:

1.  "**NearestTranscriptsToPeaks**"

    For each ChIP peak, consider only the subset of genes that are flagged as
    "significant" and which lie on the same chromosome; from that subset find
    the genes with the closest _TSS position_ to the _peak summit_ within a
    user-defined cut-off distance (default 130 kb), regardless of whether they
    lie upstream or downstream.

    (See the section below on _Input Files_ for an explanation of the
    significance flag for input gene/transcript data.)

2.  "**NearestTranscriptsToPeakEdges**"

    For each ChIP peak, consider only the subset of genes which lie on the
    same chromosome; find the genes with the smallest distance from _either
    their TSS or TES_ to the nearest _peak edge_.

    (Note that in this analysis all genes are considered, regardless of
    whether or not they are flagged as "significant".)

3.  "**NearestTSSToPeakEdges**"

    For each ChIP peak, consider only the subset of genes which lie on the
    same chromosome; find the genes with the smallest distance from their
    _TSS position_ to the nearest _peak edge_.

    (Note that in this analysis all genes are considered, regardless of
    whether or not they are flagged as "significant".)

(A fourth method "_ClosestTranscriptsToPeaksInEachDirection_", which looks
for the single nearest transcript TSS to each peak summit on each strand
and in each upstream and downstream position - is only partially
implemented.)

### Nearest ChIP peaks to genes ###

The following analysis methods are available:

1.  "**NearestPeaksToTranscripts**"

    Matches peaks to genes. For each gene transcript flagged as "significant",
    find all ChIP peaks on the same chromosome for which the _peak summit_
    lies within a specified "window" around the gene _TSS position_
    (default is 20 kb upstream and downstream from the TSS).

    (See the section below on _Input Files_ for an explanation of the
    significance flag for input gene/transcript data.)

Program usage
-------------

To run the analyses:

    RnaChipIntegrator.py [OPTIONS] <rna-seq_data> <chip-seq_data>

The available command line options set parameters for different analyses:

_Options for NearestTranscriptsToPeaks_

*   `--cutoff=<max_distance>`: sets the maximum distance that transcript
    TSS values can be from the peak summit in the
    *NearestTranscriptsToPeaks* analysis; transcripts beyond this distance
    are not considered. Defaults to 130000 bp.
*   `--number=<number_of_transcripts>`: sets the maximum number of
    closest transcripts to report. Defaults to 4.

_Options for NearestTranscriptsToPeakEdges_, _NearestTSSToPeakEdges_

*   `--number=<number_of_transcripts>`: see above.
*   `--promoter_region=<leading>,<trailing>`: defines the promoter region
    around transcript TSS. _Leading_ defines the number of bases upstream of
    the TSS where the promoter region starts, _trailing_ sets the number
    downstream where it ends. Defaults to -10000 to 2500 bp of TSS.

_Options for NearestPeaksToTranscripts_

*   `--window=<window_width>`: sets the "window size" i.e. the maximum
    distance from the transcript TSS that ChIP peak summits will be counted.
    Defaults to 20000 bp.

_General options_

*   `--project=<name>`: set a basename to be used for output files from the
    analyses; see the section on _Output files_ below.
*   `--debug`: produce screen output for each analysis.

Input Files
-----------

All inputs are text files containing tab-delimited columns, one record per
line, in the following formats.

### ChIP peak data ###

ChIP peak data files should contain at least three columns. Only the first
three columns are read; additional columns are discarded.

The peaks can be described by either a _summit_:

>    `chr | summit | summit+1`

or by a _region_:

>    `chr | start | end`

### Gene transcipts ###

Gene transcript files should contain at least 5 columns; if more columns are
supplied then only the first 6 columns are processed and the rest are
discarded:

>    `ID | chr | start | end | strand | [flag]`

If there is a 6th column then the program attempts to process this as a
_significance flag_. The only valid values for significance flags are 1 and 0.
If any other values appear in this line then the dataset as a whole is
not considered to be flagged.

(Use the `rearrange_columns.py` program to put columns into the correct order;
see appendix A.)

`ID` is an arbitrary (to the analysis program) identifier for the gene or
transcript.

Output Files
------------

The following files are produced from each run; the "basename" is set by the
`--project` command line option, and defaults to the name of the program.

1.  **NearestTranscriptsToPeaks**: `<basename>_peaks_NearestTranscriptsToPeaks.txt`

    Each line has one ChIP peak matched to a closest gene transcript,
    with the following tab-delimited fields:

    `#chr | start | ID | nearest | TSS | distance_to_TSS | distance_to_TES | strand | in_the_gene | transcripts_inbetween | transcript_ids_inbetween`

    Descriptions of each of the fields:

    > `chr`: chromosome
    >
    > `start`: start position for the peak (assumed to be the summit)
    >
    > `ID`: ID for a closest gene/transcript
    >
    > `nearest`: a string of the form "1 of 4", "2 of 3" etc, indicating how
    > many transcripts are listed for the peak, and which one of these the
    > current transcript is.
    >
    > `TSS`: the TSS position for the gene/transcript
    >
    > `distance_to_TSS`: distance from the peak start (= summit) to the gene TSS
    >
    > `distance_to_TES`: distance from the peak start (= summit) to the gene TES
    >
    > `strand`: the strand direction
    >
    > `in_the_gene`: indicates whether the peak start position (= summit) lies
    > within the gene coordinates (either `YES` or `NO`).
    >
    > `transcripts_inbetween`: for analyses involving "flagged" genes, this is
    > the number of unflagged genes which lie closer to the peak position than
    > the current gene.
    >
    > `transcript_ids_inbetween`: for analyses involving "flagged" genes, this
    > is a semi-colon separated list of the unflagged gene names between the
    > peak and the current gene (e.g. `CG12178-RA;CG2674-RE;CG2674-RB`).

2.  **NearestPeaksToTranscripts**: `<basename>_transcripts_NearestPeaksToTranscripts.txt`

    Each line has one gene transcript along with the closest ChIP peak
    summits, with the following tab-delimited fields:

    `#ID | chr_RNA | start | end | strand | differentially_expressed | number_of_peaks | chr_ChIP_1 | summit_1 | distance_1 | ...`

    Descriptions of each of the fields:

    > `ID`: gene/transcript ID
    >
    > `chr_RNA`: chromosome
    >
    > `start`: gene start position
    >
    > `end`: gene end position
    >
    > `strand`: the strand direction
    >
    > `differentially_expressed`: the value of the "significance flag" supplied
    > on input
    >
    > `number_of_peaks`: the number of ChIP peaks found within the "window"
    > distance of the gene TSS.

    Then for each peak (closest first) there are three columns:

    > `chr_ChIP_#`: chromosome (same as `chr_RNA` above)
    >
    > `summit_#`: peak start (=summit)
    >
    > `distance_#`: distance from the peak summit to the gene TSS

3.  **NearestTranscriptsToPeakEdges**: `<basename>_peaks_NearestTranscriptsToPeakEdges.txt`

    Each line has one ChIP peak matched to a closest gene transcript,
    with the following tab-delimited fields:

    `#chr | start | end | ID | strand | TSS | TES | dist_closest_edge | dist_TSS | dist_TES | overlap_transcript | overlap_promoter`

    Descriptions of each of the fields:

    > `chr`: chromosome
    >
    > `start`: peak start position
    >
    > `end`: peak end position
    >
    > `ID`: ID for a closest gene/transcript
    >
    > `strand`: the strand direction
    >
    > `TSS`: gene TSS position
    >
    > `TES`: gene TES position
    >
    > `dist_closest_edge`: closest distance between the edges of the peak
    > and gene regions.
    >
    > `dist_TSS`: distance from the closest edge to the gene TSS.
    >
    > `dist_TES`: distance from the closest edge to the gene TES.
    >
    > `overlap_transcript`: indicates whether the gene region overlaps the
    > the peak region at any point (1 indicates an overlap, 0 no overlap).
    >
    > `overlap_promoter`: indicates whether the gene promoter region overlaps
    > the peak region at any point (1 indicates an overlap, 0 no overlap).
    > The extent of the promoter region is defined using the `--promoter_region`
    > command line option when running the program.

    This file contains groups of transcripts for each peak; a second file
    called `<basename>_peaks_NearestTranscriptsToPeakEdges_summary.txt`
    contains the data for just the single nearest transcript to each gene.

4.  **NearestTSSToPeakEdges**: `<basename>_peaks_NearestTSSToPeakEdges.txt`

    The output file format and fields are the same as for the
    _NearestTranscriptsToPeakEdges_ analysis, except:

    > `dist_closest_edge`: is the distance from closest edge to the gene TSS
    > (i.e. the same as `dist_TSS`).

    There is also a second file called
    `<basename>_peaks_NearestTSSToPeakEdges_summary.txt` which contains the
    data for just the single nearest transcript to each gene.

Issues
------

1.   Not all analyses might be appropriate/useful for all data (e.g. peak
     summit data versus peak regions).

2.   "Flagged" gene/transcript data is not handled consistently (i.e. the
     flag is accounted for in some analyses but not others).

3.   There is potential ambiguity in analyses which consider peak regions, if
     the peak is partially or wholly within the gene/transcript.

4.   (Internal) The outputs could be separated from the analyses.

5.   There is not currently an option for writing analyses to XLS spreadsheets.

6.   There is not currently an option to select a subset of analyses (i.e.
     all analyses always run).

Appendix A: `rearrange_columns.py`: manipulating tab-delimited data files
-------------------------------------------------------------------------

It's possible that your files are in the wrong format for input into the
analysis program. The `rearrange_columns.py` utility can be used to create
files with the correct arrangement of columns:

    rearrange_columns.py <column_list> <input_file>

`<column_list>` is a comma-separated list of column indices specifying the
order that columns from the `<input_file>` must be written in, e.g.

    python rearrange_columns.py 5,0,1,2,4 transcripts.bed

The output is directed to stdout; redirect to a file name to create a
new file.

Appendix B: using the Python classes and functions interactively
----------------------------------------------------------------

The analyses are built around a set of Python classes and functions can be
combined to filter and sort the data in a relatively straightforward manner.

The key classes are:

*    **RNASeqData**: a container for a set of genes or transcripts
*    **RNASeqDataLine**: a single gene or transcript
*    **ChIPSeqData**: a container for a set of ChIP peaks
*    **ChIPSeqDataLine**: a single ChIP peak

To use interactively, start Python from the command line and import the module:

    >>> from RnaChipIntegrator import *

Then load some genes and peaks:

    >>> genes = RNASeqData('/path/to/gene.bed')
    >>> peaks = ChIPSeqData('/path/to/peaks.tsv')

See how many genes were loaded:

    >>> len(genes)

Get a subset of genes on a specific chromosome:

    >>> gene_subset = genes.filterByChr('chr1') # Only chr1 in gene_subset

Filter again by strand, either:

    >>> gene_subset = gene_subset.filterByStrand('+')

or

    >>> gene_subset = genes.filterByChr('chr1').filterByStrand('+')

Order the subset by TSS distance from a point:

    >>> gene_subset.sortByClosestTSSDistanceTo(10000)

or from a region:

    >>> gene_subset.sortByClosestTSSDistanceTo(10000,110000)

There is more documentation available by running [`pydoc`](http://docs.python.org/library/pydoc.html) on the Python source, e.g.:

    pydoc RnaChipIntegrator

