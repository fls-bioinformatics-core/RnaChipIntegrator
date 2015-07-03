#!/usr/bin/env python
#
#     RnaChipIntegrator.py: analyse RNA-seq and ChIP-seq data
#     Copyright (C) University of Manchester 2011-14 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
#     This code is free software; you can redistribute it and/or modify it
#     under the terms of the Artistic License 2.0 (see the file LICENSE
#     included with the distribution).
#
########################################################################
#
# RnaChipIntegrator.py
#
#########################################################################

"""RnaChipIntegrator.py

Utility for integrated analyses of RNA-seq and ChIP-seq data.

Usage: RnaChipIntegrator.py [OPTIONS] <rna-data> <chip-data>

The RNA-seq data file should have at least 5 columns of data for each
gene or transcript:

# ID  chr  start  end  strand

Optionally there can be a 6th column, indicating whether the gene was
differentially expressed (= 1) or not (= 0).

The ChIP-seq data file should have 3 columns of data for each ChIP peak:

# chr  start  end

The data can be either peak summits (in which case 'end' - 'start' = 1),
or peak regions (with 'start' and 'end' indicating the extent).

The program performs analyses to locate the nearest peaks to each
transcript and vice versa using various criteria to define "nearest".
"""

#######################################################################
# Module metadata
#######################################################################

__version__ = "0.4.4"

#######################################################################
# Import modules that this module depends on
#######################################################################

import sys
import os
import optparse
import time,datetime

# Set default logging level and output
import logging
logging.getLogger().setLevel(logging.ERROR)
logging.basicConfig(format='%(levelname)s: %(message)s')

from rnachipintegrator.RNASeq import *
from rnachipintegrator.ChIPSeq import *
from rnachipintegrator.analysis import AnalyseNearestTranscriptsToPeakEdges
from rnachipintegrator.analysis import AnalyseNearestTSSToSummits
from rnachipintegrator.analysis import AnalyseNearestPeaksToTranscripts
import rnachipintegrator.Spreadsheet as Spreadsheet

#######################################################################
# Descriptions of each method for XLS notes
#######################################################################

xls_notes_preamble = \
"""<style font=bold bgcolor=gray25>%s: integrated analyses of RNA-Seq and ChIP-Seq data</style>

The following analyses have been performed and are reported in this spreadsheet.
"""

xls_notes_credits = \
"""<style font=bold bgcolor=gray25>Credits</style>
Produced by %s on %s
Bioinformatics Core Facility, Faculty of Life Sciences, University of Manchester 
http://fls-bioinformatics-core.github.com/RnaChipIntegrator/
"""

xls_notes_for_nearest_TSS_to_summits = \
"""<style font=bold bgcolor=gray25>TSS to Summits</style>
Find the nearest transcripts (up to 4) with the smallest distance from the TSS to the nearest peak summit.

<style font=bold>Input parameters:</style>
Cutoff distance from peaks\t%d bp

<style font=bold>Description of output fields:</style>
chr\tchromosome
summit\tposition of the peak summit
geneID\tgeneID for a closest differentially expressed gene/transcript
nearest\ta string of the form "1 of 4", "2 of 3" etc, indicating how many transcripts are listed for the peak, and which one of these the current transcript is.
TSS\tthe TSS position for the gene/transcript
distance_to_TSS\tdistance from the peak summit to the gene TSS
distance_to_TES\tdistance from the peak summit to the gene TES
strand\tthe strand direction
in_the_gene\tindicates whether the peak summit lies within the gene coordinates (either `YES` or `NO`)
transcripts_inbetween\tnumber of genes (differentially expressed or not) lying between the peak and the current gene
transcript_ids_inbetween\tlist of gene names lying between the peak and the current gene
"""

xls_notes_for_nearest_transcripts_to_peak_edges = \
"""<style font=bold bgcolor=gray25>Transcripts to Peak Edges</style>
Find the nearest transcripts (up to 4) with the smallest distance from either their TSS or TES to the nearest peak edge.

<style font=bold>Input parameters:</style>
Promoter region:\t%s bp
Maximum number of transcripts to report\t%d
Cutoff distance from peaks\t%d bp

<style font=bold>Description of output fields:</style>
chr\tchromosome
start\tpeak start position
end\tpeak end position
geneID\tgeneID for a closest gene/transcript
strand\tthe strand direction
TSS\tgene TSS position
TES\tgene TES position
dist_closest_edge\tclosest distance between the edges of the peak and gene regions.
dist_TSS\tdistance from the closest edge to the gene TSS.
dist_TES\tdistance from the closest edge to the gene TES.
overlap_transcript\tindicates whether the gene region overlaps the the peak region at any point (1 indicates an overlap, 0 no overlap).
overlap_promoter\tindicates whether the gene promoter region overlaps the peak region at any point (1 indicates an overlap, 0 no overlap).

<style font=bold bgcolor=gray25>Transcripts to Peak Edges (summary)</style>
Same as "Transcripts to Peak Edges" above but lists only the single nearest transcript to each peak.
"""

xls_notes_for_nearest_tss_to_peak_edges = \
"""<style font=bold bgcolor=gray25>TSS to Peak Edges</style>
Find the nearest transcripts (up to 4) with the smallest distance from their TSS to the nearest peak edge.
The input parameters and output fields are the same as for the "Transcripts to Peak Edges" analysis above.

<style font=bold bgcolor=gray25>TSS to Peak Edges (summary)</style>
Same as "TSS To Peak Edges" above but lists only the single nearest transcript to each peak.
"""

xls_notes_for_nearest_peaks_to_transcripts = \
"""<style font=bold bgcolor=gray25>Peaks to Transcripts</style>
Find the nearest peak summits (up to 4) with the smallest distance to either the TSS or TES of each transcript.

<style font=bold>Input parameters:</style>
Window width:\t%d bp

<style font=bold>Description of output fields:</style>
geneID\tgene/transcript ID
chr_RNA\tchromosome
start\tgene start position
end\tgene end position
strand\tthe strand direction
differentially_expressed\t1 indicates gene is differentially expressed, 0 indicates no significant
number_of_peaks\tthe number of ChIP peaks found within the "window" distance of the gene TSS

Then for each peak (closest first) there are three columns:
chr_ChIP_#\tchromosome (same as `chr_RNA` above)
summit_#\tpeak summit
distance_#\tdistance from the peak summit to the gene TSS
"""

#######################################################################
# Non-core Functions
#######################################################################

def count_unique_TSS(rna_seq):
    """Count the number of unique TSS positions

    Given an RNASeqData object with a list of transcripts, returns
    the number of unique TSS positions found in the list."""
    unique = []
    for rna_data in rna_seq:
        if not rna_data.getTSS() in unique:
            unique.append(rna_data.getTSS())
    return len(unique)

#######################################################################
# Main program
#######################################################################

def main():
    """Run the RnaChipIntegrator program

    This function should be invoked from __main__ to execute the program in
    command line mode.
    """
    # Initialisations
    do_chip_analyses = False
    do_rna_analyses = False
    xls_out = None
    max_distance = 130000
    max_edge_distance = 0
    window_width = 20000
    promoter_region = (10000,2500)
    max_closest = 4

    # Set default logging level
    logging.getLogger().setLevel(logging.INFO)
    p = optparse.OptionParser(usage="%prog [options] RNA-seq_data ChIP-seq_data",
                              version="%prog "+__version__,
                              description=
                              "Perform analyses of RNA-Seq data (or any set of genomic features "
                              "or expression data) with ChIP-Seq peaks, reporting nearest peaks "
                              "to each feature (and vice versa) according to various criteria "
                              "for calculating distances between them. ChIP-centric analyses "
                              "report the nearest features to each peak; RNA-seq-centric "
                              "analyses report the nearest peaks to each feature. "
                              "Input 'RNA-seq_data' file must contain tab-delimited columns "
                              "'ID,chr,start,end,strand[,flag]' (flag indicates differential "
                              "expression, either 1=yes or 0=no). "
                              "Input 'ChIP-seq_data' file must contain tab-delimited columns "
                              "'chr,start,stop' defining either summits (start/stop differ by 1 "
                              "bp) or regions (start/stop extend over several bps). "
                              "The outputs are: one tab-delimited file from each analysis "
                              "performed (named after the appropriate input file unless "
                              "overriden by the --project option), and an XLS spreadsheet with "
                              "one worksheet per analysis.")

    # General options
    p.add_option('--chip',action="store_true",dest="do_chip_analyses",
                     help="Do ChIP-seq-centric analyses")
    p.add_option('--rna',action="store_true",dest="do_rna_analyses",
                     help="Do RNA-seq-centric analyses")
    p.add_option('--project',action="store",dest="basename",
                     help="Set basename for output files; output from each "+
                     "analysis method will use this, with the method name appended"+
                     " (defaults to the input file names)")
    p.add_option('--no-xls',action="store_false",dest="write_xls_out",default=True,
                 help="Don't write an XLS file")
    p.add_option('--debug',action="store_true",dest="debug",
                     help="Verbose output for debugging")

    # Options for NearestPeaksToTranscripts
    group = optparse.OptionGroup(p,"NearestPeaksToTranscripts (RNA-seq)",
                                 description="For each transcript, reports the peaks with summit "+
                                 "positions that lie within the specified 'window' distance of "+
                                 "the transcript TSS.")
    group.add_option('--window',action="store",dest="window_width",
                     default=window_width,type='int',
                     help="Maximum distance a peak can be from each transcript TSS "+
                     "before being omitted from analysis "+
                     "(default %d bp)" % window_width)
    p.add_option_group(group)

    # Options for NearestTSSToSummits
    group = optparse.OptionGroup(p,"NearestTSSToSummits (ChIP-seq)",
                                 description="For each ChIP peak summit, reports the "+
                                 "transcripts with TSS positions that lie within the specified "+
                                 "cut-off distance of the peak summit.")
    group.add_option('--cutoff',action="store",dest="max_distance",
                     default=max_distance,type='int',
                     help="Maximum distance a transcript TSS can be from each peak "+
                     "before being omitted from the analysis "+
                     "(default %d bp)" % max_distance)
    p.add_option_group(group)

    # Options for NearestTranscriptsToPeakEdge/NearestTSSToPeakEdge
    group = optparse.OptionGroup(p,"NearestTranscriptsToPeakEdge/NearestTSSToPeakEdge (ChIP-seq)",
                                 description="For each ChIP peak, reports the transcripts that "+
                                 "lie closest to either 'edge' of the peak region, by "+
                                 "considering the TSS alone (NearestTSSToPeakEdge) or by "+
                                 "considering both the TSS and TES positions "+
                                 "(NearestTranscriptsToPeakEdge).")
    group.add_option('--edge-cutoff',action="store",dest="max_edge_distance",
                     default=max_edge_distance,type='int',
                     help="Maximum distance a transcript edge can be from the peak "+
                     "edge before being omitted from the analysis. Set to "+
                     "zero to indicate no cut-off (default %d bp)" % max_edge_distance)
    group.add_option('--number',action="store",dest="max_closest",
                     default=max_closest,type='int',
                     help="Maximum number of transcripts per peak to report from "+
                     "from the analysis (default %d)" % max_closest)
    group.add_option('--promoter_region',action="store",dest="promoter_region",
                     default="%d,%d" % promoter_region,
                     help="Define promoter region with respect to gene TSS "+
                     "(default -%d to %d bp of TSS)" %  promoter_region)
    group.add_option('--pad',action="store_true",dest="pad_output",
                 help="Where less than MAX_CLOSEST transcripts are found for a peak "
                 "add additional lines to the output to ensure that MAX_CLOSEST lines "
                 "are still reported")
    p.add_option_group(group)

    # Process the command line
    options,arguments = p.parse_args()

    # Input files
    if len(arguments) != 2:
        p.error("incorrect number of arguments")
    else:
        rnaseq_file = arguments[0]
        chipseq_file = arguments[1]

    # Report version and authors
    p.print_version()
    print "University of Manchester"
    print "Faculty of Life Sciences"
    print "Bioinformatics Core Facility"
    print "Authors: Ian Donaldson, Leo Zeef and Peter Briggs"
    print

    # Sort out analysis settings
    do_chip_analyses = options.do_chip_analyses
    do_rna_analyses = options.do_rna_analyses
    if not (do_chip_analyses or do_rna_analyses):
        # Neither explicitly requested - do both
        do_chip_analyses = True
        do_rna_analyses = True

    # Handle options
    max_distance = options.max_distance
    window_width = options.window_width
    max_edge_distance = options.max_edge_distance
    max_closest = options.max_closest
    write_xls_out = options.write_xls_out

    # Promoter region
    promoter_region = (abs(int(options.promoter_region.split(',')[0])),
                       abs(int(int(options.promoter_region.split(',')[1]))))

    # Output basename
    if options.basename:
        output_basename = options.basename
    else:
        output_basename = os.path.basename(os.path.splitext(rnaseq_file)[0]) + \
            "_vs_" + os.path.basename(os.path.splitext(chipseq_file)[0])
    if write_xls_out:
        xls_out = output_basename + ".xls"

    # Debugging output
    if options.debug: logging.getLogger().setLevel(logging.DEBUG)

    # Report settings
    print "Input transcripts file (RNA-seq) : %s" % rnaseq_file
    print "Input peaks file (ChIP-seq)      : %s" % chipseq_file
    if do_chip_analyses:
        print ""
        print "ChIP analyses:"
        print "\tMaximum cutoff distance   : %d (bp)" % max_distance
        print "\tMaximum edge distance     : %d (bp)" % max_edge_distance
        print "\tMax no. of closest genes  : %d" % max_closest
        print "\tPad output if fewer found : %s" % options.pad_output
        print "\tPromoter region           : -%d to %d (bp from TSS)" % \
            promoter_region
    if do_rna_analyses:
        print ""
        print "RNA-seq analyses:"
        print "\tWindow width             : %d (bp)" % window_width
    print
    print "Basename for output files        : %s" % output_basename
    if xls_out:
        print "Outputting results to XLS file   : %s" % xls_out
    print

    # Initialise the data objects
    try:
        rna_seq = RNASeqData(rnaseq_file)
    except Exception, ex:
        logging.critical("Failed to read in RNA-seq data: %s" % ex)
        print "Please fix errors in input file before running again"
        sys.exit(1)
    chip_seq = ChIPSeqData(chipseq_file)

    # Check we have data
    if not len(rna_seq):
        logging.error("No RNA-seq data read in")
        sys.exit(1)
    else:
        print "%d RNA-seq records read in" % len(rna_seq)
        if rna_seq.isFlagged():
            print "RNA-seq data is flagged"
            print "\t%d gene records flagged as differentially expressed" % \
                len(rna_seq.filterByFlag(matchFlag=1))
            print "\tOnly these will be used in the analyses"
    if not len(chip_seq):
        logging.error("No ChIP-seq data read in")
        sys.exit(1)
    else:
        print "%d ChIP-seq records read in" % len(chip_seq)
        print
        if chip_seq.isSummit():
            print "ChIP data appears to be peak summits, the following analyses will be run:"
            print "\tNearestTSSToSummits"
            print "\tNearestPeaksToTranscripts"
        else:
            print "ChIP data appears to be regions, the following analyses will be run:"
            print "\tNearestTranscriptsToPeakEdges"
            print "\tNearestTranscriptsToPeakEdges (TSS only)"
        print

    if xls_out:
        # Create initial XLS document
        xls = Spreadsheet.Workbook()
        xls_notes = xls.addSheet('Notes')
        xls_notes.addText(xls_notes_preamble % p.get_version())
    else:
        xls = None

    # Analysis #1a: ChIP-seq perspective
    # NB this analysis disabled for now
    if False:
        AnalyseClosestTranscriptsToPeaksInEachDirection(chip_seq,rna_seq)

    # ChIP-seq-based analyses
    if do_chip_analyses:
        if chip_seq.isSummit():
            # "Nearest TSS to summit" analysis
            outfile = output_basename+"_TSSToSummits.txt"
            print "Running AnalyseNearestTSSToSummits"
            print "\tWriting output to %s" % outfile
            AnalyseNearestTSSToSummits(chip_seq,rna_seq,max_distance,
                                       max_closest,outfile,
                                       xls=xls,pad_output=options.pad_output)
            print "\tDone"
            if xls: xls_notes.addText(xls_notes_for_nearest_TSS_to_summits %
                                      max_distance)
        else:
            # "Nearest edge to peak region" analysis
            outfile = output_basename+"_TranscriptsToPeakEdges.txt"
            print "Running AnalyseNearestTranscriptsToPeakEdges (TSS/TES)"
            print "\tWriting output to %s" % outfile
            AnalyseNearestTranscriptsToPeakEdges(chip_seq,rna_seq,
                                                 promoter_region,
                                                 max_closest,
                                                 max_edge_distance,
                                                 TSS_only=False,
                                                 filename=outfile,
                                                 xls=xls,
                                                 pad_output=options.pad_output)
            print "\tDone"
            if xls: xls_notes.addText(
                xls_notes_for_nearest_transcripts_to_peak_edges %
                (promoter_region,max_closest,max_edge_distance))

            # "Nearest TSS to peak region" analysis
            outfile = output_basename+"_TSSToPeakEdges.txt"
            print "Running AnalyseNearestTranscriptsToPeakEdges (TSS only)"
            print "\tWriting output to %s" % outfile
            AnalyseNearestTranscriptsToPeakEdges(chip_seq,rna_seq,
                                                 promoter_region,
                                                 max_closest,
                                                 max_edge_distance,
                                                 TSS_only=True,
                                                 filename=outfile,
                                                 xls=xls,
                                                 pad_output=options.pad_output)
            print "\tDone"
            if xls: xls_notes.addText(xls_notes_for_nearest_tss_to_peak_edges)

    # RNA-seq-based analysese
    if do_rna_analyses:
        if chip_seq.isSummit():
            # "Nearest peak summits to TSS" analysis
            outfile = output_basename+"_PeaksToTranscripts.txt"
            print "Running AnalyseNearestPeaksToTranscripts"
            print "\tWriting output to %s" % outfile
            AnalyseNearestPeaksToTranscripts(rna_seq,chip_seq,window_width,
                                             filename=outfile,
                                             xls=xls)
            print "\tDone"
            if xls: xls_notes.addText(
                xls_notes_for_nearest_peaks_to_transcripts % window_width)
        else:
            # No analyses available for ChIP peak regions
            pass

    # Finish off spreadsheet output
    if xls:
        # Add the program version information etc to the spreadsheet
        xls_notes.addText(xls_notes_credits % (p.get_version(),datetime.date.today()))
        # Write the XLS file to disk
        xls.save(xls_out)

    # Finished
    print "Done"
    sys.exit()

if __name__ == "__main__":
    # Run the program in command line mode
    main()
                                      
