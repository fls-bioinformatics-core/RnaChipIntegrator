#!/bin/env python
#
#     analysis.py: functions for performing different analyses
#     Copyright (C) University of Manchester 2011-15 Peter Briggs, Leo Zeef
#     & Ian Donaldson
#
"""
analysis.py

Functions for performing different analyses between RNA-seq and ChIP-seq
datasets:

- AnalyseClosestTranscriptsToPeaksInEachDirection
- AnalyseNearestTSSToSummits
- AnalyseNearestPeaksToTranscripts
- AnalyseNearestTranscriptsToPeakEdges

Also includes a class to hold the results from the analysis:

- AnalysisResult

and a utility function to write results to a tab-delimited file.

"""
import os
import logging
from distances import regions_overlap
from ChIPSeq import PeakSet
from utils import truncate_text
try:
    import Spreadsheet
except ImportError:
    pass

#######################################################################
# Class definitions
#######################################################################

class AnalysisResult:
    """Class to hold results from RNA-seq/ChIP-seq analysis

    This class is a general container for results from an analysis.

    Each entry is a Python dictionary which associates a ChIPSeqDataLine
    and a Feature pair along with any additional arbitrary data,
    via the addResult method, e.g.

    >>> result = AnalysisResult()
    >>> result.addResult(chip,features,distance_to_TSS=12345,in_the_gene='YES')

    'distance_to_TSS' and 'in_the_gene' are arbitrary stored data.

    The results can be written to a file using the 'output' method,
    specifying a file name and a header line which names the arbitrary
    stored data items to write in each column, e.g.

    >>> result.output('analysis.out',('distance_to_TSS','in_the_gene'))

    """
    def __init__(self):
        """Create a new AnalysisResult instance.

        Use the 'addResult' method to populate."""
        # List of results
        self.results = []
        # List of fields
        self.fields = []

    def addResult(self,peak_data,feature_data,**args):
        """Add a result to the AnalysisResult instance.

        A result associates a ChIP peak with a transcript/feature,

        Arguments:
          peak_data: a Peak object
          feature_data: a Feature object
          
        In addition this method will accept arbitrary named arguments
        which can be supplementary data associated with the result.

        """
        # Create a new dictionary
        result = {'peak': peak_data,
                  'feature': feature_data }
        # Add the arbitrary items
        for arg in args:
            if arg not in self.fields:
                self.fields.append(arg)
            result[arg] = args[arg]
        self.results.append(result)

    def output(self,filename,header,pad=''):
        """Write results out to a file

        The results are written to a tab-delimited file named 'filename',
        with the items in 'header' defining the content of each line.
        These header items should be a subset of the custom items
        stored by calls to addResult.

        Arguments:
          filename: name of the file to write to
          header: list of items from each result line to be written
          pad: (optional) string to use for padding empty fields

        """
        fp = open(filename,'w')
        # Write header line
        fp.write('#'+'\t'.join(header)+'\n')
        # Write results
        for result in self.results:
            items = []
            for field in header:
                try:
                    items.append(str(result[field]))
                except KeyError:
                    items.append(pad)
            fp.write('\t'.join(items)+'\n')
        fp.close()

    def output_xls(self,xls,title,header,pad=''):
        """Write results to an XLS spreadsheet
        """
        ws = xls.addSheet(title)
        ws.addText('#'+'\t'.join(header))
        # Keep track of line number and sheet number
        row_index = 1
        nsheet = 1
        # Loop over all results
        for result in self.results:
            items = []
            for field in header:
                try:
                    item = str(result[field])
                except KeyError:
                    item = pad
                # Split too-long items into multiple items/columns
                # This is because the Spreadsheet writer truncates cells
                # that exceed the spreadsheet cell character limit
                char_limit = Spreadsheet.MAX_LEN_WORKSHEET_CELL_VALUE
                while len(item) > char_limit:
                    logging.warning("Split value across multiple cells in sheet '%s' row %d"
                                    % (ws.title,row_index))
                    # Split on ';'
                    try:
                        # Locate nearest semicolon to the character limit
                        i = item[:char_limit].rindex(';')
                        items.append(item[:i])
                        item = item[i:].strip(';')
                    except ValueError:
                        # Unable to locate semicolon so split on the
                        # character limit
                        items.append(item[:char_limit])
                        item = item[char_limit:]
                items.append(item)
            ws.addText('\t'.join(items))
            # Update line number and check if the maximum number of lines has been exceeded
            row_index += 1
            if row_index == Spreadsheet.MAX_NUMBER_ROWS_PER_WORKSHEET:
                logging.warning("Maximum number of rows in XLS sheet '%s' exceeded (%d)" %
                                (ws.title,Spreadsheet.MAX_NUMBER_ROWS_PER_WORKSHEET))
                # Make a new spreadsheet for the excess rows
                nsheet += 1
                # Create a unique title
                tmp_title = truncate_text("%s(%d)" % (title,nsheet),
                                          Spreadsheet.MAX_LEN_WORKSHEET_TITLE)
                ws = xls.addSheet(tmp_title)
                logging.warning("Created new sheet '%s' to store additional results" %
                                ws.title)
                # Add the header
                ws.addText('#'+'\t'.join(header))
                # Reset the line index counter
                row_index = 1

    def __getitem__(self,key):
        return self.results[key]

    def __len__(self):
        return len(self.results)

#######################################################################
# Analysis functions
#######################################################################

def AnalyseClosestTranscriptsToPeaksInEachDirection(peaks,features):
    """Get closest up- and down-stream transcripts for + and - strands to ChIP peak"
    """
    print "Analysis #1a: ChIP-seq perspective"
    print "Get closest upstream and downstream transcripts to each ChIP peak"
    print "on both + and - strands"
    print ""
    results = []
    for peak in peaks:
        print "\tMatching to ChIP peak: %s" % peak
        # Initialise the hits
        print "\t%d total matches in RNA-seq data" % \
            len(features.filterByChr(peak.chrom))
        # Process the subset of features with matching chromosome
        # on the + strand
        closest_upstream_on_plus_strand = None
        closest_downstream_on_plus_strand = None
        print "\t%d matches on + strand" % \
            len(features.filterByChr(peak.chrom).filterByStrand('+'))
        for feature in features.filterByChr(peak.chrom).filterByStrand('+'):
            # Get the distance from the ChIP peak
            distance = feature.getTSS() - peak.start
            # Check if this is closer than anything seen previously
            # for each case
            if distance >= 0:
                # Upstream (i.e. +ve distance)
                closest_upstream_on_plus_strand = \
                    GetNearestTranscriptToPeak(closest_upstream_on_plus_strand,\
                                               feature,peak)
            else:
                # Downstream (i.e. -ve distance)
                closest_downstream_on_plus_strand = \
                    GetNearestTranscriptToPeak(closest_downstream_on_plus_strand,\
                                               feature,peak)
        # Process the subset of RNA-seq data with matching chromosome
        # on the - strand
        closest_upstream_on_minus_strand = None
        closest_downstream_on_minus_strand = None
        print "\t%d matches on - strand" % \
            len(features.filterByChr(peak.chrom).filterByStrand('-'))
        for feature in features.filterByChr(peak.chrom).filterByStrand('-'):
            # Get the distance from the ChIP peak
            distance = feature.getTSS() - peak.start
            # Check if this is closer than anything seen previously
            # for each case
            if distance >= 0:
                # Upstream (i.e. +ve distance)
                closest_upstream_on_minus_strand =  \
                    GetNearestTranscriptToPeak(closest_upstream_on_minus_strand,
                                               feature,peak)
            else:
                # Downstream (i.e. -ve distance)
                closest_downstream_on_minus_strand = \
                    GetNearestTranscriptToPeak(closest_downstream_on_minus_strand,
                                              feature,peak)
        # Store the data about the hits to output later
        nearest_positions = {
            "upstream +": closest_upstream_on_plus_strand,
            "downstream +": closest_downstream_on_plus_strand,
            "upstream -": closest_upstream_on_minus_strand,
            "downstream +": closest_downstream_on_minus_strand,
            }
        for hit in nearest_positions.keys():
            feature = nearest_positions[hit]
            if not feature:
                # No data stored for this position
                print "*** No hit for %s ***" % hit
                continue
            # Derive the distances
            distance_to_TSS = feature.getTSS() - peak.start
            distance_to_TES = feature.getTES() - peak.start
            # ChIP peak is in the gene?
            if feature.containsPosition(peak.start):
                in_the_gene = "YES"
            else:
                in_the_gene = "NO"
            # Add to results list
            results.append([peak.chrom,
                            peak.start,
                            feature.id,
                            distance_to_TSS,
                            distance_to_TES,
                            feature.strand,
                            in_the_gene])

    # Write the results to file
    results.insert(0,["#chr","start","geneID","distance_to_TSS","distance_to_TES",
                      "strand","in_the_gene"])
    output_results("RNA-seq_to_ChIP-seq.txt",results)

def AnalyseNearestTSSToSummits(peaks,features,max_distance,
                               max_closest=4,
                               filename=None,
                               xls=None,
                               pad_output=False):
    """Find nearest features to a set of ChIP peaks

    Given a set of ChIP peaks in a PeakSet object, and a set of
    gene transcripts in a FeatureSet object, this analysis aims to
    match genes to peaks.

    The procedure is:

    - For each ChIP peak:
    -- Consider only genes that are flagged as "significant" and
       with the same chromosome
    -- Report the nearest genes (by TSS) within a cutoff distance
       of the peak position (up to 'max_closest' hits)
    -- For each hit:
    --- Report all the genes (by TSS) from the full list of genes
        which lie between the ChIP peak and the TSS of the hit

    Arguments:
      chip_seq: ChIP peaks in PeakSet object
      features: RNA-seq gene transcripts in FeatureSet object
      max_distance: cutoff distance for reporting nearest genes
        to ChIP peaks (in units of bases)
      max_closest: (optional) maximum number of peaks to report
        (defaults to 4)
      filename: (optional) if not None then specifies the file
        name to write the results to as tab-delimited data
      xls: (optional) if not None then specifies the XLS file name
        to add the results to as a new sheet
      pad_output: (optional) if True then always report max_closest
        lines of output, and pad with blank lines if fewer genes
        were actually found

    """
    # Loop over ChIP peaks and sort features into order for each,
    # based on the absolute distance of their TSS from the peak
    # After sorting the top results will be the closest
    logging.debug("Starting AnalyseNearestTranscriptsToPeaks:")

    # Check that we have summit data for ChIP peaks
    if not peaks.isSummit():
        logging.warning("The supplied ChIP data only defines regions")
        logging.warning("This analysis is intended to work with summit data")

    # Create subset of "significant" i.e. flagged features
    if features.isFlagged():
        significant_features = features.filterByFlag(matchFlag=1)
        logging.debug("%d RNA-seq records after filtering on flag" % \
                          len(significant_features))
    else:
        significant_features = features

    # Run the analysis on the significant RNA-seq transcripts
    results = AnalysisResult()
    for peak in peaks:
        logging.debug("\t%s" % peak)
        features_chr = significant_features.\
            filterByChr(peak.chrom).\
            filterByTSS(peak.start-max_distance,\
                            peak.start+max_distance).\
                            sortByDistanceFrom(peak.start)
        features_chr_full = features.filterByChr(peak.chrom)
        closest = features_chr[:max_closest]
        i = 0
        for feature in closest:
            logging.debug("\t\t%s\t%d" % (feature,
                                          (feature.getTSS()-peak.start)))
            # Get the transcripts from the full list which
            # lie between the transcript TSS and the ChIP peak
            # These can be on either strand
            transcripts = features_chr_full.filterByTSS(peak.start,
                                                        feature.getTSS(),
                                                        exclude_limits=True)
            # Put them into distance order from the peak
            transcripts.sortByDistanceFrom(peak.start)
            # Make a list of just the ids
            transcript_ids = []
            for transcript in transcripts:
                transcript_ids.append(transcript.id)
                logging.debug("\t\t\t%s" % transcript.id)
            # Derive the distances
            distance_to_TSS = feature.getTSS() - peak.start
            distance_to_TES = feature.getTES() - peak.start
            # ChIP peak is in the gene?
            if feature.containsPosition(peak.start):
                in_the_gene = "YES"
            else:
                in_the_gene = "NO"
            # "Nearest" column reports "1 of 4", "2 of 4" etc
            i += 1
            nearest = "%d of %d" % (i,len(closest))
            # Add to the analysis result
            results.addResult(peak,feature,
                              chr=peak.chrom,
                              start=peak.start,
                              geneID=feature.id,
                              nearest=nearest,
                              TSS=feature.getTSS(),
                              distance_to_TSS=distance_to_TSS,
                              distance_to_TES=distance_to_TES,
                              strand=feature.strand,
                              in_the_gene=in_the_gene,
                              transcripts_inbetween=len(transcripts),
                              transcript_ids_inbetween=\
                                  ';'.join(transcript_ids))
        n_hits = len(closest)
        if n_hits == 0:
            # Report peaks with no significant genes in the cut-off region
            results.addResult(peak,None,
                              chr=peak.chrom,
                              start=peak.start)
            n_hits = 1
            logging.debug("\t\tNo transcripts found")
        if pad_output:
            # Pad with blank lines
            for i in range(n_hits,max_closest):
                results.addResult(peak,None)
        logging.debug("")
    # Write the results to file
    if filename:
        results.output(filename,('chr','start','geneID','nearest','TSS',
                                 'distance_to_TSS',
                                 'distance_to_TES',
                                 'strand','in_the_gene',
                                 'transcripts_inbetween',
                                 'transcript_ids_inbetween'),pad='.')
    # Write the results to a spreadsheet
    if xls: results.output_xls(xls,'TSSToSummits',
                               ('chr','start','geneID','nearest','TSS',
                                'distance_to_TSS',
                                'distance_to_TES',
                                'strand','in_the_gene',
                                'transcripts_inbetween',
                                'transcript_ids_inbetween'),pad='.')
    # Return results object
    return results

def AnalyseNearestPeaksToTranscripts(features,peaks,window_width,
                                     filename=None,xls=None):
    """Find nearest ChIP peaks to RNA-seq gene transcripts/features

    Given a set of ChIP peaks in a PeakSet object, and a set of
    gene transcripts in a FeatureSet object, this analysis aims to
    match peaks to genes.

    The procedure is:

    - For each RNA-seq gene transcript:
    -- If the gene is flagged as "significant" then:
    --- Report the ChIP peaks on the same chromosome which lie
        within +/- window_width of the gene TSS

    Note that all genes are reported in the output file.

    Arguments:
      features: RNA-seq gene transcripts in FeatureSet object
      chip_seq: ChIP peaks in PeakSet object
      window_width: cutoff distance for reporting nearest ChIP peaks
        to RNA-seq gene TSS positions (in units of bases)
      filename: (optional) if not None then specifies the file
        name to write the results to as tab-delimited data
      xls: (optional) if not None then specifies the XLS file name
        to add the results to as a new sheet
    """
    logging.debug("Starting AnalyseNearestPeaksToTranscripts:")
    # Check that we have summit data for ChIP peaks
    if not peaks.isSummit():
        logging.warning("The supplied ChIP data only defines regions")
        logging.warning("This analysis is intended to work with summit data")
    # Get flag status of feature data
    features_is_flagged = features.isFlagged()
    # Do the analysis
    results = AnalysisResult()
    max_peaks = 0
    for feature in features:
        logging.debug("\t%s" % feature)
        # Define a "window" either side of the TSS
        window = (feature.getTSS()-window_width,
                  feature.getTSS()+window_width)
        logging.debug("\t\tWindow: %d - %d" % (window[0],window[1]))
        # Get ChIP peaks on the matching chromosome that lie in this window
        # Only do this for the significant (i.e. flagged) genes, unless
        # the data is not flagged
        if not features_is_flagged or feature.flag == 1:
            filtered_peaks = peaks.\
                filterByChr(feature.chrom).filterByPosition(window[0],window[1])
            filtered_peaks.sortByDistanceFrom(feature.getTSS())
        else:
            # No hits - empty result object
            filtered_peaks = PeakSet()
            logging.debug("\t\tNo peaks found")
        # Store results
        result = { 'geneID': feature.id,
                   'chr_RNA': feature.chrom,
                   'start': feature.start,
                   'end': feature.end,
                   'strand': feature.strand,
                   'differentially_expressed': feature.flag,
                   'number_of_peaks': len(filtered_peaks) }
        for i in range(len(filtered_peaks)):
            peak = filtered_peaks[i]
            # Derive the distance
            distance = feature.getTSS() - peak.start
            # Add to results for output
            result['chr_ChIP_'+str(i+1)] = peak.chrom
            result['summit_'+str(i+1)] = peak.start
            result['distance_'+str(i+1)] = distance
            logging.debug("\t\t%s (%d)" % (peak,distance))
        logging.debug("")
        # Add line to results to be written to the output file
        results.addResult(filtered_peaks,feature,**result)
        # Keep track of the maximum number of peaks for output
        max_peaks = max(len(filtered_peaks),max_peaks)
    # Construct header line for output
    if filename or xls:
        header = ["geneID","chr_RNA","start","end","strand",
                  "differentially_expressed","number_of_peaks"]
        for i in range(max_peaks):
            header.extend(["chr_ChIP_"+str(i+1),
                           "summit_"+str(i+1),
                           "distance_"+str(i+1)])
    # Write to output file
    if filename: results.output(filename,header,pad='---')
    # Write the results to a spreadsheet
    if xls: results.output_xls(xls,'PeaksToTranscripts',
                               header,pad='---')
    # Finished
    return results

def AnalyseNearestTranscriptsToPeakEdges(peaks,features,
                                         promoter_region=(10000,2500),
                                         max_closest=4,
                                         max_distance=0,
                                         TSS_only=False,
                                         filename=None,
                                         xls=None,
                                         pad_output=False):
    """Find nearest RNA transcripts to a set of ChIP peaks based on edges

    Given a set of ChIP peaks in a PeakSet object, and a set of
    gene transcripts in a FeatureSet object, this analysis aims to
    match genes to peaks.

    The procedure is:

    - For each ChIP peak:
    -- Consider only genes on the same chromosome
    -- Report the nearest genes by finding either the nearest TSS or TES, or
       just the nearest TSS, to the peak edges (start/stop) (up to four hits)
       (Both TSS and TES are used if the 'TSS_only' flag is set to False;
       otherwise only TSS is used).
    -- For each hit:
    --- Report whether the peak overlaps the transcript region
    --- Report whether the peak overlaps the promoter region

    Arguments:
      peaks: ChIP peaks in PeakSet object
      features: RNA-seq gene transcripts in FeatureSet object
      promoter_region: (optional) a tuple (leading,trailing) which
        specifies the promoter region by the "leading" (upstream) and
        "trailing" (downstream) edges of the region with respect to the
        gene TSS. Defaults to (10000,2500).
      max_closest: (optional) maximum number of peaks to report
        (defaults to 4).
      max_distance: (optional) set cutoff gene-to-edge distance, beyond
        which genes are not reported. Defaults to zero (don't apply a
        cutoff and report all genes regardless of distance).
      TSS_only: if set to True then determine closest distance from TSS
        only, otherwise use both TSS and TES (the default).
      filename: (optional) if specified then write output to file with
        with this name; also write a 'summary' file with just the top
        hit for each ChIP peak.
      xls: (optional) if not None then specifies the XLS file name
        to add the results to as a new sheet
      pad_output: (optional) if True then always report max_closest
        lines of output, and pad with blank lines if fewer genes
        were actually found
    """
    logging.debug("Starting AnalyseNearestTranscriptsToPeakEdges:")
    # Check that we have region data for ChIP peaks
    if peaks.isSummit():
        logging.warning("The supplied ChIP data only defines summits")
        logging.warning("This analysis is intended to work with regions")
    # Create subset of "significant" i.e. flagged RNA-seq data
    if features.isFlagged():
        significant_features = features.filterByFlag(matchFlag=1)
        logging.debug("%d RNA-seq records after filtering on flag" % \
                          len(significant_features))
    else:
        significant_features = features
    # Check mode of operation
    if not TSS_only:
        logging.debug("Use both TSS and TES in analysis (TSS_only = %s)" %
                      TSS_only)
    else:
        logging.debug("Use only TSS in analysis (TSS_only = %s)" % TSS_only)
    results = AnalysisResult()
    # Loop over all peaks
    for peak in peaks:
        logging.debug("\t%s" % peak)
        # Sort transcripts into order of the nearest TSS or TES to
        # the edge of the peak/binding region
        # Note that this method does not discriminate between situations
        # where the transcript lies partially or wholly within the
        # binding region
        features_chr = significant_features.filterByChr(peak.chrom)
        if not TSS_only:
            features_chr = features_chr.sortByClosestEdgeTo(peak.start,
                                                            peak.end)
        else:
            features_chr = features_chr.sortByClosestTSSTo(peak.start,
                                                           peak.end)
        # Get the closest peaks (filter by distance and max number)
        closest = []
        for feature in features_chr:
            # Get closest edge distances
            distance = feature.getClosestEdgeDistanceTo(peak.start,
                                                        peak.end)
            # Apply distance cutoff
            if max_distance > 0 and distance > max_distance:
                # Break out of the loop
                logging.debug("Exceeded maximum distance, stopping")
                break
            distance_to_TSS = feature.getClosestTSSDistanceTo(\
                peak.start,peak.end)
            distance_to_TES = feature.getClosestTESDistanceTo(\
                peak.start,peak.end)
            # Determine if transcript and peak overlap
            binding_region = (peak.start,peak.end)
            if regions_overlap(binding_region,(feature.getTSS(),
                                               feature.getTES())):
                overlap_transcript = 1
            else:
                overlap_transcript = 0
            # Determine if promoter region and peak overlap?
            promoter = feature.getPromoterRegion(*promoter_region)
            if regions_overlap(binding_region,promoter):
                overlap_promoter = 1
            else:
                overlap_promoter = 0
            # Add to results
            results.addResult(peak,feature,
                              chr=peak.chrom,
                              start=peak.start,
                              end=peak.end,
                              geneID=feature.id,
                              TSS=feature.getTSS(),
                              TES=feature.getTES(),
                              strand=feature.strand,
                              dist_closest_edge=distance,
                              dist_TSS=distance_to_TSS,
                              dist_TES=distance_to_TES,
                              overlap_transcript=overlap_transcript,
                              overlap_promoter=overlap_promoter)
            # Report
            logging.debug("\t\t%s (%d,%s,%s)" % \
                             (feature,
                              distance,
                              overlap_transcript,
                              overlap_promoter))
            # Check whether we've reached the limit to report
            closest.append(feature)
            if len(closest) == max_closest:
                logging.debug("Found %d closest" % max_closest)
                break
        # Finished loop
        n_hits = len(closest)
        if n_hits == 0:
            logging.debug("\t\tNo transcripts found")
        # Pad with blank lines, if requested
        if pad_output:
            for i in range(n_hits,max_closest):
                results.addResult(peak,None,
                                  chr=peak.chrom,
                                  start=peak.start,
                                  end=peak.end)
        logging.debug("")
    # Construct header line and summary results for output
    if filename or xls:
        # Header line
        header = ["chr","start","end","geneID","strand","TSS","TES",
                  "dist_closest_edge",
                  "dist_TSS","dist_TES",
                  "overlap_transcript",
                  "overlap_promoter"]
        # Summary results: one gene per ChIP peak (i.e. top hit for
        # each peak)
        summary_results = AnalysisResult()
        last_peak = None
        for result in results:
            if result['peak'] != last_peak and result['feature'] is not None:
                summary_results.addResult(result['peak'],
                                          result['feature'],
                                          **result)
                last_peak = result['peak']
            else:
                pass
    # Write to output file
    if filename:
        # Full results
        results.output(filename,header,pad='.')
        # Summary
        summary_filename = os.path.splitext(filename)[0]+"_summary"+\
            os.path.splitext(filename)[1]
        summary_results.output(summary_filename,header)
    # Write the results to a spreadsheet
    if xls:
        # Full results
        if not TSS_only:
            xls_title = 'TranscriptsToPeakEdges'
        else:
            xls_title = 'TSSToPeakEdges'
        results.output_xls(xls,xls_title,header,pad='.')
        # Summary results
        xls_title += '(summary)'
        summary_results.output_xls(xls,xls_title,header)
            
    # Return (full) analysis results
    return results

#######################################################################
# Helper Functions
#######################################################################

def output_results(filen,results):
    """Write list of results to tab-delimited file

    Given a (Python) list of results, each item is written to the
    specified file as a set of tab-delimited fields.

    If the file already exists then it will be overwritten, and
    there are no checks on the integrity of the data.
    """
    fp = open(filen,'w')
    for result in results:
        items = [str(x) for x in result]
        fp.write('\t'.join(items)+'\n')
    fp.close()

