#######################################################################
# Descriptions of each method for XLS notes
#######################################################################

preamble = \
"""<style font=bold bgcolor=gray25>%s: integrated analyses of RNA-Seq and ChIP-Seq data</style>

The following analyses have been performed and are reported in this spreadsheet.
"""

credits = \
"""<style font=bold bgcolor=gray25>Credits</style>
Produced by %s on %s
Bioinformatics Core Facility, Faculty of Life Sciences, University of Manchester 
http://fls-bioinformatics-core.github.com/RnaChipIntegrator/
"""

nearest_TSS_to_summits = \
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

nearest_transcripts_to_peak_edges = \
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

nearest_tss_to_peak_edges = \
"""<style font=bold bgcolor=gray25>TSS to Peak Edges</style>
Find the nearest transcripts (up to 4) with the smallest distance from their TSS to the nearest peak edge.
The input parameters and output fields are the same as for the "Transcripts to Peak Edges" analysis above.

<style font=bold bgcolor=gray25>TSS to Peak Edges (summary)</style>
Same as "TSS To Peak Edges" above but lists only the single nearest transcript to each peak.
"""

nearest_peaks_to_transcripts = \
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
