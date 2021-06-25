########################################################################
#
# Inlined test data for unit tests
#
#########################################################################
#
import io
#
# Transcripts-ex1.txt
transcripts_ex1 = \
"""CG9130-RB	chr3L	1252012	1255989	+	1
CG8616-RA	chr3L	7231642	7233354	+	0
CG32847-RB	chr3L	15114722	15115217	+	0
CG14448-RA	chr3L	22781539	22782354	+	1
CG32065-RA	chr3L	10428770	10441014	+	1
CG10541-RA	chr3L	5787500	5789502	-	0
CG10583-RA	chr3L	5652794	5655342	-	0
CG34099-RB	chr3L	182522	184423	-	0
CG40225-RA	chr3LHet	2276035	2276749	-	1
CG13051-RA	chr3L	16257914	16258166	-	0"""
#
# Transcripts-ex2.txt ("bad" RNA transcripts file)
transcripts_ex2 = \
"""CG9130-RB	chr3L	1252012	1.2E+06	+	1
CG8616-RA	chr3L	7231642	7.2E+06	+	0
CG32847-RB	chr3L	15114722	1.5E+07	+	0
CG14448-RA	chr3L	22781539	2.2E+07	+	1
CG32065-RA	chr3L	10428770	1.0E+07	+	1
CG10541-RA	chr3L	5787500	5.7E+06	-	0
CG10583-RA	chr3L	5652794	5.6E+06	-	0
CG34099-RB	chr3L	182522	1.8E+05	-	0
CG40225-RA	chr3LHet	2276035	2.2E+06	-	1
CG13051-RA	chr3L	16257914	1.6E+0	-	0"""
#
# Transcripts-ex2a.txt ("bad" RNA transcripts file)
transcripts_ex2a = \
"""CG9130-RB	chr3L	1252012	1200000	+	1
CG8616-RA	chr3L	7231642	7200000	+	0
CG32847-RB	chr3L	15114722	15000000	+	0
CG14448-RA	chr3L	22781539	22000000	+	1
CG32065-RA	chr3L	10428770	10000000	+	1
CG10541-RA	chr3L	5787500	5700000	-	0
CG10583-RA	chr3L	5652794	5600000	-	0
CG34099-RB	chr3L	182522	180000	-	0
CG40225-RA	chr3LHet	2276035	2200000	-	1
CG13051-RA	chr3L	16257914	16000000	-	0"""
#
# Transcripts-ex3.txt 
transcripts_ex3 = \
"""Probe_name	chromosome	start	stop	Strand	DESeq_1200
CG31973-RD	chr2L	25402	59243	-	1
CG31973-RA	chr2L	25402	59243	-	0
CG31973-RE	chr2L	25402	59243	-	0
CG31973-RB	chr2L	25402	59243	-	0
CG2674-RC	chr2L	107926	114433	+	1
CG2674-RG	chr2L	109608	114433	+	1
CG2674-RJ	chr2L	108094	114434	+	1
CG2674-RD	chr2L	107926	114433	+	0
CG2674-RE	chr2L	106903	114433	+	0
CG2674-RH	chr2L	107926	114433	+	0
CG2674-RB	chr2L	107760	114433	+	0
CG2674-RA	chr2L	107760	114433	+	0
CG2674-RF	chr2L	108089	114434	+	0
CG2674-RI	chr2L	108089	114434	+	0
CG3625-RB	chr2L	283385	285777	-	1
CG3625-RC	chr2L	283385	291011	-	1
CG3625-RA	chr2L	283385	286528	-	0
CG2851-RA	chr2L	583540	594688	-	1
CG2851-RB	chr2L	583540	594811	-	0
CG2762-RA	chr2L	523467	540542	+	1
CG17941-RA	chr2L	640021	714969	-	1"""
#
# Transcripts-ex4.txt 
transcripts_ex4 = \
"""Probe_name	chromosome	start	stop	Strand	DESeq_1200
CG31973-RD	chr2L	25402	59243	-	1
CG31973-RA	chr2L	25402	59243	-	0
CG14026-RA	chr2L	5218996	5271354	-	1
CG18024-RA	chr2L	8825625	8829671	+	1"""
#
# ChIP_peaks-ex1.txt 
chip_peaks_ex1 = \
"""chr3L	4252919	4252920
chr3L	9502640	9502641
chr3L	12139192	12139193
chr3L	14983597	14983598
chr3L	17004143	17004144"""
#
# ChIP_peaks_multi_columns-ex1.txt
chip_peaks_multi_columns_ex1 = \
"""peak1	chr3L	66.1	4252919	4252920
peak2	chr3L	49.9	9502640	9502641
peak3	chr3L	50.0	12139192	12139193
peak4	chr3L	89.2	14983597	14983598
peak5	chr3L	63.4	17004143	17004144"""
#
# ChIP_peaks-ex2.txt 
chip_peaks_ex2 = \
"""chr2L	66811	66812
chr2L	91509	91510
chr2L	249177	249178
chr2L	481213	481214
chr2L	518665	518666
chr3L	47666	47667
chr3L	162061	162062
chr3L	258189	258190
chr3L	264677	264678
chr3L	319804	319805"""
#
# ChIP_peaks-ex3.txt
chip_peaks_ex3 = \
"""#chr	SummitMidpoint 	SummitMidpoint+1
chr2L	66811	66812
chr2L	249177	249178
chr2L	605950	605951"""
#
# ChIP_peaks-ex3.txt
chip_peaks_binding_region_ex3 = \
"""#chr	Summit-100	Summit+100
chr2L	66711	66911
chr2L	249077	249277
chr2L	605850	606050"""
#
# ChIP_peaks-ex4.txt 
chip_peaks_ex4 = \
"""chr2L	66811	66812
chr2L	8825470	8825471
chr2L	8832523	8832524
chr2L	8842036	8842037"""
#
# ChIP_peaks-ex5.txt
chip_peaks_ex5 = \
""""chr"	"start"	"end"
"chr8"	6179401	6180699	
"chr3"	86512837	86514328
"chr4"	14116403	14118581
"chr4"	107209680	107211432
"chr12"	28921259	28923030"""
#
# ChIP_peaks-ex6.txt (same start/end position)
chip_peaks_ex6 = \
"""chr2L	66811	66812
chr2L	8825470	8825470
chr2L	8832523	8832524
chr2L	8842036	8842037"""
#
# ChIP_peaks-ex7.txt (end position before start)
chip_peaks_ex7 = \
"""chr2L	66811	66812
chr2L	8825470	8825471
chr2L	8832524	8832523
chr2L	8842036	8842037"""
#
# NearestTranscriptToPeak-ex3.txt
nearest_transcript_to_peak_ex3 = \
"""#chr	start	geneID	nearest	TSS	distance_to_TSS	distance_to_TES	strand	in_the_gene	transcripts_inbetween	transcript_ids_inbetween
chr2L	66811	CG31973-RD	1 of 4	59243	-7568	-41409	-	NO	0	
chr2L	66811	CG2674-RC	2 of 4	107926	41115	47622	+	NO	3	CG2674-RE;CG2674-RB;CG2674-RA
chr2L	66811	CG2674-RJ	3 of 4	108094	41283	47623	+	NO	8	CG2674-RE;CG2674-RB;CG2674-RA;CG2674-RC;CG2674-RD;CG2674-RH;CG2674-RF;CG2674-RI
chr2L	66811	CG2674-RG	4 of 4	109608	42797	47622	+	NO	9	CG2674-RE;CG2674-RB;CG2674-RA;CG2674-RC;CG2674-RD;CG2674-RH;CG2674-RF;CG2674-RI;CG2674-RJ
chr2L	249177	CG3625-RB	1 of 2	285777	36600	34208	-	NO	0	
chr2L	249177	CG3625-RC	2 of 2	291011	41834	34208	-	NO	2	CG3625-RB;CG3625-RA
chr2L	605950	CG2851-RA	1 of 3	594688	-11262	-22410	-	NO	1	CG2851-RB
chr2L	605950	CG2762-RA	2 of 3	523467	-82483	-65408	+	NO	2	CG2851-RB;CG2851-RA
chr2L	605950	CG17941-RA	3 of 3	714969	109019	34071	-	NO	0	"""
#
# NearestTranscriptToPeakEdge-ex3.txt
nearest_transcript_to_peak_edge_ex3 = \
"""#chr	start	end	geneID	strand	TSS	TES	dist_closest_edge	dist_TSS	dist_TES	overlap_transcript	overlap_promoter
chr2L	66711	66911	CG31973-RD	-	59243	25402	7468	7468	41309	0	1
chr2L	66711	66911	CG31973-RA	-	59243	25402	7468	7468	41309	0	1
chr2L	66711	66911	CG31973-RE	-	59243	25402	7468	7468	41309	0	1
chr2L	66711	66911	CG31973-RB	-	59243	25402	7468	7468	41309	0	1
chr2L	249077	249277	CG3625-RB	-	285777	283385	34108	36500	34108	0	0
chr2L	249077	249277	CG3625-RC	-	291011	283385	34108	41734	34108	0	0
chr2L	249077	249277	CG3625-RA	-	286528	283385	34108	37251	34108	0	0
chr2L	249077	249277	CG2674-RJ	+	108094	114434	134643	140983	134643	0	0
chr2L	605850	606050	CG2851-RB	-	594811	583540	11039	11039	22310	0	0
chr2L	605850	606050	CG2851-RA	-	594688	583540	11162	11162	22310	0	0
chr2L	605850	606050	CG17941-RA	-	714969	640021	33971	108919	33971	0	0
chr2L	605850	606050	CG2762-RA	+	523467	540542	65308	82383	65308	0	0"""
#
# NearestTranscriptTSSToPeakEdge-ex3.txt
nearest_transcript_tss_to_peak_edge_ex3 = \
"""#chr	start	end	geneID	strand	TSS	TES	dist_closest_edge	dist_TSS	dist_TES	overlap_transcript	overlap_promoter
chr2L	66711	66911	CG31973-RD	-	59243	25402	7468	7468	41309	0	1
chr2L	66711	66911	CG31973-RA	-	59243	25402	7468	7468	41309	0	1
chr2L	66711	66911	CG31973-RE	-	59243	25402	7468	7468	41309	0	1
chr2L	66711	66911	CG31973-RB	-	59243	25402	7468	7468	41309	0	1
chr2L	249077	249277	CG3625-RB	-	285777	283385	34108	36500	34108	0	0
chr2L	249077	249277	CG3625-RA	-	286528	283385	34108	37251	34108	0	0
chr2L	249077	249277	CG3625-RC	-	291011	283385	34108	41734	34108	0	0
chr2L	249077	249277	CG2674-RG	+	109608	114433	134644	139469	134644	0	0
chr2L	605850	606050	CG2851-RB	-	594811	583540	11039	11039	22310	0	0
chr2L	605850	606050	CG2851-RA	-	594688	583540	11162	11162	22310	0	0
chr2L	605850	606050	CG2762-RA	+	523467	540542	65308	82383	65308	0	0
chr2L	605850	606050	CG17941-RA	-	714969	640021	33971	108919	33971	0	0"""
#
# NearestTranscriptToPeakEdge_diff_expression-ex3.txt
nearest_transcript_to_peak_edge_diff_expression_ex3 = \
"""#chr	start	end	geneID	strand	TSS	TES	dist_closest_edge	dist_TSS	dist_TES	overlap_transcript	overlap_promoter
chr2L	66711	66911	CG31973-RD	-	59243	25402	7468	7468	41309	0	1
chr2L	66711	66911	CG2674-RC	+	107926	114433	41015	41015	47522	0	0
chr2L	66711	66911	CG2674-RJ	+	108094	114434	41183	41183	47523	0	0
chr2L	66711	66911	CG2674-RG	+	109608	114433	42697	42697	47522	0	0
chr2L	249077	249277	CG3625-RB	-	285777	283385	34108	36500	34108	0	0
chr2L	249077	249277	CG3625-RC	-	291011	283385	34108	41734	34108	0	0
chr2L	249077	249277	CG2674-RJ	+	108094	114434	134643	140983	134643	0	0
chr2L	249077	249277	CG2674-RC	+	107926	114433	134644	141151	134644	0	0
chr2L	605850	606050	CG2851-RA	-	594688	583540	11162	11162	22310	0	0
chr2L	605850	606050	CG17941-RA	-	714969	640021	33971	108919	33971	0	0
chr2L	605850	606050	CG2762-RA	+	523467	540542	65308	82383	65308	0	0
chr2L	605850	606050	CG3625-RC	-	291011	283385	314839	314839	322465	0	0"""
#
# NearestTranscriptTSSToPeakEdge_diff_expression-ex3.txt
nearest_transcript_tss_to_peak_edge_diff_expression_ex3 = \
"""#chr	start	end	geneID	strand	TSS	TES	dist_closest_edge	dist_TSS	dist_TES	overlap_transcript	overlap_promoter
chr2L	66711	66911	CG31973-RD	-	59243	25402	7468	7468	41309	0	1
chr2L	66711	66911	CG2674-RC	+	107926	114433	41015	41015	47522	0	0
chr2L	66711	66911	CG2674-RJ	+	108094	114434	41183	41183	47523	0	0
chr2L	66711	66911	CG2674-RG	+	109608	114433	42697	42697	47522	0	0
chr2L	249077	249277	CG3625-RB	-	285777	283385	34108	36500	34108	0	0
chr2L	249077	249277	CG3625-RC	-	291011	283385	34108	41734	34108	0	0
chr2L	249077	249277	CG2674-RG	+	109608	114433	134644	139469	134644	0	0
chr2L	249077	249277	CG2674-RJ	+	108094	114434	134643	140983	134643	0	0
chr2L	605850	606050	CG2851-RA	-	594688	583540	11162	11162	22310	0	0
chr2L	605850	606050	CG2762-RA	+	523467	540542	65308	82383	65308	0	0
chr2L	605850	606050	CG17941-RA	-	714969	640021	33971	108919	33971	0	0
chr2L	605850	606050	CG3625-RC	-	291011	283385	314839	314839	322465	0	0"""
#
# NearestPeakToTranscript-ex4.txt
nearest_peak_to_transcript_ex4 = \
"""#geneID	chr_RNA	start	end	strand	differentially_expressed	number_of_peaks	chr_ChIP	summit	distance	chr_ChIP	summit	distance	chr_ChIP	summit	distance
CG31973-RD	chr2L	25402	59243	-	1	1	chr2L	66811	-7568	---	---	---	---	---	---
CG31973-RA	chr2L	25402	59243	-	0	0	---	---	---	---	---	---	---	---	---
CG14026-RA	chr2L	5218996	5271354	-	1	0	---	---	---	---	---	---	---	---	---
CG18024-RA	chr2L	8825625	8829671	+	1	3	chr2L	8825470	155	chr2L	8832523	-6898	chr2L	8842036	-16411"""

########################################################################
#
# Utility functions for tests
#
#########################################################################

import os

def create_test_file(name,data):
    # Writes data to a file to be used in testing
    fp = io.open(name,'wt')
    for line in data.split('\n'):
        fp.write(u"%s\n" % line)
    fp.close()

def delete_test_file(name):
    # Deletes named test file
    if os.path.exists(name):
        os.remove(name)
