#     test_RnaChipIntegrator.py: unit tests for RnaChipIntegrator.py
#     Copyright (C) University of Manchester 2011 Peter Briggs
#
########################################################################
#
# test_RnaChipIntegrator.py
#
#########################################################################

import unittest, os
from RnaChipIntegrator import *

########################################################################
#
# Inlined test data
#
#########################################################################
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
# NearestTranscriptToPeak-ex3.txt
nearest_transcript_to_peak_ex3 = \
"""#chr	start	ID	nearest	TSS	distance_to_TSS	distance_to_TES	strand	in_the_gene	transcripts_inbetween	transcript_ids_inbetween
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
"""#chr	start	end	ID	strand	TSS	TES	dist_closest_edge	dist_TSS	dist_TES	overlap_transcript	overlap_promoter
chr2L	66711	66911	CG31973-RD	-	59243	25402	7468	7468	41309	0	0
chr2L	66711	66911	CG31973-RA	-	59243	25402	7468	7468	41309	0	0
chr2L	66711	66911	CG31973-RE	-	59243	25402	7468	7468	41309	0	0
chr2L	66711	66911	CG31973-RB	-	59243	25402	7468	7468	41309	0	0
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
"""#chr	start	end	ID	strand	TSS	TES	dist_closest_edge	dist_TSS	dist_TES	overlap_transcript	overlap_promoter
chr2L	66711	66911	CG31973-RD	-	59243	25402	7468	7468	41309	0	0
chr2L	66711	66911	CG31973-RA	-	59243	25402	7468	7468	41309	0	0
chr2L	66711	66911	CG31973-RE	-	59243	25402	7468	7468	41309	0	0
chr2L	66711	66911	CG31973-RB	-	59243	25402	7468	7468	41309	0	0
chr2L	249077	249277	CG3625-RB	-	285777	283385	34108	36500	34108	0	0
chr2L	249077	249277	CG3625-RA	-	286528	283385	34108	37251	34108	0	0
chr2L	249077	249277	CG3625-RC	-	291011	283385	34108	41734	34108	0	0
chr2L	249077	249277	CG2674-RG	+	109608	114433	134644	139469	134644	0	0
chr2L	605850	606050	CG2851-RB	-	594811	583540	11039	11039	22310	0	0
chr2L	605850	606050	CG2851-RA	-	594688	583540	11162	11162	22310	0	0
chr2L	605850	606050	CG2762-RA	+	523467	540542	65308	82383	65308	0	0
chr2L	605850	606050	CG17941-RA	-	714969	640021	33971	108919	33971	0	0"""
#
# NearestPeakToTranscript-ex4.txt
nearest_peak_to_transcript_ex4 = \
"""#ID	chr_RNA	start	end	strand	differentially_expressed	number_of_peaks	chr_ChIP	summit	distance	chr_ChIP	summit	distance	chr_ChIP	summit	distance
CG31973-RD	chr2L	25402	59243	-	1	1	chr2L	66811	-7568	---	---	---	---	---	---
CG31973-RA	chr2L	25402	59243	-	0	0	---	---	---	---	---	---	---	---	---
CG14026-RA	chr2L	5218996	5271354	-	1	0	---	---	---	---	---	---	---	---	---
CG18024-RA	chr2L	8825625	8829671	+	1	3	chr2L	8825470	155	chr2L	8832523	-6898	chr2L	8842036	-16411"""

########################################################################
#
# TestRNASeqDataLine
#
#########################################################################
class TestRNASeqDataLine(unittest.TestCase):

    def setUp(self):
        # Setup RNASeqDataLine objects to be used in the tests
        # Forward strand example
        self.rna_data = \
            RNASeqDataLine('CG9130-RB','chr3L','1252012','1255989','+')
        # Reverse strand example
        self.rna_data_2 = \
            RNASeqDataLine('CG13051-RA','chr3L','16257914','16258166','-')

    def test_contains_position(self):
        position = 1253000
        self.assertTrue(self.rna_data.containsPosition(position),
                        "Transcript should contain position")
        position = 4250000
        self.assertFalse(self.rna_data.containsPosition(position),
                         "Transcript should not contain position")
        position = 10000
        self.assertFalse(self.rna_data.containsPosition(position),
                         "Transcript should not contain position")

    def test_get_closest_edge_distance(self):
        # Single position
        position = 1200000
        # Single reference position
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position),
                         abs(self.rna_data.start-position),
                         "Incorrect closest edge distance #1")
        position = 1300000
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position),
                         abs(self.rna_data.end-position),
                         "Incorrect closest edge distance #2")

    def test_get_closest_edge_distance_outside_region(self):
        # Reference region (2 positions)
        position1 = 1100000
        position2 = 1200000
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2),
                         abs(self.rna_data.start-position2),
                         "Incorrect closest edge distance (region #1)")
        position1 = 1300000
        position2 = 1400000
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2),
                         abs(self.rna_data.end-position1),
                         "Incorrect closest edge distance (region #2)")

    def test_get_closest_edge_distance_partially_inside_region(self):
        # Partially inside reference region
        position1 = 1200000
        position2 = 1255000
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2),
                         abs(self.rna_data.end-position2),
                         "Incorrect closest edge distance (inside region #1)")
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2,
                                                                zero_inside_region=True),
                         0,
                         "Incorrect closest edge distance (inside region #2)")

    def test_get_closest_edge_distance_completely_inside_region(self):
        # Completely inside reference region
        position1 = 1250000
        position2 = 1300000
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2),
                         abs(self.rna_data.start-position1),
                         "Incorrect closest edge distance (inside region #3)")
        self.assertEqual(self.rna_data.getClosestEdgeDistanceTo(position1,
                                                                position2,
                                                                zero_inside_region=True),
                         0,
                         "Incorrect closest edge distance (inside region #4)")

    def test_get_promoter_region(self):
        leading = 10000
        trailing = 2500
        promoter = self.rna_data.getPromoterRegion(leading,trailing)
        self.assertEqual(promoter,
                         (self.rna_data.getTSS()-leading,
                          self.rna_data.getTSS()+trailing),
                         "Incorrect promoter region for + strand example")
        promoter = self.rna_data_2.getPromoterRegion(10000,2500)
        self.assertEqual(promoter,
                         (self.rna_data_2.getTSS()+leading,
                          self.rna_data_2.getTSS()-trailing),
                         "Incorrect promoter region for - strand example")

########################################################################
#
# TestRNASeqDataLine
#
#########################################################################
class TestRNASeqData(unittest.TestCase):

    def setUp(self):
        # Create input files for tests
        create_test_file('Transcripts-ex1.txt',transcripts_ex1)
        create_test_file('Transcripts-ex2.txt',transcripts_ex2)

    def tearDown(self):
        # Remove input files
        delete_test_file('Transcripts-ex1.txt')
        delete_test_file('Transcripts-ex2.txt')

    def test_reading_in_RNAseq_data(self):
        rna_seq = RNASeqData('Transcripts-ex1.txt')
        self.assertEqual(len(rna_seq),10,
                         "Wrong number of lines from RNA-seq file")
        self.assertTrue(rna_seq.isFlagged(),
                        "Data should be flagged")

    def test_filter_on_chromosome(self):
        rna_seq = RNASeqData('Transcripts-ex1.txt')
        rna_chr = rna_seq.filterByChr('chr3LHet')
        self.assertEqual(len(rna_chr),1,
                         "Wrong number of chromosomes")
        for rna_data in rna_chr:
            self.assertEqual(rna_data.chr,'chr3LHet',
                             "Wrong chromosome filtered")

    def test_filter_on_strand(self):
        rna_seq = RNASeqData('Transcripts-ex1.txt')
        rna_plus = rna_seq.filterByStrand('+')
        self.assertEqual(len(rna_plus),5,
                         "Wrong number of + strands")
        rna_minus = rna_seq.filterByStrand('-')
        self.assertEqual(len(rna_minus),5,
                         "Wrong number of - strands")

    def test_filter_on_flag(self):
        rna_seq = RNASeqData('Transcripts-ex1.txt')
        rna_flagged = rna_seq.filterByFlag(1)
        self.assertEqual(len(rna_flagged),4,
                         "Wrong number of flagged data lines")

    def test_getTSS(self):
        rna_seq = RNASeqData('Transcripts-ex1.txt')
        rna_plus = rna_seq.filterByStrand('+')
        for rna_data in rna_plus:
            self.assertTrue((rna_data.strand == '+' and
                             rna_data.start == rna_data.getTSS()),
                            "Incorrect TSS on + strand")
        rna_minus = rna_seq.filterByStrand('-')
        for rna_data in rna_minus:
            self.assertTrue((rna_data.strand == '-' and
                             rna_data.end == rna_data.getTSS()),
                            "Incorrect TSS on - strand")

    def test_filter_on_TSS(self):
        rna_seq = RNASeqData('Transcripts-ex1.txt')
        lower,upper = 5000000,10000000
        rna_tss = rna_seq.filterByTSS(upper,lower)
        self.assertEqual(len(rna_tss),3,
                         "Wrong number of transcripts filtered on TSS")
        for rna_data in rna_tss:
            self.assertTrue((rna_data.getTSS() >= lower and
                             rna_data.getTSS() <= upper),
                            "Transcript outside range")

    def test_sort_by_distance(self):
        rna_sort = RNASeqData('Transcripts-ex1.txt')
        position = 4250000
        # Do sort on distance
        # Sort is done in place, so assignment is not required
        # however the sort function should return a reference to
        # the initial object
        result = rna_sort.sortByDistanceFrom(position)
        self.assertEqual(result,rna_sort,
                         "Returned object doesn't match subject")
        # Check that each distance is greater than the previous one
        last_rna_data = None
        for rna_data in rna_sort:
            if not last_rna_data:
                last_rna_data = rna_data
        else:
            self.assertTrue((abs(rna_data.getTSS() - position) >=
                             abs(last_rna_data.getTSS() - position)),
                             "Sort by distance failed")

    def test_sort_by_closest_distance_to_edge(self):
        rna_sort = RNASeqData('Transcripts-ex1.txt')
        position = 4250000
        # Do sort
        # Sort is done in place, so assignment is not required
        # however the sort function should return a reference to
        # the initial object
        result = rna_sort.sortByClosestEdgeTo(position)
        self.assertEqual(result,rna_sort,
                         "Returned object doesn't match subject")
        # Check that the closest distances are in ascending order
        last_rna_data = None
        for rna_data in rna_sort:
            if not last_rna_data:
                last_rna_data = rna_data
            else:
                self.assertTrue((min(abs(rna_data.getTSS() - position),
                                     abs(rna_data.getTES() - position)) >=
                                 min(abs(last_rna_data.getTSS() - position),
                                     abs(last_rna_data.getTES() - position))),
                                "Sort by closest distance to edge failed")

    def test_sort_by_closest_TSS_to_edge(self):
        rna_sort = RNASeqData('Transcripts-ex1.txt')
        position = (16000000,17500000)
        # Do sort
        # Sort is done in place, so assignment is not required
        # however the sort function should return a reference to
        # the initial object
        result = rna_sort.sortByClosestTSSTo(*position)
        self.assertEqual(result,rna_sort,
                         "Returned object doesn't match subject")
        # Check that the closest distances are in ascending order
        last_rna_data = None
        for rna_data in rna_sort:
            if not last_rna_data:
                last_rna_data = rna_data
            else:
                self.assertTrue((min(abs(rna_data.getTSS() - position[0]),
                                     abs(rna_data.getTSS() - position[1])) >=
                                 min(abs(last_rna_data.getTSS() - position[0]),
                                     abs(last_rna_data.getTSS() - position[1]))),
                                "Sort by closest TSS to edge failed")
        

    def test_reading_bad_file(self):
        rna_bad = RNASeqData('Transcripts-ex2.txt')
        self.assertEqual(len(rna_bad),0,
                         "No transcripts should be read from bad input file")

########################################################################
#
# TestChIPSeqData
#
#########################################################################
class TestChIPSeqData(unittest.TestCase):

    def setUp(self):
        # Create input files for tests
        create_test_file('ChIP_peaks-ex1.txt',chip_peaks_ex1)
        create_test_file('ChIP_peaks-ex2.txt',chip_peaks_ex2)
        create_test_file('ChIP_peaks-ex5.txt',chip_peaks_ex5)

    def tearDown(self):
        # Remove input files
        delete_test_file('ChIP_peaks-ex1.txt')
        delete_test_file('ChIP_peaks-ex2.txt')
        delete_test_file('ChIP_peaks-ex5.txt')

    def test_reading_in_ChIPseq_data(self):
        chip_seq = ChIPSeqData('ChIP_peaks-ex1.txt')
        self.assertEqual(len(chip_seq),5,
                         "Wrong number of lines read from ChIP-seq file")

    def test_is_summit_data(self):
        chip_seq = ChIPSeqData('ChIP_peaks-ex1.txt')
        self.assertTrue(chip_seq.isSummit(),
                        "ChIP data are summits")
        chip_seq = ChIPSeqData('ChIP_peaks-ex5.txt')
        self.assertFalse(chip_seq.isSummit(),
                        "ChIP data are not summits")

    def test_filter_on_chromosome(self):
        chip_chr = ChIPSeqData('ChIP_peaks-ex2.txt')
        chromosome = 'chr2L'
        self.assertEqual(len(chip_chr),10,
                         "Wrong number of lines from ChIP-seq file")
        chip_chr = chip_chr.filterByChr(chromosome)
        self.assertEqual(len(chip_chr), 5,
                         "Wrong number of lines from ChIP-seq after chr filter")
        for chip_data in chip_chr:
            self.assertTrue((chip_data.chr == chromosome),
                            "Wrong chromosome name filtered by chr")

    def test_filter_on_peak_position(self):
        chip_seq = ChIPSeqData('ChIP_peaks-ex1.txt')
        lower,upper = 12000000,15000000
        chip_pos = chip_seq.filterByPosition(upper,lower)
        self.assertEqual(len(chip_pos),2,
                         "Wrong number of peaks filtered")
        for chip_data in chip_pos:
            self.assertTrue((chip_data.start >= lower and
                             chip_data.start <= upper),
                            "Peak should have been filtered out")

    def test_sort_by_distance_from(self):
        chip_sort = ChIPSeqData('ChIP_peaks-ex1.txt')
        position = 12000000
        # Do the sorting
        chip_sort.sortByDistanceFrom(position)
        # Check that each distance is greater than the previous one
        last_chip_data = None
        for chip_data in chip_sort:
            if not last_chip_data:
                last_chip_data = chip_data
            else:
                self.assertTrue((abs(chip_data.start - position) >=
                                 abs(last_chip_data.start - position)),
                                "Sort by distance failed")

    def test_ChIP_peak_inside_region(self):
        chip_seq = ChIPSeqData('ChIP_peaks-ex1.txt')
        lower,upper = 4252000,4254000
        self.assertTrue(chip_seq[0].insideRegion(upper,lower),
                        "ChIP peak should be in region")
        self.assertTrue(chip_seq[0].insideRegion(lower,upper),
                        "ChIP peak should be in region (reversed limits)")
        upper,lower = 4252000,4242000
        self.assertFalse(chip_seq[0].insideRegion(upper,lower),
                         "ChIP peak should not be inside region")
        self.assertFalse(chip_seq[0].insideRegion(lower,upper),
                         "ChIP peak should not be inside region (reversed limits")

    def test_handling_quoted_chr(self):
        chip_seq = ChIPSeqData('ChIP_peaks-ex5.txt')
        self.assertEqual(len(chip_seq),5,
                         "Wrong number of lines read from ChIP-seq file")
        chip_chr = chip_seq.filterByChr("chr4")
        self.assertEqual(len(chip_chr),2,
                         "Wrong number of ChIP-seq records filtered")

########################################################################
#
# TestRNASeqDataWithChIPSeqData
#
#########################################################################
class TestRNASeqDataWithChIPSeqData(unittest.TestCase):

    def setUp(self):
        # Create input files for tests
        create_test_file('Transcripts-ex1.txt',transcripts_ex1)
        create_test_file('ChIP_peaks-ex1.txt',chip_peaks_ex1)

    def tearDown(self):
        # Remove input files
        delete_test_file('Transcripts-ex1.txt')
        delete_test_file('ChIP_peaks-ex1.txt')

    def test_closest_transcript_to_peak(self):
        rna_seq = RNASeqData('Transcripts-ex1.txt')
        rna_data1 = rna_seq[1]
        rna_data2 = rna_seq[2]
        chip_seq = ChIPSeqData('ChIP_peaks-ex1.txt')
        chip_peak = chip_seq[0]
        nearest = GetNearestTranscriptToPeak(rna_data1,rna_data2,chip_peak)
        # 2nd peak should be closer than first
        self.assertEqual(nearest,rna_data1,
                         "Wrong transcript selected as nearest")
        # test when only one is set
        nearest = GetNearestTranscriptToPeak(None,rna_data2,chip_peak)
        self.assertEqual(nearest,rna_data2,
                         "Wrong transcript selected as nearest")

########################################################################
#
# TestAnalyseNearestTSSToSummits
#
#########################################################################
class TestAnalyseNearestTSSToSummits(unittest.TestCase):

    def setUp(self):
        # Create input files for tests
        create_test_file('Transcripts-ex3.txt',transcripts_ex3)
        create_test_file('ChIP_peaks-ex3.txt',chip_peaks_ex3)

    def tearDown(self):
        # Remove input files
        delete_test_file('Transcripts-ex3.txt')
        delete_test_file('ChIP_peaks-ex3.txt')

    def test_AnalyseNearestTSSToSummits(self):
        # Initialise data
        rna_seq = RNASeqData('Transcripts-ex3.txt')
        chip_seq = ChIPSeqData('ChIP_peaks-ex3.txt')
        max_distance = 130000
        # Run the analysis
        results = AnalyseNearestTSSToSummits(chip_seq,
                                             rna_seq,
                                             max_distance)
        # Verify the results are as expected
        self.assertEqual(len(results),
                         nearest_transcript_to_peak_ex3.count('\n'),
                         "Didn't get expected number of results")
        for i in range(len(results)):
            # Get the expected line from the reference output
            expected = nearest_transcript_to_peak_ex3.split('\n')[i+1]
            # Build a line from the actual results
            items = []
            for field in ('chr','start','ID','nearest','TSS',
                          'distance_to_TSS','distance_to_TES',
                          'strand','in_the_gene','transcripts_inbetween',
                          'transcript_ids_inbetween'):
                items.append(str(results[i][field]))
            actual = '\t'.join(items)
            self.assertEqual(actual,expected,
                             "Result differs for line %s\nRef:%s\nAct:%s" % \
                                 (i,expected,actual))

########################################################################
#
# TestAnalyseNearestTranscriptsToPeakEdges
#
#########################################################################
class TestAnalyseNearestTranscriptsToPeakEdges(unittest.TestCase):

    def setUp(self):
        # Create input files for tests
        create_test_file('Transcripts-ex3.txt',transcripts_ex3)
        create_test_file('ChIP_peaks_binding_region-ex3.txt',
                         chip_peaks_binding_region_ex3)

    def tearDown(self):
        # Remove input files
        delete_test_file('Transcripts-ex3.txt')
        delete_test_file('ChIP_peaks_binding_region-ex3.txt')

    def test_AnalyseNearestTranscriptsToPeakEdges(self):
        # Initialise data
        rna_seq = RNASeqData('Transcripts-ex3.txt')
        chip_seq = ChIPSeqData('ChIP_peaks_binding_region-ex3.txt')
        promoter_region = (10000,2500)
        max_closest=4
        # Run the analysis
        results = AnalyseNearestTranscriptsToPeakEdges(chip_seq,
                                                       rna_seq,
                                                       promoter_region,
                                                       max_closest)
        # Verify the results are as expected
        self.assertEqual(len(results),
                         nearest_transcript_to_peak_edge_ex3.count('\n'),
                         "Didn't get expected number of results")
        for i in range(len(results)):
            # Get the expected line from the reference output
            expected = nearest_transcript_to_peak_edge_ex3.split('\n')[i+1]
            # Build a line from the actual results
            items = []
            for field in ('chr','start','end','ID','strand','TSS','TES',
                          'dist_closest_edge','dist_TSS','dist_TES',
                          'overlap_transcript','overlap_promoter'):
                items.append(str(results[i][field]))
            actual = '\t'.join(items)
            self.assertEqual(actual,expected,
                             "Result differs for line %s\nRef:%s\nAct:%s" % \
                                 (i,expected,actual))

    def test_AnalyseNearestTranscriptTSSToPeakEdges(self):
        # Initialise data
        rna_seq = RNASeqData('Transcripts-ex3.txt')
        chip_seq = ChIPSeqData('ChIP_peaks_binding_region-ex3.txt')
        promoter_region = (10000,2500)
        max_closest=4
        # Run the analysis
        results = AnalyseNearestTranscriptsToPeakEdges(chip_seq,
                                                       rna_seq,
                                                       promoter_region,
                                                       max_closest,
                                                       TSS_only=True)
        # Verify the results are as expected
        self.assertEqual(len(results),
                         nearest_transcript_tss_to_peak_edge_ex3.count('\n'),
                         "Didn't get expected number of results")
        for i in range(len(results)):
            # Get the expected line from the reference output
            expected = nearest_transcript_tss_to_peak_edge_ex3.split('\n')[i+1]
            # Build a line from the actual results
            items = []
            for field in ('chr','start','end','ID','strand','TSS','TES',
                          'dist_closest_edge','dist_TSS','dist_TES',
                          'overlap_transcript','overlap_promoter'):
                items.append(str(results[i][field]))
            actual = '\t'.join(items)
            self.assertEqual(actual,expected,
                             "Result differs for line %s\nRef:%s\nAct:%s" % \
                                 (i,expected,actual))

########################################################################
#
# TestAnalyseNearestPeaksToTranscripts
#
#########################################################################
class TestAnalyseNearestPeaksToTranscripts(unittest.TestCase):

    def setUp(self):
        # Create input files for tests
        create_test_file('Transcripts-ex4.txt',transcripts_ex4)
        create_test_file('ChIP_peaks-ex4.txt',chip_peaks_ex4)

    def tearDown(self):
        # Remove input files
        delete_test_file('Transcripts-ex4.txt')
        delete_test_file('ChIP_peaks-ex4.txt')

    def test_AnalyseNearestTranscriptsToPeaks(self):
        # Initialise data
        rna_seq = RNASeqData('Transcripts-ex4.txt')
        chip_seq = ChIPSeqData('ChIP_peaks-ex4.txt')
        window_width = 20000
        # Run the analysis
        results = AnalyseNearestPeaksToTranscripts(rna_seq,
                                                   chip_seq,
                                                   window_width)
        # Verify the results are as expected
        self.assertEqual(len(results),
                         nearest_peak_to_transcript_ex4.count('\n'),
                         "Didn't get expected number of results")
        for i in range(len(results)):
            # Get the expected line from the reference output
            expected = nearest_peak_to_transcript_ex4.split('\n')[i+1]
            # Build a line from the actual results
            items = []
            for field in ('ID','chr_RNA','start','end','strand',
                          'differentially_expressed','number_of_peaks',
                          'chr_ChIP_1','summit_1','distance_1',
                          'chr_ChIP_2','summit_2','distance_2',
                          'chr_ChIP_3','summit_3','distance_3'):
                try:
                    items.append(str(results[i][field]))
                except KeyError:
                    items.append('---')
            actual = '\t'.join(items)
            self.assertEqual(actual,expected,
                             "Result differs for line %s\nRef:'%s'\nAct:'%s'" \
                                 % (i,expected,actual))

########################################################################
#
# TestRegionsOverlap
#
#########################################################################
class TestRegionsOverlap(unittest.TestCase):
    
    def test_regions_dont_overlap(self):
        region1 = (1000,2000)
        region2 = (3000,4000)
        self.assertFalse(regions_overlap(region1,region2),
                         "Regions should not overlap")
        self.assertFalse(regions_overlap(region2,region1),
                         "Regions should not overlap (reversed)")

    def test_regions_partially_overlap(self):
        region1 = (1000,2000)
        region2 = (1500,2500)
        self.assertTrue(regions_overlap(region1,region2),
                        "Regions should overlap")
        self.assertTrue(regions_overlap(region2,region1),
                        "Regions should overlap (reversed)")

    def test_one_region_inside_the_other(self):
        region1 = (1000,2000)
        region2 = (1250,1750)
        self.assertTrue(regions_overlap(region1,region2),
                        "Regions should overlap")
        self.assertTrue(regions_overlap(region2,region1),
                        "Regions should overlap (reversed)")

########################################################################
#
# TestClosestDistanceToRegion
#
#########################################################################
class TestClosestDistanceToRegion(unittest.TestCase):

    def test_get_distance_to_point(self):
        # Single position
        reference = 1240000
        position = 1200000
        # Single reference position
        self.assertEqual(closestDistanceToRegion(reference,position),
                         abs(reference-position),
                         "Incorrect distance #1")
        position = 1300000
        self.assertEqual(closestDistanceToRegion(reference,position),
                         abs(reference-position),
                         "Incorrect distance #2")

    def test_get_closest_distance_outside_region(self):
        # Outside reference region (2 positions)
        reference = 1240000
        position1 = 1100000
        position2 = 1200000
        self.assertEqual(closestDistanceToRegion(reference,position1,position2),
                         abs(reference-position2),
                         "Incorrect closest distance (outside region #1)")
        position1 = 1300000
        position2 = 1400000
        self.assertEqual(closestDistanceToRegion(reference,position1,position2),
                         abs(reference-position1),
                         "Incorrect closest distance (outside region #2)")

    def test_get_closest_distance_inside_region(self):
        # Inside reference region
        reference = 1240000
        position1 = 1200000
        position2 = 1255000
        self.assertEqual(closestDistanceToRegion(reference,position1,position2),
                         abs(reference-position2),
                         "Incorrect closest distance (inside region #1)")
        self.assertEqual(closestDistanceToRegion(reference,position1,position2,
                                                 zero_inside_region=True),
                         0,
                         "Incorrect closest distance (inside region #2)")

########################################################################
#
# Utility functions for tests
#
#########################################################################

def create_test_file(name,data):
    # Writes data to a file to be used in testing
    fp = open(name,'w')
    for line in data.split('\n'):
        fp.write(line+'\n')
    fp.close()

def delete_test_file(name):
    # Deletes named test file
    if os.path.exists(name):
        os.remove(name)

########################################################################
#
# Main: test runner
#
#########################################################################
if __name__ == "__main__":
    # Run tests
    unittest.main()
