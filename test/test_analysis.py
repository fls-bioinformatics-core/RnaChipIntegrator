#
#     test_analysis.py: unit tests for analysis module
#     Copyright (C) University of Manchester 2011-5 Peter Briggs

import unittest
from rnachipintegrator.Features import FeatureSet,Feature
from rnachipintegrator.Peaks import PeakSet,Peak
from rnachipintegrator.analysis import find_nearest_features
from rnachipintegrator.analysis import find_nearest_peaks

# Data
from common import *

summit_data = """chr2L	66811	66812
chr2L	249177	249178
chr2L	605950	605951
chr3L	2258189	2258190
chr3L	19498040	19498041"""

peak_data = """chr2L	66711	66911
chr2L	249077	249277
chr2L	605850	606050
chr3L	2258089	2258290
chr3L	19497940	19498140"""

feature_data = """CG31973	chr2L	25402	59243	-	1
CG2674-RC	chr2L	107926	114433	+	1
CG2674-RE	chr2L	106903	114433	+	0
CG2674-RA	chr2L	107760	114433	+	0
CG3625-RA	chr2L	283385	286528	-	0
CG3625-RB	chr2L	283385	285777	-	1
CG3625-RC	chr2L	283385	291011	-	1
CG2762-RA	chr2L	523467	540542	+	1
CG17941-RA	chr2L	640021	714969	-	0
CG9130-RB	chr3L	1252012	1255989	+	1
CG14448-RA	chr3L	22781539	22782354	+	1
CG13051-RA	chr3L	16257914	16258166	-	0"""

# find_nearest_features
class TestFindNearestFeaturesForSummits(unittest.TestCase):

    def setUp(self):
        # Set up data for tests
        create_test_file('features.txt',feature_data)
        create_test_file('summits.txt',summit_data)
        self.features = FeatureSet('features.txt')
        self.summits = PeakSet('summits.txt')

    def tearDown(self):
        # Remove input files
        delete_test_file('features.txt')
        delete_test_file('summits.txt')

    def test_find_nearest_features_doesnt_change_input_features(self):
        # Check that the input FeatureSet is not altered by the analysis
        features_ref = FeatureSet('features.txt')
        find_nearest_features(self.summits,self.features)
        for f,fref in zip(self.features,features_ref):
            self.assertEqual(f,fref)

    def test_find_nearest_features_doesnt_change_input_peaks(self):
        # Check that the input PeakSet is not altered by the analysis
        summits_ref = PeakSet('summits.txt')
        find_nearest_features(self.summits,self.features)
        for p,pref in zip(self.summits,summits_ref):
            self.assertEqual(p,pref)

    def test_find_nearest_features_summits(self):
        # Run the analysis
        results = list(find_nearest_features(self.summits,
                                             self.features))
        # Correct number of results
        self.assertEqual(len(results),5)
        # Closest distances for 1st peak
        # CG31973      7568
        # CG2674-RE   40092
        # CG2674-RA   40948
        # CG2674-RC   41114
        # CG3625-RB  216573  (218965)
        # CG3625-RA  216573  (219716)
        # CG3625-RC  216573  (224199)
        # CG2762-RA  456655
        # CG17941-RA 573209
        peak,features = results[0]
        self.assertEqual(peak,Peak('chr2L','66811','66812'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG31973')
        self.assertEqual(features[1].id,'CG2674-RE')
        self.assertEqual(features[2].id,'CG2674-RA')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG3625-RB')
        self.assertEqual(features[5].id,'CG3625-RA')
        self.assertEqual(features[6].id,'CG3625-RC')
        self.assertEqual(features[7].id,'CG2762-RA')
        self.assertEqual(features[8].id,'CG17941-RA')
        # Closest distances for 2nd peak
        # CG3625-RB   34207  (36599)
        # CG3625-RA   34207  (37350)
        # CG3625-RC   34207  (41833)
        # CG2674-RC  134744  (141251)
        # CG2674-RA  134744  (141417)
        # CG2674-RE  134744  (142274)
        # CG31973    189934
        # CG2762-RA  274289
        # CG17941-RA 390843
        peak,features = results[1]
        self.assertEqual(peak,Peak('chr2L','249177','249178'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG3625-RB')
        self.assertEqual(features[1].id,'CG3625-RA')
        self.assertEqual(features[2].id,'CG3625-RC')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG2674-RA')
        self.assertEqual(features[5].id,'CG2674-RE')
        self.assertEqual(features[6].id,'CG31973')
        self.assertEqual(features[7].id,'CG2762-RA')
        self.assertEqual(features[8].id,'CG17941-RA')
        # Closest distances for 3rd peak
        # CG17941-RA  34070
        # CG2762-RA   65408
        # CG3625-RC  314939
        # CG3625-RA  319422
        # CG3625-RB  320173
        # CG2674-RC  491517  (498024)
        # CG2674-RA  491517  (498190)
        # CG2674-RE  491517  (499047)
        # CG31973    545707
        peak,features = results[2]
        self.assertEqual(peak,Peak('chr2L','605950','605951'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG17941-RA')
        self.assertEqual(features[1].id,'CG2762-RA')
        self.assertEqual(features[2].id,'CG3625-RC')
        self.assertEqual(features[3].id,'CG3625-RA')
        self.assertEqual(features[4].id,'CG3625-RB')
        self.assertEqual(features[5].id,'CG2674-RC')
        self.assertEqual(features[6].id,'CG2674-RA')
        self.assertEqual(features[7].id,'CG2674-RE')
        self.assertEqual(features[8].id,'CG31973')
        # Closest distances for 4th peak
        # CG9130-RB   1002200
        # CG13051-RA 13999724
        # CG14448-RA 20523349
        peak,features = results[3]
        self.assertEqual(peak,Peak('chr3L','2258189','2258190'))
        self.assertEqual(len(features),3)
        self.assertEqual(features[0].id,'CG9130-RB')
        self.assertEqual(features[1].id,'CG13051-RA')
        self.assertEqual(features[2].id,'CG14448-RA')
        # Closest distances for 5th peak
        # CG13051-RA  3239874
        # CG14448-RA  3283489
        # CG9130-RB  18242051
        peak,features = results[4]
        self.assertEqual(peak,Peak('chr3L','19498040','19498041'))
        self.assertEqual(len(features),3)
        self.assertEqual(features[0].id,'CG13051-RA')
        self.assertEqual(features[1].id,'CG14448-RA')
        self.assertEqual(features[2].id,'CG9130-RB')

    def test_find_nearest_features_summits_tss_only(self):
        # Run the analysis
        results = list(find_nearest_features(self.summits,
                                             self.features,
                                             tss_only=True))
        # Correct number of results
        self.assertEqual(len(results),5)
        # TSS distances for 1st peak
        # CG31973      7568
        # CG2674-RE   40091
        # CG2674-RA   40948
        # CG2674-RC   41114
        # CG3625-RB  218965
        # CG3625-RA  219716
        # CG3625-RC  224199
        # CG2762-RA  456655
        # CG17941-RA 648157
        peak,features = results[0]
        self.assertEqual(peak,Peak('chr2L','66811','66812'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG31973')
        self.assertEqual(features[1].id,'CG2674-RE')
        self.assertEqual(features[2].id,'CG2674-RA')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG3625-RB')
        self.assertEqual(features[5].id,'CG3625-RA')
        self.assertEqual(features[6].id,'CG3625-RC')
        self.assertEqual(features[7].id,'CG2762-RA')
        self.assertEqual(features[8].id,'CG17941-RA')
        # TSS distances for 2nd peak
        # CG3625-RB   36599
        # CG3625-RA   37350
        # CG3625-RC   41833
        # CG2674-RC  141251
        # CG2674-RA  141417
        # CG2674-RE  142274
        # CG31973    189934
        # CG2762-RA  274289
        # CG17941-RA 465791
        peak,features = results[1]
        self.assertEqual(peak,Peak('chr2L','249177','249178'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG3625-RB')
        self.assertEqual(features[1].id,'CG3625-RA')
        self.assertEqual(features[2].id,'CG3625-RC')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG2674-RA')
        self.assertEqual(features[5].id,'CG2674-RE')
        self.assertEqual(features[6].id,'CG31973')
        self.assertEqual(features[7].id,'CG2762-RA')
        self.assertEqual(features[8].id,'CG17941-RA')
        # TSS distances for 3rd peak
        # CG2762-RA   82483
        # CG17941-RA 109018
        # CG3625-RC  314939
        # CG3625-RA  319422
        # CG3625-RB  320173
        # CG2674-RC  498024
        # CG2674-RA  498190
        # CG2674-RE  499047
        # CG31973    545707
        peak,features = results[2]
        self.assertEqual(peak,Peak('chr2L','605950','605951'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG2762-RA')
        self.assertEqual(features[1].id,'CG17941-RA')
        self.assertEqual(features[2].id,'CG3625-RC')
        self.assertEqual(features[3].id,'CG3625-RA')
        self.assertEqual(features[4].id,'CG3625-RB')
        self.assertEqual(features[5].id,'CG2674-RC')
        self.assertEqual(features[6].id,'CG2674-RA')
        self.assertEqual(features[7].id,'CG2674-RE')
        self.assertEqual(features[8].id,'CG31973')
        # TSS distances for 4th peak
        # CG9130-RB   1006177
        # CG13051-RA 13999976
        # CG14448-RA 20523349
        peak,features = results[3]
        self.assertEqual(peak,Peak('chr3L','2258189','2258190'))
        self.assertEqual(len(features),3)
        self.assertEqual(features[0].id,'CG9130-RB')
        self.assertEqual(features[1].id,'CG13051-RA')
        self.assertEqual(features[2].id,'CG14448-RA')
        # TSS distances for 5th peak
        # CG13051-RA  3239874
        # CG14448-RA  3283489
        # CG9130-RB  18246028
        peak,features = results[4]
        self.assertEqual(peak,Peak('chr3L','19498040','19498041'))
        self.assertEqual(len(features),3)
        self.assertEqual(features[0].id,'CG13051-RA')
        self.assertEqual(features[1].id,'CG14448-RA')
        self.assertEqual(features[2].id,'CG9130-RB')

    def test_find_nearest_features_summits_differentially_expressed(self):
        # Run the analysis
        results = list(find_nearest_features(self.summits,
                                             self.features,
                                             only_differentially_expressed=True))
        # Correct number of results
        self.assertEqual(len(results),5)
        # Closest distances for 1st peak, non-DE features omitted
        # CG31973      7568
        # CG2674-RC   41114
        # CG3625-RB  216573  (218965)
        # CG3625-RC  216573  (224199)
        # CG2762-RA  456655
        peak,features = results[0]
        self.assertEqual(peak,Peak('chr2L','66811','66812'))
        self.assertEqual(len(features),5)
        self.assertEqual(features[0].id,'CG31973')
        self.assertEqual(features[1].id,'CG2674-RC')
        self.assertEqual(features[2].id,'CG3625-RB')
        self.assertEqual(features[3].id,'CG3625-RC')
        self.assertEqual(features[4].id,'CG2762-RA')
        # Closest distances for 2nd peak
        # CG3625-RB   34207  (36599)
        # CG3625-RC   34207  (41833)
        # CG2674-RC  134744
        # CG31973    189934
        # CG2762-RA  274289
        peak,features = results[1]
        self.assertEqual(peak,Peak('chr2L','249177','249178'))
        self.assertEqual(len(features),5)
        self.assertEqual(features[0].id,'CG3625-RB')
        self.assertEqual(features[1].id,'CG3625-RC')
        self.assertEqual(features[2].id,'CG2674-RC')
        self.assertEqual(features[3].id,'CG31973')
        self.assertEqual(features[4].id,'CG2762-RA')
        # Closest distances for 3rd peak
        # CG2762-RA   65408
        # CG3625-RC  314939
        # CG3625-RB  320173
        # CG2674-RC  491517
        # CG31973    545707
        peak,features = results[2]
        self.assertEqual(peak,Peak('chr2L','605950','605951'))
        self.assertEqual(len(features),5)
        self.assertEqual(features[0].id,'CG2762-RA')
        self.assertEqual(features[1].id,'CG3625-RC')
        self.assertEqual(features[2].id,'CG3625-RB')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG31973')
        # Closest distances for 4th peak
        # CG9130-RB   1002200
        # CG14448-RA 20523349
        peak,features = results[3]
        self.assertEqual(peak,Peak('chr3L','2258189','2258190'))
        self.assertEqual(len(features),2)
        self.assertEqual(features[0].id,'CG9130-RB')
        self.assertEqual(features[1].id,'CG14448-RA')
        # Closest distances for 5th peak
        # CG14448-RA  3283498
        # CG9130-RB  18242051
        peak,features = results[4]
        self.assertEqual(peak,Peak('chr3L','19498040','19498041'))
        self.assertEqual(len(features),2)
        self.assertEqual(features[0].id,'CG14448-RA')
        self.assertEqual(features[1].id,'CG9130-RB')

    def test_find_nearest_features_summits_distance_cutoff(self):
        # Run the analysis
        results = list(find_nearest_features(self.summits,
                                             self.features,
                                             distance=250000))
        # Correct number of results
        self.assertEqual(len(results),5)
        # Closest distances for 1st peak
        # CG31973      7568
        # CG2674-RE   40092
        # CG2674-RA   40948
        # CG2674-RC   41114
        # CG3625-RB  216573  (218965)
        # CG3625-RA  216573  (219716)
        # CG3625-RC  216573  (224199)
        peak,features = results[0]
        self.assertEqual(peak,Peak('chr2L','66811','66812'))
        self.assertEqual(len(features),7)
        self.assertEqual(features[0].id,'CG31973')
        self.assertEqual(features[1].id,'CG2674-RE')
        self.assertEqual(features[2].id,'CG2674-RA')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG3625-RB')
        self.assertEqual(features[5].id,'CG3625-RA')
        self.assertEqual(features[6].id,'CG3625-RC')
        # Closest distances for 2nd peak
        # CG3625-RB   34207  (36599)
        # CG3625-RA   34207  (37350)
        # CG3625-RC   34207  (41833)
        # CG2674-RC  134744  (141251)
        # CG2674-RA  134744  (141417)
        # CG2674-RE  134744  (142274)
        # CG31973    189934
        peak,features = results[1]
        self.assertEqual(peak,Peak('chr2L','249177','249178'))
        self.assertEqual(len(features),7)
        self.assertEqual(features[0].id,'CG3625-RB')
        self.assertEqual(features[1].id,'CG3625-RA')
        self.assertEqual(features[2].id,'CG3625-RC')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG2674-RA')
        self.assertEqual(features[5].id,'CG2674-RE')
        self.assertEqual(features[6].id,'CG31973')
        # Closest distances for 3rd peak
        # CG17941-RA  34070
        # CG2762-RA   65408
        peak,features = results[2]
        self.assertEqual(peak,Peak('chr2L','605950','605951'))
        self.assertEqual(len(features),2)
        self.assertEqual(features[0].id,'CG17941-RA')
        self.assertEqual(features[1].id,'CG2762-RA')
        # Closest distances for 4th peak
        # No matches
        peak,features = results[3]
        self.assertEqual(peak,Peak('chr3L','2258189','2258190'))
        self.assertEqual(len(features),1)
        self.assertEqual(features[0],None)
        # Closest distances for 5th peak
        # No matches
        peak,features = results[4]
        self.assertEqual(peak,Peak('chr3L','19498040','19498041'))
        self.assertEqual(len(features),1)
        self.assertEqual(features[0],None)

# find_nearest_features
class TestFindNearestFeaturesForRegions(unittest.TestCase):

    def setUp(self):
        # Set up data for tests
        create_test_file('features.txt',feature_data)
        create_test_file('peaks.txt',peak_data)
        self.features = FeatureSet('features.txt')
        self.peaks = PeakSet('peaks.txt')

    def tearDown(self):
        # Remove input files
        delete_test_file('features.txt')
        delete_test_file('peaks.txt')

    def test_find_nearest_features_doesnt_change_input_features(self):
        # Check that the input FeatureSet is not altered by the analysis
        features_ref = FeatureSet('features.txt')
        find_nearest_features(self.peaks,self.features)
        for f,fref in zip(self.features,features_ref):
            self.assertEqual(f,fref)

    def test_find_nearest_features_doesnt_change_input_peaks(self):
        # Check that the input PeakSet is not altered by the analysis
        peaks_ref = PeakSet('peaks.txt')
        find_nearest_features(self.peaks,self.features)
        for p,pref in zip(self.peaks,peaks_ref):
            self.assertEqual(p,pref)

    def test_find_nearest_features_regions(self):
        # Run the analysis
        results = list(find_nearest_features(self.peaks,
                                             self.features))
        # Correct number of results
        self.assertEqual(len(results),5)
        # Closest distances for 1st peak
        # CG31973      7468
        # CG2674-RE   39992
        # CG2674-RA   40849
        # CG2674-RC   41015
        # CG3625-RB  216474  (218866)
        # CG3625-RA  216474  (219617)
        # CG3625-RC  216474  (224100)
        # CG2762-RA  456556
        # CG17941-RA 573110
        peak,features = results[0]
        self.assertEqual(peak,Peak('chr2L','66711','66911'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG31973')
        self.assertEqual(features[1].id,'CG2674-RE')
        self.assertEqual(features[2].id,'CG2674-RA')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG3625-RB')
        self.assertEqual(features[5].id,'CG3625-RA')
        self.assertEqual(features[6].id,'CG3625-RC')
        self.assertEqual(features[7].id,'CG2762-RA')
        self.assertEqual(features[8].id,'CG17941-RA')
        # Closest distances for 2nd peak
        # CG3625-RB   34108  (36500)
        # CG3625-RA   34108  (37251)
        # CG3625-RC   34108  (41734)
        # CG2674-RC  134644  (141151)
        # CG2674-RA  134644  (141317)
        # CG2674-RE  134644  (142174)
        # CG31973    189834
        # CG2762-RA  274190
        # CG17941-RA 390744
        peak,features = results[1]
        self.assertEqual(peak,Peak('chr2L','249077','249277'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG3625-RB')
        self.assertEqual(features[1].id,'CG3625-RA')
        self.assertEqual(features[2].id,'CG3625-RC')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG2674-RA')
        self.assertEqual(features[5].id,'CG2674-RE')
        self.assertEqual(features[6].id,'CG31973')
        self.assertEqual(features[7].id,'CG2762-RA')
        self.assertEqual(features[8].id,'CG17941-RA')
        # Closest distances for 3rd peak
        # CG17941-RA  33971
        # CG2762-RA   65308
        # CG3625-RC  314839
        # CG3625-RA  319322
        # CG3625-RB  320073
        # CG2674-RC  491417  (497924)
        # CG2674-RA  491417  (498090)
        # CG2674-RE  491417  (498947)
        # CG31973    546607
        peak,features = results[2]
        self.assertEqual(peak,Peak('chr2L','605850','606050'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG17941-RA')
        self.assertEqual(features[1].id,'CG2762-RA')
        self.assertEqual(features[2].id,'CG3625-RC')
        self.assertEqual(features[3].id,'CG3625-RA')
        self.assertEqual(features[4].id,'CG3625-RB')
        self.assertEqual(features[5].id,'CG2674-RC')
        self.assertEqual(features[6].id,'CG2674-RA')
        self.assertEqual(features[7].id,'CG2674-RE')
        self.assertEqual(features[8].id,'CG31973')
        # Closest distances for 4th peak
        # CG9130-RB   1002100
        # CG13051-RA 13999624
        # CG14448-RA 20523249
        peak,features = results[3]
        self.assertEqual(peak,Peak('chr3L','2258089','2258290'))
        self.assertEqual(len(features),3)
        self.assertEqual(features[0].id,'CG9130-RB')
        self.assertEqual(features[1].id,'CG13051-RA')
        self.assertEqual(features[2].id,'CG14448-RA')
        # Closest distances for 5th peak
        # CG13051-RA  3239774
        # CG14448-RA  3283399
        # CG9130-RB  18241951
        peak,features = results[4]
        self.assertEqual(peak,Peak('chr3L','19497940','19498140'))
        self.assertEqual(len(features),3)
        self.assertEqual(features[0].id,'CG13051-RA')
        self.assertEqual(features[1].id,'CG14448-RA')
        self.assertEqual(features[2].id,'CG9130-RB')

    def test_find_nearest_features_regions_tss_only(self):
        # Run the analysis
        results = list(find_nearest_features(self.peaks,
                                             self.features,
                                             tss_only=True))
        # Correct number of results
        self.assertEqual(len(results),5)
        # Closest distances for 1st peak
        # CG31973      7468
        # CG2674-RE   39992
        # CG2674-RA   40849
        # CG2674-RC   41015
        # CG3625-RB  218866
        # CG3625-RA  219617
        # CG3625-RC  224100
        # CG2762-RA  456556
        # CG17941-RA 648058
        peak,features = results[0]
        self.assertEqual(peak,Peak('chr2L','66711','66911'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG31973')
        self.assertEqual(features[1].id,'CG2674-RE')
        self.assertEqual(features[2].id,'CG2674-RA')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG3625-RB')
        self.assertEqual(features[5].id,'CG3625-RA')
        self.assertEqual(features[6].id,'CG3625-RC')
        self.assertEqual(features[7].id,'CG2762-RA')
        self.assertEqual(features[8].id,'CG17941-RA')
        # Closest distances for 2nd peak
        # CG3625-RB   36500
        # CG3625-RA   37251
        # CG3625-RC   41734
        # CG2674-RC  141151
        # CG2674-RA  141317
        # CG2674-RE  142174
        # CG31973    189834
        # CG2762-RA  274190
        # CG17941-RA 465692
        peak,features = results[1]
        self.assertEqual(peak,Peak('chr2L','249077','249277'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG3625-RB')
        self.assertEqual(features[1].id,'CG3625-RA')
        self.assertEqual(features[2].id,'CG3625-RC')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG2674-RA')
        self.assertEqual(features[5].id,'CG2674-RE')
        self.assertEqual(features[6].id,'CG31973')
        self.assertEqual(features[7].id,'CG2762-RA')
        self.assertEqual(features[8].id,'CG17941-RA')
        # Closest distances for 3rd peak
        # CG2762-RA   82382
        # CG17941-RA 108919
        # CG3625-RC  314839
        # CG3625-RA  319322
        # CG3625-RB  320073
        # CG2674-RC  497924
        # CG2674-RA  498090
        # CG2674-RE  498947
        # CG31973    546607
        peak,features = results[2]
        self.assertEqual(peak,Peak('chr2L','605850','606050'))
        self.assertEqual(len(features),9)
        self.assertEqual(features[0].id,'CG2762-RA')
        self.assertEqual(features[1].id,'CG17941-RA')
        self.assertEqual(features[2].id,'CG3625-RC')
        self.assertEqual(features[3].id,'CG3625-RA')
        self.assertEqual(features[4].id,'CG3625-RB')
        self.assertEqual(features[5].id,'CG2674-RC')
        self.assertEqual(features[6].id,'CG2674-RA')
        self.assertEqual(features[7].id,'CG2674-RE')
        self.assertEqual(features[8].id,'CG31973')
        # Closest distances for 4th peak
        # CG9130-RB   1006077
        # CG13051-RA 13999876
        # CG14448-RA 20523249
        peak,features = results[3]
        self.assertEqual(peak,Peak('chr3L','2258089','2258290'))
        self.assertEqual(len(features),3)
        self.assertEqual(features[0].id,'CG9130-RB')
        self.assertEqual(features[1].id,'CG13051-RA')
        self.assertEqual(features[2].id,'CG14448-RA')
        # Closest distances for 5th peak
        # CG13051-RA  3239774
        # CG14448-RA  3283399
        # CG9130-RB  18245928
        peak,features = results[4]
        self.assertEqual(peak,Peak('chr3L','19497940','19498140'))
        self.assertEqual(len(features),3)
        self.assertEqual(features[0].id,'CG13051-RA')
        self.assertEqual(features[1].id,'CG14448-RA')
        self.assertEqual(features[2].id,'CG9130-RB')

    def test_find_nearest_features_regions_differentially_expressed(self):
        # Run the analysis
        results = list(find_nearest_features(self.peaks,
                                             self.features,
                                             only_differentially_expressed=True))
        # Correct number of results
        self.assertEqual(len(results),5)
        # Closest distances for 1st peak
        # CG31973      7468
        # CG2674-RC   41015
        # CG3625-RB  216474  (218866)
        # CG3625-RC  216474  (224100)
        # CG2762-RA  456556
        peak,features = results[0]
        self.assertEqual(peak,Peak('chr2L','66711','66911'))
        self.assertEqual(len(features),5)
        self.assertEqual(features[0].id,'CG31973')
        self.assertEqual(features[1].id,'CG2674-RC')
        self.assertEqual(features[2].id,'CG3625-RB')
        self.assertEqual(features[3].id,'CG3625-RC')
        self.assertEqual(features[4].id,'CG2762-RA')
        # Closest distances for 2nd peak
        # CG3625-RB   34108  (36500)
        # CG3625-RC   34108  (41734)
        # CG2674-RC  134644
        # CG31973    189834
        # CG2762-RA  274190
        peak,features = results[1]
        self.assertEqual(peak,Peak('chr2L','249077','249277'))
        self.assertEqual(len(features),5)
        self.assertEqual(features[0].id,'CG3625-RB')
        self.assertEqual(features[1].id,'CG3625-RC')
        self.assertEqual(features[2].id,'CG2674-RC')
        self.assertEqual(features[3].id,'CG31973')
        self.assertEqual(features[4].id,'CG2762-RA')
        # Closest distances for 3rd peak
        # CG2762-RA   65308
        # CG3625-RC  314839
        # CG3625-RB  320073
        # CG2674-RC  491417
        # CG31973    546607
        peak,features = results[2]
        self.assertEqual(peak,Peak('chr2L','605850','606050'))
        self.assertEqual(len(features),5)
        self.assertEqual(features[0].id,'CG2762-RA')
        self.assertEqual(features[1].id,'CG3625-RC')
        self.assertEqual(features[2].id,'CG3625-RB')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG31973')
        # Closest distances for 4th peak
        # CG9130-RB   1002100
        # CG14448-RA 20523249
        peak,features = results[3]
        self.assertEqual(peak,Peak('chr3L','2258089','2258290'))
        self.assertEqual(len(features),2)
        self.assertEqual(features[0].id,'CG9130-RB')
        self.assertEqual(features[1].id,'CG14448-RA')
        # Closest distances for 5th peak
        # CG14448-RA  3283399
        # CG9130-RB  18241951
        peak,features = results[4]
        self.assertEqual(peak,Peak('chr3L','19497940','19498140'))
        self.assertEqual(len(features),2)
        self.assertEqual(features[0].id,'CG14448-RA')
        self.assertEqual(features[1].id,'CG9130-RB')

    def test_find_nearest_features_regions_distance_cutoff(self):
        # Run the analysis
        results = list(find_nearest_features(self.peaks,
                                             self.features,
                                             distance=250000))
        # Correct number of results
        self.assertEqual(len(results),5)
        # Closest distances for 1st peak
        # CG31973      7468
        # CG2674-RE   39992
        # CG2674-RA   40849
        # CG2674-RC   41015
        # CG3625-RB  216474  (218866)
        # CG3625-RA  216474  (219617)
        # CG3625-RC  216474  (224100)
        peak,features = results[0]
        self.assertEqual(peak,Peak('chr2L','66711','66911'))
        self.assertEqual(len(features),7)
        self.assertEqual(features[0].id,'CG31973')
        self.assertEqual(features[1].id,'CG2674-RE')
        self.assertEqual(features[2].id,'CG2674-RA')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG3625-RB')
        self.assertEqual(features[5].id,'CG3625-RA')
        self.assertEqual(features[6].id,'CG3625-RC')
        # Closest distances for 2nd peak
        # CG3625-RB   34108  (36500)
        # CG3625-RA   34108  (37251)
        # CG3625-RC   34108  (41734)
        # CG2674-RC  134644  (141151)
        # CG2674-RA  134644  (141317)
        # CG2674-RE  134644  (142174)
        # CG31973    189834
        peak,features = results[1]
        self.assertEqual(peak,Peak('chr2L','249077','249277'))
        self.assertEqual(len(features),7)
        self.assertEqual(features[0].id,'CG3625-RB')
        self.assertEqual(features[1].id,'CG3625-RA')
        self.assertEqual(features[2].id,'CG3625-RC')
        self.assertEqual(features[3].id,'CG2674-RC')
        self.assertEqual(features[4].id,'CG2674-RA')
        self.assertEqual(features[5].id,'CG2674-RE')
        self.assertEqual(features[6].id,'CG31973')
        # Closest distances for 3rd peak
        # CG17941-RA  33971
        # CG2762-RA   65308
        # CG3625-RC  314839
        # CG3625-RA  319322
        # CG3625-RB  320073
        # CG2674-RC  491417  (497924)
        # CG2674-RA  491417  (498090)
        # CG2674-RE  491417  (498947)
        # CG31973    546607
        peak,features = results[2]
        self.assertEqual(peak,Peak('chr2L','605850','606050'))
        self.assertEqual(len(features),2)
        self.assertEqual(features[0].id,'CG17941-RA')
        self.assertEqual(features[1].id,'CG2762-RA')
        # Closest distances for 4th peak
        # CG9130-RB   1002100
        # CG13051-RA 13999624
        # CG14448-RA 20523249
        peak,features = results[3]
        self.assertEqual(peak,Peak('chr3L','2258089','2258290'))
        self.assertEqual(len(features),1)
        self.assertEqual(features[0],None)
        # Closest distances for 5th peak
        # CG13051-RA  3239774
        # CG14448-RA  3283399
        # CG9130-RB  18241951
        peak,features = results[4]
        self.assertEqual(peak,Peak('chr3L','19497940','19498140'))
        self.assertEqual(len(features),1)
        self.assertEqual(features[0],None)

# find_nearest_peaks
class TestFindNearestPeaksForSummits(unittest.TestCase):

    def setUp(self):
        # Set up data for tests
        create_test_file('features.txt',feature_data)
        create_test_file('summits.txt',summit_data)
        self.features = FeatureSet('features.txt')
        self.summits = PeakSet('summits.txt')

    def tearDown(self):
        # Remove input files
        delete_test_file('features.txt')
        delete_test_file('summits.txt')

    def test_find_nearest_peaks_summits(self):
        # Run the analysis
        results = list(find_nearest_peaks(self.features,
                                          self.summits))
        # Convenience variables for summits
        summits = (Peak('chr2L','66811','66812'),
                   Peak('chr2L','249177','249178'),
                   Peak('chr2L','605950','605951'))
        # Correct number of results
        self.assertEqual(len(results),12)
        # Closest distances for 1st feature
        # #0    7568
        # #1  189934
        # #2  546707
        feature,peaks = results[0]
        self.assertEqual(feature,Feature('CG31973',
                                         'chr2L','25402','59243','-'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 2nd feature
        # #0   41114
        # #1  134744
        # #2  491517
        feature,peaks = results[1]
        self.assertEqual(feature,Feature('CG2674-RC',
                                         'chr2L','107926','114433','+'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 3rd feature
        # #0   40092
        # #1  134744
        # #2  491517
        feature,peaks = results[2]
        self.assertEqual(feature,Feature('CG2674-RE',
                                         'chr2L','106903','114433','+'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 4th feature
        # #0   40948
        # #1  134744
        # #2  491517
        feature,peaks = results[3]
        self.assertEqual(feature,Feature('CG2674-RA',
                                         'chr2L','107760','114433','+'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 5th feature
        # #1   34207
        # #0  216573
        # #2  319422
        feature,peaks = results[4]
        self.assertEqual(feature,Feature('CG3625-RA',
                                         'chr2L','283385','286528','-'))
        expected = PeakSet()
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 6th feature
        # #1   34207
        # #0  216573
        # #2  320173
        feature,peaks = results[5]
        self.assertEqual(feature,Feature('CG3625-RB',
                                         'chr2L','283385','285777','-'))
        expected = PeakSet()
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 7th feature
        # #1    34207
        # #0   216573
        # #2   314939
        feature,peaks = results[6]
        self.assertEqual(feature,Feature('CG3625-RC',
                                         'chr2L','283385','291011','-'))
        expected = PeakSet()
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 8th feature
        # #2    65408
        # #1   274289
        # #0   456665
        feature,peaks = results[7]
        self.assertEqual(feature,Feature('CG2762-RA',
                                         'chr2L','523467','540542','+'))
        expected = PeakSet()
        expected.addPeak(summits[2])
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        self.assertEqual(peaks,expected)
        # Closest distances for 9th feature
        # #2    34070
        # #1   390843
        # #0   573209
        feature,peaks = results[8]
        self.assertEqual(feature,Feature('CG17941-RA',
                                         'chr2L','640021','714969','-'))
        expected = PeakSet()
        expected.addPeak(summits[2])
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        self.assertEqual(peaks,expected)

    def test_find_nearest_peaks_summits_tss_only(self):
        # Run the analysis
        results = list(find_nearest_peaks(self.features,
                                          self.summits,
                                          tss_only=True))
        # Convenience variables for summits
        summits = (Peak('chr2L','66811','66812'),
                   Peak('chr2L','249177','249178'),
                   Peak('chr2L','605950','605951'))
        # Correct number of results
        self.assertEqual(len(results),12)
        # Closest distances for 1st feature
        # #0    7568
        # #1  189934
        # #2  546707
        feature,peaks = results[0]
        self.assertEqual(feature,Feature('CG31973',
                                         'chr2L','25402','59243','-'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 2nd feature
        # #0   41114
        # #1  141251
        # #2  498024
        feature,peaks = results[1]
        self.assertEqual(feature,Feature('CG2674-RC',
                                         'chr2L','107926','114433','+'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 3rd feature
        # #0   40091
        # #1  142274
        # #2  499047
        feature,peaks = results[2]
        self.assertEqual(feature,Feature('CG2674-RE',
                                         'chr2L','106903','114433','+'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 4th feature
        # #0   40948
        # #1  141417
        # #2  498190
        feature,peaks = results[3]
        self.assertEqual(feature,Feature('CG2674-RA',
                                         'chr2L','107760','114433','+'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 5th feature
        # #1   37251
        # #0  219617
        # #2  319322
        feature,peaks = results[4]
        self.assertEqual(feature,Feature('CG3625-RA',
                                         'chr2L','283385','286528','-'))
        expected = PeakSet()
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 6th feature
        # #1   36500
        # #0  218866
        # #2  320073
        feature,peaks = results[5]
        self.assertEqual(feature,Feature('CG3625-RB',
                                         'chr2L','283385','285777','-'))
        expected = PeakSet()
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 7th feature
        # #1    41734
        # #0   224100
        # #2   314839
        feature,peaks = results[6]
        self.assertEqual(feature,Feature('CG3625-RC',
                                         'chr2L','283385','291011','-'))
        expected = PeakSet()
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 8th feature
        # #2    82382
        # #1   27490
        # #0   456665
        feature,peaks = results[7]
        self.assertEqual(feature,Feature('CG2762-RA',
                                         'chr2L','523467','540542','+'))
        expected = PeakSet()
        expected.addPeak(summits[2])
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        self.assertEqual(peaks,expected)
        # Closest distances for 9th feature
        # #2   108919
        # #1   465692
        # #0   648058
        feature,peaks = results[8]
        self.assertEqual(feature,Feature('CG17941-RA',
                                         'chr2L','640021','714969','-'))
        expected = PeakSet()
        expected.addPeak(summits[2])
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        self.assertEqual(peaks,expected)

    def test_find_nearest_peaks_summits_differentially_expressed(self):
        # Run the analysis
        results = list(find_nearest_peaks(self.features,
                                          self.summits,
                                          only_differentially_expressed=True))
        # Convenience variables for summits
        summits = (Peak('chr2L','66811','66812'),
                   Peak('chr2L','249177','249178'),
                   Peak('chr2L','605950','605951'))
        # Correct number of results
        self.assertEqual(len(results),7)
        # Closest distances for 1st feature
        feature,peaks = results[0]
        self.assertEqual(feature,Feature('CG31973',
                                         'chr2L','25402','59243','-'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 2nd feature
        feature,peaks = results[1]
        self.assertEqual(feature,Feature('CG2674-RC',
                                         'chr2L','107926','114433','+'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 3rd feature
        feature,peaks = results[2]
        self.assertEqual(feature,Feature('CG3625-RB',
                                         'chr2L','283385','285777','-'))
        expected = PeakSet()
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 4th feature
        feature,peaks = results[3]
        self.assertEqual(feature,Feature('CG3625-RC',
                                         'chr2L','283385','291011','-'))
        expected = PeakSet()
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 5th feature
        feature,peaks = results[4]
        self.assertEqual(feature,Feature('CG2762-RA',
                                         'chr2L','523467','540542','+'))
        expected = PeakSet()
        expected.addPeak(summits[2])
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        self.assertEqual(peaks,expected)

    def test_find_nearest_peaks_summits_distance_cutoff(self):
        # Run the analysis
        results = list(find_nearest_peaks(self.features,
                                          self.summits,
                                          distance=250000))
        # Convenience variables for summits
        summits = (Peak('chr2L','66811','66812'),
                   Peak('chr2L','249177','249178'),
                   Peak('chr2L','605950','605951'))
        # Correct number of results
        self.assertEqual(len(results),12)
        # Closest distances for 1st feature
        # #0    7568
        # #1  189934
        feature,peaks = results[0]
        self.assertEqual(feature,Feature('CG31973',
                                         'chr2L','25402','59243','-'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        self.assertEqual(peaks,expected)
        # Closest distances for 2nd feature
        # #0   41114
        # #1  134744
        feature,peaks = results[1]
        self.assertEqual(feature,Feature('CG2674-RC',
                                         'chr2L','107926','114433','+'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        self.assertEqual(peaks,expected)
        # Closest distances for 3rd feature
        # #0   40092
        # #1  134744
        feature,peaks = results[2]
        self.assertEqual(feature,Feature('CG2674-RE',
                                         'chr2L','106903','114433','+'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        self.assertEqual(peaks,expected)
        # Closest distances for 4th feature
        # #0   40948
        # #1  134744
        # #2  491517
        feature,peaks = results[3]
        self.assertEqual(feature,Feature('CG2674-RA',
                                         'chr2L','107760','114433','+'))
        expected = PeakSet()
        expected.addPeak(summits[0])
        expected.addPeak(summits[1])
        self.assertEqual(peaks,expected)
        # Closest distances for 5th feature
        # #1   34207
        # #0  216573
        feature,peaks = results[4]
        self.assertEqual(feature,Feature('CG3625-RA',
                                         'chr2L','283385','286528','-'))
        expected = PeakSet()
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        self.assertEqual(peaks,expected)
        # Closest distances for 6th feature
        # #1   34207
        # #0  216573
        feature,peaks = results[5]
        self.assertEqual(feature,Feature('CG3625-RB',
                                         'chr2L','283385','285777','-'))
        expected = PeakSet()
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        self.assertEqual(peaks,expected)
        # Closest distances for 7th feature
        # #1    34207
        # #0   216573
        feature,peaks = results[6]
        self.assertEqual(feature,Feature('CG3625-RC',
                                         'chr2L','283385','291011','-'))
        expected = PeakSet()
        expected.addPeak(summits[1])
        expected.addPeak(summits[0])
        self.assertEqual(peaks,expected)
        # Closest distances for 8th feature
        # #2    65408
        feature,peaks = results[7]
        self.assertEqual(feature,Feature('CG2762-RA',
                                         'chr2L','523467','540542','+'))
        expected = PeakSet()
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 9th feature
        # #2    34070
        feature,peaks = results[8]
        self.assertEqual(feature,Feature('CG17941-RA',
                                         'chr2L','640021','714969','-'))
        expected = PeakSet()
        expected.addPeak(summits[2])
        self.assertEqual(peaks,expected)

class TestFindNearestPeaksForRegions(unittest.TestCase):

    def setUp(self):
        # Set up data for tests
        create_test_file('features.txt',feature_data)
        create_test_file('peaks.txt',peak_data)
        self.features = FeatureSet('features.txt')
        self.peaks = PeakSet('peaks.txt')

    def tearDown(self):
        # Remove input files
        delete_test_file('features.txt')
        delete_test_file('peaks.txt')

    def test_find_nearest_peaks_regions(self):
        # Run the analysis
        results = list(find_nearest_peaks(self.features,
                                          self.peaks))
        # Convenience variables for summits
        peaks_ = (Peak('chr2L','66711','66911'),
                  Peak('chr2L','249077','249277'),
                  Peak('chr2L','605850','606050'))
        # Correct number of results
        self.assertEqual(len(results),12)
        # Closest distances for 1st feature
        # #0    7568
        # #1  189834
        # #2  546607
        feature,peaks = results[0]
        self.assertEqual(feature,Feature('CG31973',
                                         'chr2L','25402','59243','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 2nd feature
        # #0   41015
        # #1  134644
        # #2  491417
        feature,peaks = results[1]
        self.assertEqual(feature,Feature('CG2674-RC',
                                         'chr2L','107926','114433','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 3rd feature
        # #0   40092
        # #1  134744
        # #2  491517
        feature,peaks = results[2]
        self.assertEqual(feature,Feature('CG2674-RE',
                                         'chr2L','106903','114433','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 4th feature
        # #0   39992
        # #1  134644
        # #2  491417
        feature,peaks = results[3]
        self.assertEqual(feature,Feature('CG2674-RA',
                                         'chr2L','107760','114433','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 5th feature
        # #1   34108
        # #0  216474
        # #2  319322
        feature,peaks = results[4]
        self.assertEqual(feature,Feature('CG3625-RA',
                                         'chr2L','283385','286528','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 6th feature
        # #1   34108
        # #0  216474
        # #2  320073
        feature,peaks = results[5]
        self.assertEqual(feature,Feature('CG3625-RB',
                                         'chr2L','283385','285777','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 7th feature
        # #1    34108
        # #0   216474
        # #2   314839
        feature,peaks = results[6]
        self.assertEqual(feature,Feature('CG3625-RC',
                                         'chr2L','283385','291011','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 8th feature
        # #2    65308
        # #1   274190
        # #0   456556
        feature,peaks = results[7]
        self.assertEqual(feature,Feature('CG2762-RA',
                                         'chr2L','523467','540542','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[2])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        self.assertEqual(peaks,expected)
        # Closest distances for 9th feature
        # #2    33971
        # #1   390744
        # #0   573110
        feature,peaks = results[8]
        self.assertEqual(feature,Feature('CG17941-RA',
                                         'chr2L','640021','714969','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[2])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        self.assertEqual(peaks,expected)

    def test_find_nearest_peaks_regions_tss_only(self):
        # Run the analysis
        results = list(find_nearest_peaks(self.features,
                                          self.peaks,
                                          tss_only=True))
        # Convenience variables for summits
        peaks_ = (Peak('chr2L','66711','66911'),
                  Peak('chr2L','249077','249277'),
                  Peak('chr2L','605850','606050'))
        # Correct number of results
        self.assertEqual(len(results),12)
        # Closest TSS distances for 1st feature
        # #0    7568
        # #1  189834
        # #2  546607
        feature,peaks = results[0]
        self.assertEqual(feature,Feature('CG31973',
                                         'chr2L','25402','59243','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 2nd feature
        # #0   41015
        # #1  141151
        # #2  497924
        feature,peaks = results[1]
        self.assertEqual(feature,Feature('CG2674-RC',
                                         'chr2L','107926','114433','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 3rd feature
        # #0   39992
        # #1  142174
        # #2  498947
        feature,peaks = results[2]
        self.assertEqual(feature,Feature('CG2674-RE',
                                         'chr2L','106903','114433','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 4th feature
        # #0   40849
        # #1  141317
        # #2  498090
        feature,peaks = results[3]
        self.assertEqual(feature,Feature('CG2674-RA',
                                         'chr2L','107760','114433','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 5th feature
        # #1   37251
        # #0  219617
        # #2  319322
        feature,peaks = results[4]
        self.assertEqual(feature,Feature('CG3625-RA',
                                         'chr2L','283385','286528','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 6th feature
        # #1   36500
        # #0  218866
        # #2  320073
        feature,peaks = results[5]
        self.assertEqual(feature,Feature('CG3625-RB',
                                         'chr2L','283385','285777','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 7th feature
        # #1    41734
        # #0   224100
        # #2   314839
        feature,peaks = results[6]
        self.assertEqual(feature,Feature('CG3625-RC',
                                         'chr2L','283385','291011','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 8th feature
        # #2    82383
        # #1   274190
        # #0   456556
        feature,peaks = results[7]
        self.assertEqual(feature,Feature('CG2762-RA',
                                         'chr2L','523467','540542','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[2])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        self.assertEqual(peaks,expected)
        # Closest distances for 9th feature
        # #2   108919
        # #1   465692
        # #0   648058
        feature,peaks = results[8]
        self.assertEqual(feature,Feature('CG17941-RA',
                                         'chr2L','640021','714969','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[2])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        self.assertEqual(peaks,expected)

    def test_find_nearest_peaks_regions_differentially_expressed(self):
        # Run the analysis
        results = list(find_nearest_peaks(self.features,
                                          self.peaks,
                                          only_differentially_expressed=True))
        # Convenience variables for summits
        peaks_ = (Peak('chr2L','66711','66911'),
                  Peak('chr2L','249077','249277'),
                  Peak('chr2L','605850','606050'))
        # Correct number of results
        self.assertEqual(len(results),7)
        # Closest distances for 1st feature
        # #0    7568
        # #1  189834
        # #2  546607
        feature,peaks = results[0]
        self.assertEqual(feature,Feature('CG31973',
                                         'chr2L','25402','59243','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 2nd feature
        # #0   41015
        # #1  134644
        # #2  491417
        feature,peaks = results[1]
        self.assertEqual(feature,Feature('CG2674-RC',
                                         'chr2L','107926','114433','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 3rd feature
        # #1   34108
        # #0  216474
        # #2  320073
        feature,peaks = results[2]
        self.assertEqual(feature,Feature('CG3625-RB',
                                         'chr2L','283385','285777','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 4th feature
        # #1    34108
        # #0   216474
        # #2   314839
        feature,peaks = results[3]
        self.assertEqual(feature,Feature('CG3625-RC',
                                         'chr2L','283385','291011','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 5th feature
        # #2    65308
        # #1   274190
        # #0   456556
        feature,peaks = results[4]
        self.assertEqual(feature,Feature('CG2762-RA',
                                         'chr2L','523467','540542','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[2])
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        self.assertEqual(peaks,expected)

    def test_find_nearest_peaks_regions_distance_cutoff(self):
        # Run the analysis
        results = list(find_nearest_peaks(self.features,
                                          self.peaks,
                                          distance=250000))
        # Convenience variables for summits
        peaks_ = (Peak('chr2L','66711','66911'),
                  Peak('chr2L','249077','249277'),
                  Peak('chr2L','605850','606050'))
        # Correct number of results
        self.assertEqual(len(results),12)
        # Closest distances for 1st feature
        # #0    7568
        # #1  189834
        feature,peaks = results[0]
        self.assertEqual(feature,Feature('CG31973',
                                         'chr2L','25402','59243','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        self.assertEqual(peaks,expected)
        # Closest distances for 2nd feature
        # #0   41015
        # #1  134644
        feature,peaks = results[1]
        self.assertEqual(feature,Feature('CG2674-RC',
                                         'chr2L','107926','114433','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        self.assertEqual(peaks,expected)
        # Closest distances for 3rd feature
        # #0   40092
        # #1  134744
        feature,peaks = results[2]
        self.assertEqual(feature,Feature('CG2674-RE',
                                         'chr2L','106903','114433','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        self.assertEqual(peaks,expected)
        # Closest distances for 4th feature
        # #0   39992
        # #1  134644
        feature,peaks = results[3]
        self.assertEqual(feature,Feature('CG2674-RA',
                                         'chr2L','107760','114433','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[0])
        expected.addPeak(peaks_[1])
        self.assertEqual(peaks,expected)
        # Closest distances for 5th feature
        # #1   34108
        # #0  216474
        feature,peaks = results[4]
        self.assertEqual(feature,Feature('CG3625-RA',
                                         'chr2L','283385','286528','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        self.assertEqual(peaks,expected)
        # Closest distances for 6th feature
        # #1   34108
        # #0  216474
        feature,peaks = results[5]
        self.assertEqual(feature,Feature('CG3625-RB',
                                         'chr2L','283385','285777','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        self.assertEqual(peaks,expected)
        # Closest distances for 7th feature
        # #1    34108
        # #0   216474
        feature,peaks = results[6]
        self.assertEqual(feature,Feature('CG3625-RC',
                                         'chr2L','283385','291011','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[1])
        expected.addPeak(peaks_[0])
        self.assertEqual(peaks,expected)
        # Closest distances for 8th feature
        # #2    65308
        feature,peaks = results[7]
        self.assertEqual(feature,Feature('CG2762-RA',
                                         'chr2L','523467','540542','+'))
        expected = PeakSet()
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
        # Closest distances for 9th feature
        # #2    33971
        feature,peaks = results[8]
        self.assertEqual(feature,Feature('CG17941-RA',
                                         'chr2L','640021','714969','-'))
        expected = PeakSet()
        expected.addPeak(peaks_[2])
        self.assertEqual(peaks,expected)
