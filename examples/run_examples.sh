#!/bin/sh -e
#
# Run examples for RnaChipIntegrator
#
TEST_DIR=$(pwd)
REF_DATA=$TEST_DIR/ref-data
if [ ! -d test-output ] ; then
    mkdir test-output
else
    rm -f test-output/*
fi
cd test-output
#
# Summits
SUMMIT_OUTPUTS="summits_features_per_peak.txt summits_peaks_per_feature.txt"
RnaChipIntegrator --name=summits \
    --cutoff=130000 \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_summits.txt
for f in $SUMMIT_OUTPUTS ; do
    diff -f $REF_DATA/ref_$f $f
done
# Regions
REGION_OUTPUTS="regions_features_per_peak.txt regions_peaks_per_feature.txt"
RnaChipIntegrator --name=regions \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_regions.txt
for f in $REGION_OUTPUTS ; do
    diff -f $REF_DATA/ref_$f $f
done
##
#
