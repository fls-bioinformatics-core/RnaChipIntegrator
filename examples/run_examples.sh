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
SUMMIT_OUTPUTS="summits_PeaksToTranscripts.txt summits_TSSToSummits.txt"
RnaChipIntegrator.py --project=summits \
    --window=130000 \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_summits.txt
for f in $SUMMIT_OUTPUTS ; do
    diff -f $REF_DATA/ref_$f $f
done
# Regions
REGION_OUTPUTS="regions_TranscriptsToPeakEdges_summary.txt regions_TranscriptsToPeakEdges.txt regions_TSSToPeakEdges_summary.txt regions_TSSToPeakEdges.txt"
RnaChipIntegrator.py --project=regions \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_regions.txt
for f in $REGION_OUTPUTS ; do
    diff -f $REF_DATA/ref_$f $f
done
##
#
