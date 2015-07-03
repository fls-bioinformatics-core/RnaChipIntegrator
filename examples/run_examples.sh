#!/bin/sh -e
#
# Run examples for RnaChipIntegrator
#
TEST_DIR=$(pwd)
if [ ! -d test-output ] ; then
    mkdir test-output
fi
cd test-output
#
# Summits
RnaChipIntegrator.py --project=summits \
    --window=130000 \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_summits.txt

# Regions
RnaChipIntegrator.py --project=regions \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_regions.txt
##
#
