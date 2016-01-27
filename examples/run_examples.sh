#!/bin/sh
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
SUMMIT_OUTPUTS="summits_peak_centric.txt summits_gene_centric.txt"
RnaChipIntegrator --name=summits \
    --cutoff=130000 \
    --number=4 \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_summits.txt
for f in $SUMMIT_OUTPUTS ; do
    diff -q $REF_DATA/ref_$f $f
    if [ $? -ne 0 ] ; then
	echo "$f: doesn't match reference data:"
	diff $REF_DATA/ref_$f $f
	exit 1
    fi
done
# Regions
REGION_OUTPUTS="regions_peak_centric.txt regions_gene_centric.txt"
RnaChipIntegrator --name=regions \
    --number=4 \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_regions.txt
for f in $REGION_OUTPUTS ; do
    diff -q $REF_DATA/ref_$f $f
    if [ $? -ne 0 ] ; then
	echo "$f: doesn't match reference data:"
	diff $REF_DATA/ref_$f $f
	exit 1
    fi
done
exit 0
##
#
