#!/bin/bash
#
# Run examples for RnaChipIntegrator
#
# Assertion functions for tests
function assert_equal {
    # Check two files are the same
    if [ ! -e $1 ] ; then
	echo "$1: missing reference data"
	return 1
    elif [ ! -e $2 ] ; then
	echo "$2: missing"
	return 1
    fi
    diff -q $1 $2
    if [ $? -ne 0 ] ; then
	echo "$2: doesn't match reference data:"
	diff $1 $2
	return 1
    else
	return 0
    fi
}
#
# Initialise and set up dir for test outputs
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
    assert_equal $REF_DATA/ref_$f $f
    if [ $? -ne 0 ] ; then
	echo "'Summit' test: FAILED"
	exit 1
    fi
done
echo "'Summit' test: OK"
# Regions
REGION_OUTPUTS="regions_peak_centric.txt regions_gene_centric.txt"
RnaChipIntegrator --name=regions \
    --number=4 \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_regions.txt
for f in $REGION_OUTPUTS ; do
    assert_equal $REF_DATA/ref_$f $f
    if [ $? -ne 0 ] ; then
	echo "'Region' test: FAILED"
	exit 1
    fi
done
echo "'Region' test: OK"
# Check XLSX file is produced by --xlsx
RnaChipIntegrator --name=test_xls \
    --xlsx --compact \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_regions.txt
if [ $? -ne 0 ] || [ ! -f test_xls.xlsx ] ; then
    echo "XLS test: FAILED"
    exit 1
fi
echo "XLSX test: OK"
# Check single line output from --compact
RnaChipIntegrator --name=compact \
    --compact \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_regions.txt
COMPACT_OUTPUTS="compact_peak_centric.txt compact_gene_centric.txt"
for f in $COMPACT_OUTPUTS ; do
    assert_equal $REF_DATA/ref_$f $f
    if [ $? -ne 0 ] ; then
	echo "'Compact' test: FAILED"
	exit 1
    fi
done
echo "'Compact' test: OK"
# Check --cutoff=0 works
RnaChipIntegrator --name=zero_cutoff \
    --cutoff=0 \
    --xlsx \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_regions.txt
if [ $? -ne 0 ] ; then
    echo "'Zero cutoff' test: FAILED"
    exit 1
fi
ZERO_CUTOFF_OUTPUTS="zero_cutoff_peak_centric.txt zero_cutoff_gene_centric.txt"
for f in $ZERO_CUTOFF_OUTPUTS ; do
    assert_equal $REF_DATA/ref_$f $f
    if [ $? -ne 0 ] ; then
	echo "'Zero cutoff' test: FAILED"
	exit 1
    fi
done
echo "'Zero cutoff' test: OK"
# Check using --analyses=peak_centric
RnaChipIntegrator --name=test_peakcentric \
    --analyses=peak_centric \
    --number=4 \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_regions.txt
PEAK_CENTRIC_OUTPUTS="test_peakcentric_peak_centric.txt"
for f in $PEAK_CENTRIC_OUTPUTS ; do
    assert_equal $REF_DATA/ref_$f $f
    if [ $? -ne 0 ] ; then
	echo "Peak-centric-only test: FAILED"
	exit 1
    fi
done
echo "Peak-centric-only test: OK"
# Check using --analyses=gene_centric
RnaChipIntegrator --name=test_genecentric \
    --analyses=gene_centric \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_regions.txt
GENE_CENTRIC_OUTPUTS="test_genecentric_gene_centric.txt"
for f in $PEAK_CENTRIC_OUTPUTS ; do
    assert_equal $REF_DATA/ref_$f $f
    if [ $? -ne 0 ] ; then
	echo "Gene-centric-only test: FAILED"
	exit 1
    fi
done
echo "Gene-centric-only test: OK"
exit 0
##
#
