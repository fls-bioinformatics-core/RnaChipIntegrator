#!/bin/bash
#
# Run examples for RnaChipIntegrator
#
# Assertion functions for tests
function assert_equal {
    # Check two files are the same
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
# Check XLS file is produced by --xls
RnaChipIntegrator --name=test_xls \
    --xls --number=4 \
    $TEST_DIR/ExpressionData.txt \
    $TEST_DIR/ChIP_regions.txt
if [ $? -ne 0 ] || [ ! -f test_xls.xls ] ; then
    echo "XLS test: FAILED"
    exit 1
fi
echo "XLS test: OK"
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
exit 0
##
#
