#!/bin/bash
#
# Run examples for RnaChipIntegrator
#
# Assertion functions for tests
_PASSED=0
_FAILED=0
function report_tests {
    # Report summary of tests passed and failed
    local n_tests=$((_PASSED+_FAILED))
    echo "---------------------------------------------------------"
    echo "Ran $n_tests tests: $_PASSED passed, $_FAILED failed"
    if [ $_FAILED -ne 0 ] ; then
	return 1
    else
	return 0
    fi
}
function run_test {
    # Run a command and check outputs
    # Takes following arguments:
    # --command CMD: specify the command to execute
    # --expected FILES: list of file names which should
    #                have reference versions to compare
    #                against
    # --must_exist FILES: list of file names which should
    #                exist after the command has run
    # --status INT: exit code to check (if command was
    #                run externally)
    local test_name=$1
    local command=
    local expected_outputs=
    local check_exists=
    local exit_status=
    local working_dir=
    local test_status=
    # Collect arguments
    shift
    while [ $# -gt 0 ] ; do
	case $1 in
	    --command)
		command=$2
		shift
		;;
	    --expected)
		expected_outputs=$2
		shift
		;;
	    --must_exist)
		check_exists=$2
		shift
		;;
	    --status)
		exit_status=$2
		shift
		;;
	    *)
		echo "$test_name: SKIPPED (unrecognised test argument '$1')"
		return
		;;
	esac
	shift
    done
    echo "---------------------------------------------------------"
    echo test_name: $test_name
    echo command: $command
    echo expected_outputs: $expected_outputs
    echo check_exists: $check_exists
    echo exit_status: $exit_status
    echo PWD: $(pwd)
    # If command supplied then run it
    if [ ! -z "$command" ] ; then
	working_dir=$(mktemp -d --tmpdir=$(pwd))
	echo working_dir: $working_dir
	cd $working_dir
	echo "Running command"
	$command 1>STDOUT 2>STDERR
	exit_status=$?
	echo "Exit status $exit_status"
    fi
    # Check exit status
    if [ ! -z "$exit_status" ] ; then
	if [ $exit_status -ne 0 ] ; then
	    echo Failed exit status check
	    test_status=FAILED
	fi
    fi
    # Compare expected outputs
    for f in $expected_outputs ; do
	assert_equal $REF_DATA/ref_$f $f
	if [ $? -ne 0 ] ; then
	    echo Failed output comparison check
	    test_status=FAILED
	fi
    done
    # Check existence
    for f in $check_exists ; do
	if [ ! -e $f ] ; then
	    echo "$f: missing"
	    echo Failed output existence check
	    test_status=FAILED
	fi
    done
    # Set test status if no failures
    if [ -z "$test_status" ] ; then
	test_status=OK
    fi
    echo test_status: $test_status
    # Report logs from failed job
    if [ $test_status == FAILED ] ; then
	for f in STDOUT STDERR ; do
	    if [ -e $f ] ; then
		echo "===== $test_name: $f ====="
		cat $f
	    fi
	done
    fi
    # Clean up any working area
    if [ ! -z "$working_dir" ] ; then
	cd ..
	#rm -rf $working_dir
    fi
    # Test counts
    case $test_status in
	OK)
	    _PASSED=$((_PASSED+1))
	    ;;
	FAILED)
	    _FAILED=$((_FAILED+1))
	    ;;
    esac
    # Finish
    echo "---------------------------------------------------------"
    echo "TEST: $test_name: $test_status"
}
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
    rm -rf test-output/*
fi
cd test-output
#
# Summits
SUMMIT_OUTPUTS="summits_peak_centric.txt summits_gene_centric.txt"
run_test "Input summits" \
	 --expected "$SUMMIT_OUTPUTS" \
	 --command "RnaChipIntegrator --name=summits --cutoff=130000 --number=4 $TEST_DIR/ExpressionData.txt $TEST_DIR/ChIP_summits.txt"
#
# Regions
REGION_OUTPUTS="regions_peak_centric.txt regions_gene_centric.txt"
run_test "Input regions" \
	 --expected "$REGION_OUTPUTS" \
	 --command "RnaChipIntegrator --name=regions --number=4 $TEST_DIR/ExpressionData.txt $TEST_DIR/ChIP_regions.txt"
#
# Check XLSX file is produced by --xlsx
run_test "XLSX output" \
	 --must_exist "test_xls.xlsx" \
	 --command "RnaChipIntegrator --name=test_xls --xlsx --compact $TEST_DIR/ExpressionData.txt $TEST_DIR/ChIP_regions.txt"
#
# Check single line output from --compact
COMPACT_OUTPUTS="compact_peak_centric.txt compact_gene_centric.txt"
run_test "Compact output" \
	 --expected "$COMPACT_OUTPUTS" \
	 --command "RnaChipIntegrator --name=compact --compact $TEST_DIR/ExpressionData.txt $TEST_DIR/ChIP_regions.txt"
#
# Check --cutoff=0 works
ZERO_CUTOFF_OUTPUTS="zero_cutoff_peak_centric.txt zero_cutoff_gene_centric.txt"
run_test "Zero cutoff" \
	 --expected "$ZERO_CUTOFF_OUTPUTS" \
	 --command "RnaChipIntegrator --name=zero_cutoff --cutoff=0 --xlsx $TEST_DIR/ExpressionData.txt $TEST_DIR/ChIP_regions.txt"
#
# Check inclusion of peak IDs
PEAK_ID_OUTPUTS="peak_id_peak_centric.txt peak_id_gene_centric.txt"
run_test "Peak IDs" \
	 --expected "$PEAK_ID_OUTPUTS" \
	 --command "RnaChipIntegrator --name=peak_id --number=4 --peak_id=1 --peak_cols=2,3,4 --number=4 $TEST_DIR/ExpressionData.txt $TEST_DIR/ChIP_regions_with_IDs.txt"
#
# Check using --analyses=peak_centric
PEAK_CENTRIC_OUTPUTS="test_peakcentric_peak_centric.txt"
run_test "Peak-centric-only" \
	 --expected "$PEAK_CENTRIC_OUTPUTS" \
	 --command "RnaChipIntegrator --name=test_peakcentric --number=4 --analyses=peak_centric --number=4 $TEST_DIR/ExpressionData.txt $TEST_DIR/ChIP_regions.txt"
#
# Check using --analyses=gene_centric
GENE_CENTRIC_OUTPUTS="test_genecentric_gene_centric.txt"
run_test "Gene-centric-only" \
	 --expected "$GENE_CENTRIC_OUTPUTS" \
	 --command "RnaChipIntegrator --name=test_genecentric --number=4 --analyses=gene_centric $TEST_DIR/ExpressionData.txt $TEST_DIR/ChIP_regions.txt"
#
# Batch mode with multiple cutoffs
BATCH_MULTIPLE_CUTOFFS="test_batch_multi_cutoff_peak_centric.txt"
run_test "Batch mode: multiple cutoffs" \
	 --expected "$BATCH_MULTIPLE_CUTOFFS" \
	 --command "RnaChipIntegrator-batch --name=test_batch_multi_cutoff --number=4 --cutoffs=50000,100000,150000 $TEST_DIR/ExpressionData.txt $TEST_DIR/ChIP_regions.txt"
#
# Batch mode with multiple peak files (peak-centric only)
run_test "Batch mode: multiple peak files" \
	 --expected "$BATCH_MULTIPLE_PEAKS" \
	 --command "RnaChipIntegrator-batch --name=test_batch_multi_peaks --number=4 $TEST_DIR/ExpressionData.txt $TEST_DIR/peaks1.txt $TEST_DIR/peaks2.txt"
#
# Finished
report_tests
exit $?
##
#
