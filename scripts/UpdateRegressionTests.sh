#!/bin/bash

scriptname=$0

function usage {
    echo "usage: $scriptname <first> <last>"
    echo "  first      Number of the first regression test to update"
    echo "  expected    Number of the last regression test to update"
    exit 1
}

if [[ $# -ne 2 ]]; then
    echo "ERROR: Invalid ($#) number of arguments"
    usage
    exit 1
fi

FIRST=$1
LAST=$2
if [[ ! $FIRST -le $LAST ]]; then
    echo "ERROR: $FIRST should be less than or equal to $LAST"
    usage
    exit 1
fi

WORKING_DATA_DIR=./data
WORKING_OUTPUTS_DIR=./outputs
REGRESSION_DIR=./test/regression
EXEPATH=./build/bin/Covid19EERAModel
SETUP_SCRIPT=./scripts/SetUpModelRun.sh
RUN_SCRIPT=./scripts/RunModel.sh

for i in $(seq $FIRST $LAST) ; do   
    echo ""
    echo "***************** Updating regression test #$i *****************"

    # Tests 1-6 and 13-18 are against the original model
    if ([[ $i -ge 1 ]] && [[ $i -le 6 ]]) || ([[ $i -ge 13 ]] && [[ $i -le 18 ]]); then
        STRUCT_FLAG=$( echo "-s original" )
    # Tests 7-12 and 18-24 are against the irish model
    else
        STRUCT_FLAG=$( echo "-s irish" )
    fi

    # Tests 1-12 are inference mode
    if [[ $i -ge 1 ]] && [[ $i -le 12 ]]; then
        MODE_FLAG=$( echo "-m inference" )
    # Tests 13-24 are prediction mode
    else 
        MODE_FLAG=$( echo "-m prediction" )
    fi

    FLAGS=$( echo "$STRUCT_FLAG $MODE_FLAG" )
    regression_test_dir=$REGRESSION_DIR/run$i
    
    if [[ ! -d $regression_test_dir ]]; then
        echo "WARNING: Cannot find regression tests directory $regression_test_dir"
        continue
    fi

    $SETUP_SCRIPT $WORKING_DATA_DIR $regression_test_dir/data
    setup=$?

    $RUN_SCRIPT $EXEPATH $WORKING_OUTPUTS_DIR $FLAGS
    run=$?

    output_dir=$regression_test_dir/outputs
    for filepath in $output_dir/*.txt ; do
        rm $filepath
    done
    
    for filepath in $WORKING_OUTPUTS_DIR/*.txt ; do
        # Remove the timestamp from the file before putting it in
        # the regression test outputs directory
        filename=$(basename -- "$filepath")
        cp $filepath $output_dir/${filename::(-24)}.txt
    done
    
    echo "==========================================================================="

done
