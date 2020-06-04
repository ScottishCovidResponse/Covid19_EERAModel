#!/bin/bash

scriptname=$0

WORKING_DATA_DIR=./data
WORKING_OUTPUTS_DIR=./outputs
REGRESSION_DIR=./test/regression
EXEPATH=./build/bin/Covid19EERAModel
SETUP_SCRIPT=./scripts/SetUpModelRun.sh
RUN_SCRIPT=./scripts/RunModel.sh
CHECK_SCRIPT=./scripts/CheckRunOutputs.sh

function usage {
    echo "usage: $scriptname"
    exit 1
}

if [[ $# -ne 0 ]]; then
  echo "ERROR: Invalid ($#) number of arguments"
  usage
  exit 1
fi

failures=0
successes=0

echo "Running regression tests..."
echo ""
echo "==========================================================================="

for i in $(seq 1 6) ; do   
  echo ""
  echo "***************** Running regression test #$i *****************"
  
  regression_test_dir=$REGRESSION_DIR/run$i
  
  $SETUP_SCRIPT $WORKING_DATA_DIR $regression_test_dir/data
  setup=$?

  $RUN_SCRIPT $EXEPATH $WORKING_OUTPUTS_DIR
  run=$?

  $CHECK_SCRIPT $WORKING_OUTPUTS_DIR $regression_test_dir/outputs
  check=$?
  
  if [[ setup -ne 0 ]] || [[ run -ne 0 ]] || [[ check -ne 0 ]]; then
    failures=$((failures+1))
  else
    successes=$((successes+1))
  fi
  
  echo "==========================================================================="

done

echo "REGRESSION TESTS COMPLETE: $successes PASSED; $failures FAILED"

if [[ $failures -ne 0 ]]; then
  exit 1
fi
