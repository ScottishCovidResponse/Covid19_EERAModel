#!/bin/bash

PIPELINE_DIR=./test/datapipeline
WORKING_PIPELINE_DIR=./datapipeline
SETUP_DP_SCRIPT=./scripts/SetUpDPModelRun.sh
UNIT_TESTS=./build/bin/Covid19EERAModel-unit_tests

# Setup the data pipeline directory
echo "Setting up data pipeline data..."
$SETUP_DP_SCRIPT $WORKING_PIPELINE_DIR $PIPELINE_DIR

if [[ $? -ne 0 ]]; then
   echo "Error: Problem downloading data pipeline files"
   exit 1
else
   echo "Running unit tests..."
   $UNIT_TESTS
   exit $?
fi