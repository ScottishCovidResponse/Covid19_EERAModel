#!/bin/bash

function usage {
	echo "usage: $scriptname <working_dir> <inputs_dir>"
    echo "  working_dir     working directory where the model executable expects to find the datapipeline directory"
    echo "  inputs_dir      directory where the test data pipeline is stored"
    exit 1
}

WORKING_DIR=$1
INPUTS_DIR=$2

if [[ $# -ne 2 ]]; then
  echo "ERROR: Invalid ($#) number of arguments"
  usage
  exit 1
fi

echo "Setting up model run..."

if [[ ! -d $INPUTS_DIR ]]; then
  echo "ERROR: Unable to locate test data pipeline directory $INPUTS_DIR"
  usage
else
  if [[ ! -d $WORKING_DIR ]]; then
    # No datapipeline directory - create and download the data
    mkdir -p $WORKING_DIR
    cp $INPUTS_DIR/config.yaml $INPUTS_DIR/parameters.ini $WORKING_DIR
    python3 -m data_pipeline_api.registry.download --config $WORKING_DIR/config.yaml
  else
    # Directory exists - assume still valid and just delete results
    rm -rf $WORKING_DIR/data/outputs $WORKING_DIR/access*.yaml
  fi

  echo "Configured data pipeline contents of $WORKING_DIR from $INPUTS_DIR"
fi
