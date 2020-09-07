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
    mkdir -p $WORKING_DIR
  else
    rm -rf $WORKING_DIR/data $WORKING_DIR/*.yaml $WORKING_DIR/*.ini
  fi
  cp -r $INPUTS_DIR/* $WORKING_DIR
  rm -rf $WORKDIR_DIR/data/outputs
  echo "Copied contents of $INPUTS_DIR to $WORKING_DIR"
fi
