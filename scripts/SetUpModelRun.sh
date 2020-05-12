#!/bin/bash

function usage {
	echo "usage: $scriptname <working_dir> <inputs_dir>"
    echo "  working_dir     working directory where the model executable expects to find inputs"
    echo "  inputs_dir      directory where the inputs for the regression test are stored"
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

if [[ ! -d $WORKING_DIR ]]; then
  echo "ERROR: Unable to locate working data directory $WORKING_DIR"
  usage
elif [[ ! -d $INPUTS_DIR ]]; then
  echo "ERROR: Unable to locate input data directory $INPUTS_DIR"
  usage
else
  rm -rf $WORKING_DIR/*.csv $WORKING_DIR/*.ini
  cp $INPUTS_DIR/* $WORKING_DIR
  echo "Copied contents of $INPUTS_DIR to $WORKING_DIR"
fi
