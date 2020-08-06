#!/bin/bash

function usage {
	echo "usage: $scriptname <exepath> <output_dir>"
    echo "  exepath     Path to model executable"
    echo "  output_dir  Directory where output files are stored"
    exit 1
}

EXEPATH=$1
OUTPUT_DIR=$2
FLAGS=${@:3}

if [[ $# -le 2 || $# -gt 8 ]]; then
  echo "ERROR: Invalid ($#) number of arguments"
  usage
  exit 1
fi

echo "Initiating model run..."

if [ ! -f $EXEPATH ]; then
  echo "ERROR: Unable to locate executable $EXEPATH"
  exit 1
else
  if [ ! -d $OUTPUT_DIR ]; then
    echo "Warning: Outputs directory $OUTPUT_DIR does not exist; creating it now"
	mkdir -p $OUTPUT_DIR
  else
    rm $OUTPUT_DIR/* 2> /dev/null
  fi
  
  $EXEPATH $FLAGS
  if [ $? -ne 0 ]; then
    echo "ERROR: Running command $EXEPATH $FLAGS failed"
	exit 1
  else
    echo "Model run completed"
  fi
fi
