#!/bin/bash

scriptname=$0

function usage {
    echo "usage: $scriptname <actual> <expected>"
    echo "  actual      directory containing actual outputs of the model run"
    echo "  expected    directory containing expected outputs of the model run"
    exit 1
}

ACTUAL_OUTPUTS_DIR=$1
EXPECTED_OUTPUTS_DIR=$2

if [[ $# -ne 2 ]]; then
  echo "ERROR: Invalid ($#) number of arguments"
  usage
  exit 1
fi

echo "Checking model outputs..."

if [[ ! -d $ACTUAL_OUTPUTS_DIR ]]; then
  echo "ERROR: Cannot find directory $ACTUAL_OUTPUTS_DIR"
  usage
elif [[ ! -d $EXPECTED_OUTPUTS_DIR ]]; then
  echo "ERROR: Cannot find directory $EXPECTED_OUTPUTS_DIR"
  usage
else
  ACTUAL_FILE_COUNT=$(find $ACTUAL_OUTPUTS_DIR -maxdepth 1 -type f -name "*.txt" | wc -l)
  EXPECTED_FILE_COUNT=$(find $EXPECTED_OUTPUTS_DIR -maxdepth 1 -type f -name "*.txt" | wc -l)
  if [[ $ACTUAL_FILE_COUNT -ne $EXPECTED_FILE_COUNT ]]; then
    echo "ERROR: Found $ACTUAL_FILE_COUNT files in $ACTUAL_OUTPUTS_DIR; expected $EXPECTED_FILE_COUNT"
	exit 1
  else
    echo "Found $ACTUAL_FILE_COUNT files in directory $ACTUAL_OUTPUTS_DIR as expected"
	missing_file_count=0
	difference_count=0
	for expected_filepath in $EXPECTED_OUTPUTS_DIR/*.txt
	do
	  actual_filepath=$ACTUAL_OUTPUTS_DIR/$(basename $expected_filepath)
	  if [[ ! -f $actual_filepath ]]; then
	    echo "ERROR: Failed to find expected file $actual_filepath"
		missing_file_count=$((missing_file_count+1))
	  else
	    cmp --silent $expected_filepath $actual_filepath
		if [[ $? -ne 0 ]]; then
		  echo "ERROR: Files $expected_filepath and $actual_filepath are different"
		  difference_count=$((difference_count+1))
		else
          echo "Files $expected_filepath and $actual_filepath are identical"
		fi
	  fi
	done
	if [ $missing_file_count -ne 0 ] || [ $difference_count -ne 0 ]; then
	  echo "SUMMARY - FAILURE: $missing_file_count files missing; $difference_count files different"
	  exit 1
	else
	  echo "SUMMARY - SUCCESS: Model outputs match expected outputs"
	  exit 0
	fi
  fi
fi
  