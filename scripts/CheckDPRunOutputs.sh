#!/bin/bash

scriptname=$0
DPCHECK_SCRIPT=$( dirname "$scriptname" )/dpcheck.py

function usage {
    echo "usage: $scriptname <actual dp> <expected>"
    echo "  actual dp   directory containing data pipeline directory used for the model run"
    echo "  expected    directory containing expected outputs of the model run"
    exit 1
}

ACTUALDP_OUTPUTS_DIR=$1
EXPECTED_OUTPUTS_DIR=$2

if [[ $# -ne 2 ]]; then
  echo "ERROR: Invalid ($#) number of arguments"
  usage
  exit 1
fi

echo "Checking model outputs..."

if [[ ! -d $ACTUALDP_OUTPUTS_DIR ]]; then
  echo "ERROR: Cannot find directory $ACTUALDP_OUTPUTS_DIR"
  usage
elif [[ ! -d $EXPECTED_OUTPUTS_DIR ]]; then
  echo "ERROR: Cannot find directory $EXPECTED_OUTPUTS_DIR"
  usage
else
  # Find the data pipeline output file hopefully there should be only 1 present...
  DPH5FILE=$( find $ACTUALDP_OUTPUTS_DIR/data/outputs -name '*.h5' )
  EXPECTED_FILE_COUNT=$(find $EXPECTED_OUTPUTS_DIR -maxdepth 1 -type f -name "*.txt" | wc -l)

  echo "Checking for $EXPECTED_FILE_COUNT outputs in $DPH5FILE"
  difference_count=0
  for expected_filepath in $EXPECTED_OUTPUTS_DIR/*.txt
  do
	filename=$(basename $expected_filepath)
    component=$(
		echo $filename |
		awk -F_ '{ if ($2 == "prediction") { match($3, /([a-z]*).txt/, arr); printf("%s", arr[1]); } else { match($4, /step([0-9]*)/, arr); printf("steps/%s/%s", arr[1], $3); } }'
	)

	echo "Checking component: $component" 
	$DPCHECK_SCRIPT "$DPH5FILE" "$component" "$expected_filepath"
	if [[ $? -ne 0 ]]; then
	  echo "ERROR: Files $expected_filepath and $component are different"
	  difference_count=$((difference_count+1))
	else
      echo "Files $expected_filepath and $component are identical"
	fi
  done
	if [ $difference_count -ne 0 ]; then
	  echo "SUMMARY - FAILURE: $difference_count components different"
	  exit 1
	else
	  echo "SUMMARY - SUCCESS: Model outputs match expected outputs"
	  exit 0
	fi
  fi
fi
  