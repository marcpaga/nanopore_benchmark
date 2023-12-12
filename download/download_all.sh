#!/bin/bash

# Script to download all the benchmarking Nanopore data
# usage:
# bash download_all.sh OUTPUT_DIR

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# make the output dir if does not exists
if [ ! -d $1 ]
then
     mkdir $1
else
     echo "Download directory exists"
fi

source "${SCRIPT_DIR}/download_jain.sh" "$1/jain"
source "${SCRIPT_DIR}/download_verm.sh" "$1/verm"
source "${SCRIPT_DIR}/download_wick.sh" "$1/wick"

