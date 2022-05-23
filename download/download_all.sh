#!/bin/bash

# Script to download all the benchmarking Nanopore data
# usage:
# bash download_all.sh OUTPUT_DIR

# make the output dir if does not exists
if [ ! -d $1 ]
then
     mkdir $1
else
     echo "Download directory exists"
fi

source download_jain.sh "$1/jain"
source download_verm.sh "$1/verm"
source download_wick.sh "$1/wick"