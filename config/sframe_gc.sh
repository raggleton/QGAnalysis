#!/bin/bash

# set -euo pipefail
START=$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_STORED:$LD_LIBRARY_PATH
export PATH=$PATH_STORED:$PATH
# export CMSSW_BASE=$CMSSW_BASE
cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`
cd $START
cd $SFRAME_DIR
SFRAME_DIR=""
source setup.sh
cd $START
printenv | sort
which sframe_main

cd $(dirname $FILENAME)
XML=$(basename $FILENAME)
# when running locally, no need to do any setup, all env vars are there
nice -n 10 sframe_main "$XML"
