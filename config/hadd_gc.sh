#!/bin/bash

# set -euo pipefail
START=$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_STORED:$LD_LIBRARY_PATH
export PATH=$PATH_STORED:$PATH
# export CMSSW_BASE=$CMSSW_BASE
cd $CMSSW_BASE/src
eval `scramv1 runtime -sh`
cd $STARTDIR
pwd

$CMD
