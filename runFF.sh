#!/bin/bash
source utils/setup_python.sh
export LCG_RELEASE=96
source utils/setup_cvmfs_sft.sh
source utils/bashFunctionCollection.sh

logandrun python fake-factor-derivation/FF.py

