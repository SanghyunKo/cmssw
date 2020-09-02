#!/bin/sh

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /u/user/sako/ModHEEP/CMSSW_8_0_32/src/ZprimeTo4l
cmsenv

root -q -b 'analysis.cpp("'$1'")'
