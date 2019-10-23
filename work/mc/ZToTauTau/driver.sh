#!/bin/bash

cmsrun() {
  echo "cmsRun $1"
  cmsRun "$1" 2>&1 | tee ${1%.*}.log
  return ${PIPESTATUS[0]}
}

set -xe

cmsDriver.py ZTT_Tauola_All_hadronic_13TeV_TuneCUETP8M1_cfi.py \
             --fileout file:ZToTauTau_13TeV_GEN-SIM.root \
             -s GEN,SIM \
             --mc --datatier GEN-SIM \
             --conditions auto:phase1_2017_realistic \
             --eventcontent RAWSIM \
             --era Run2_2017 \
             --python_filename ZToTauTau_13TeV_GEN-SIM.py \
             -n 1000 \
             --no_exec

cmsDriver.py -s DIGI,L1,DIGI2RAW,HLT \
             --datatier GEN-SIM-RAW \
             --conditions auto:phase1_2017_realistic \
             --eventcontent RAWSIM \
             --era Run2_2017 \
             --filein file:ZToTauTau_13TeV_GEN-SIM.root \
             --fileout file:ZToTauTau_13TeV_GEN-SIM-RAW.root \
             --python_filename ZToTauTau_13TeV_GEN-SIM-RAW.py \
             -n -1 \
             --mc \
             --no_exec

cmsDriver.py -s RAW2DIGI,L1Reco,RECO \
             --datatier RECO \
             --conditions auto:phase1_2017_realistic \
             --eventcontent AODSIM \
             --era Run2_2017 \
             --filein file:ZToTauTau_13TeV_GEN-SIM-RAW.root \
             --fileout file:ZToTauTau_13TeV_GEN-SIM-RAW-RECO.root \
             --python_filename ZToTauTau_13TeV_GEN-SIM-RAW-RECO.py \
             -n -1 \
             --mc \
             --no_exec

set +x

# cmsrun ZToTauTau_13TeV_GEN-SIM.py 
# cmsrun ZToTauTau_13TeV_GEN-SIM-RAW.py
# cmsrun ZToTauTau_13TeV_GEN-SIM-RAW-RECO.py
