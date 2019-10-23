#!/bin/bash

set -xe

cmsRun -j ZToTauTau_13TeV_GEN-SIM.xml PSet.py
cmsRun -j ZToTauTau_13TeV_GEN-SIM-RAW.xml ZToTauTau_13TeV_GEN-SIM-RAW.py
cmsRun -j ZToTauTau_13TeV_GEN-SIM-RAW-RECO.xml ZToTauTau_13TeV_GEN-SIM-RAW-RECO.py

./fwjr-merge ZToTauTau_13TeV_GEN-SIM{,-RAW{,-RECO}}.xml > FrameworkJobReport.xml
