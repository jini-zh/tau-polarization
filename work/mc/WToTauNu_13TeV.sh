#!/bin/bash

set -xe

cmsRun -j WToTauNu_13TeV_GEN-SIM.xml PSet.py
cmsRun -j WToTauNu_13TeV_GEN-SIM-RAW.xml WToTauNu_13TeV_GEN-SIM-RAW.py
cmsRun -j WToTauNu_13TeV_GEN-SIM-RAW-RECO.xml WToTauNu_13TeV_GEN-SIM-RAW-RECO.py

./fwjr-merge WToTauNu_13TeV_GEN-SIM{,-RAW{,-RECO}}.xml > FrameworkJobReport.xml
