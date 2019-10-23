#!/bin/bash

set -xe

cmsRun -j WRToTauNu_13TeV_GEN-SIM.xml PSet.py
cmsRun -j WRToTauNu_13TeV_GEN-SIM-RAW.xml WRToTauNu_13TeV_GEN-SIM-RAW.py
cmsRun -j WRToTauNu_13TeV_GEN-SIM-RAW-RECO.xml WRToTauNu_13TeV_GEN-SIM-RAW-RECO.py

./fwjr-merge WRToTauNu_13TeV_GEN-SIM{,-RAW{,-RECO}}.xml > FrameworkJobReport.xml
