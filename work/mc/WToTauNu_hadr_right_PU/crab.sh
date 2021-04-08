#!/bin/bash

set -xe

cmsRun -j WToTauNu_hadr_right_PU_GEN-SIM-RAW.xml PSet.py
cmsRun -j WToTauNu_hadr_right_PU_GEN-SIM-RAW-RECO.{xml,py}
cmsRun -j WToTauNu_hadr_right_PU_MINIAODSIM.{xml,py}

./fwjr-merge WToTauNu_hadr_right_PU_{GEN-SIM-RAW{,-RECO},MINIAODSIM}.xml > FrameworkJobReport.xml || true
