#!/bin/bash

set -xe

cmsRun -j TTbarLepton_Tau_PU_GEN-SIM-RAW.xml PSet.py
cmsRun -j TTbarLepton_Tau_PU_GEN-SIM-RAW-RECO.{xml,py}
cmsRun -j TTbarLepton_Tau_PU_MINIAODSIM.{xml,py}

./fwjr-merge TTbarLepton_Tau_PU_{GEN-SIM-RAW{,-RECO},MINIAODSIM}.xml > FrameworkJobReport.xml || true
