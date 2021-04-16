#!/bin/bash

set -xe

cmsRun -j TTbarLepton_TauRight_GEN-SIM-RAW.xml PSet.py
cmsRun -j TTbarLepton_TauRight_GEN-SIM-RAW-RECO.{xml,py}
cmsRun -j TTbarLepton_TauRight_MINIAODSIM.{xml,py}

./fwjr-merge TTbarLepton_TauRight_{GEN-SIM-RAW{,-RECO},MINIAODSIM}.xml > FrameworkJobReport.xml || true
