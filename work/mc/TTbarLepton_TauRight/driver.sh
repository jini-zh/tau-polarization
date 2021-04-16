set -xe

cmsDriver.py TTbarLepton_TauRight_cff.py \
             --fileout file:TTbarLepton_TauRight_GEN-SIM-RAW.root \
             -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT \
             --mc \
             --conditions auto:phase1_2017_realistic \
             --eventcontent RAWSIM \
             --era Run2_2017 \
             --python_filename TTbarLepton_TauRight_GEN-SIM-RAW.py \
             -n 10 \
             --no_exec

cmsDriver.py -s RAW2DIGI,L1Reco,RECO \
             --conditions auto:phase1_2017_realistic \
             --eventcontent AODSIM \
             --era Run2_2017 \
             --filein file:TTbarLepton_TauRight_GEN-SIM-RAW.root \
             --fileout file:TTbarLepton_TauRight_GEN-SIM-RAW-RECO.root \
             --python_filename TTbarLepton_TauRight_GEN-SIM-RAW-RECO.py \
             -n -1 \
             --mc \
             --no_exec

cmsDriver.py miniAOD-prod \
             -s PAT \
             --eventcontent MINIAODSIM \
             --mc \
             --conditions auto:phase1_2017_realistic \
             --era Run2_2017 \
             --filein file:TTbarLepton_TauRight_GEN-SIM-RAW-RECO.root \
             --fileout TTbarLepton_TauRight_MINIAODSIM.root \
             --python_filename TTbarLepton_TauRight_MINIAODSIM.py \
             --runUnscheduled \
             -n -1 \
             --no_exec
