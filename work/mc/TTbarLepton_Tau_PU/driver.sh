set -xe

cmsDriver.py TTbarLepton_Tau_cff.py \
             --fileout file:TTbarLepton_Tau_PU_GEN-SIM-RAW.root \
             -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT \
             --mc \
             --conditions auto:phase1_2017_realistic \
             --eventcontent RAWSIM \
             --era Run2_2017 \
             --python_filename TTbarLepton_Tau_PU_GEN-SIM-RAW.py \
             -n 10 \
             --no_exec \
             --pileup AVE_35_BX_25ns \
             --pileup_input das:/RelValMinBias_13/CMSSW_9_4_14_UL-94X_mc2017_realistic_v17-v1/GEN-SIM-DIGI-RAW

cmsDriver.py -s RAW2DIGI,L1Reco,RECO \
             --conditions auto:phase1_2017_realistic \
             --eventcontent AODSIM \
             --era Run2_2017 \
             --filein file:TTbarLepton_Tau_PU_GEN-SIM-RAW.root \
             --fileout file:TTbarLepton_Tau_PU_GEN-SIM-RAW-RECO.root \
             --python_filename TTbarLepton_Tau_PU_GEN-SIM-RAW-RECO.py \
             -n -1 \
             --mc \
             --no_exec

cmsDriver.py miniAOD-prod \
             -s PAT \
             --eventcontent MINIAODSIM \
             --mc \
             --conditions auto:phase1_2017_realistic \
             --era Run2_2017 \
             --filein file:TTbarLepton_Tau_PU_GEN-SIM-RAW-RECO.root \
             --fileout TTbarLepton_Tau_PU_MINIAODSIM.root \
             --python_filename TTbarLepton_Tau_PU_MINIAODSIM.py \
             --runUnscheduled \
             -n -1 \
             --no_exec
