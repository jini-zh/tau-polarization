from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.General.requestName       = 'ZToTauTau_13TeV'
config.General.workArea          = 'crab'
config.General.transferOutputs   = True
config.General.transferLogs      = True

config.Data.outputPrimaryDataset = 'ZToTauTau_13TeV'
config.Data.splitting            = 'EventBased'
config.Data.unitsPerJob          = 1000
config.Data.totalUnits           = 50000
config.Data.outLFNDirBase        = '/store/user/ezhemchu/crab/ZToTauTau_13TeV'
config.Data.publication          = False

config.JobType.pluginName        = 'PrivateMC'
config.JobType.psetName          = 'ZToTauTau_13TeV_GEN-SIM.py'
config.JobType.maxJobRuntimeMin  = 900
config.JobType.inputFiles        = [
    'ZToTauTau_13TeV_GEN-SIM-RAW.py',
    'ZToTauTau_13TeV_GEN-SIM-RAW-RECO.py',
    'fwjr-merge'
]
config.JobType.scriptExe         = 'ZToTauTau_13TeV.sh'
config.JobType.outputFiles       = [ 
    'ZToTauTau_13TeV_GEN-SIM-RAW.root',
    'ZToTauTau_13TeV_GEN-SIM-RAW-RECO.root'
]

config.JobType.allowUndistributedCMSSW = True

config.Site.storageSite          = 'T3_RU_MEPhI'
