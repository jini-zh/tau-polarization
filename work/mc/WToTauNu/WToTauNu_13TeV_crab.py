from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.General.requestName       = 'WToTauNu_13TeV'
config.General.workArea          = 'crab'
config.General.transferOutputs   = True
config.General.transferLogs      = True

config.Data.outputPrimaryDataset = 'WToTauNu_13TeV'
config.Data.splitting            = 'EventBased'
config.Data.unitsPerJob          = 1000
config.Data.totalUnits           = 40000
config.Data.outLFNDirBase        = '/store/user/ezhemchu/crab/WToTauNu_13TeV'
config.Data.publication          = False

config.JobType.pluginName        = 'PrivateMC'
config.JobType.psetName          = 'WToTauNu_13TeV_GEN-SIM.py'
config.JobType.maxJobRuntimeMin  = 900
config.JobType.inputFiles        = [
    'WToTauNu_13TeV_GEN-SIM-RAW.py',
    'WToTauNu_13TeV_GEN-SIM-RAW-RECO.py',
    'fwjr-merge'
]
config.JobType.scriptExe         = 'WToTauNu_13TeV.sh'
config.JobType.outputFiles       = [ 
    'WToTauNu_13TeV_GEN-SIM-RAW.root',
    'WToTauNu_13TeV_GEN-SIM-RAW-RECO.root'
]

config.JobType.allowUndistributedCMSSW = True

config.Site.storageSite          = 'T3_RU_MEPhI'
