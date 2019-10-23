from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()

config.General.requestName       = 'WRToTauNu_13TeV'
config.General.workArea          = 'crab'
config.General.transferOutputs   = True
config.General.transferLogs      = True
config.General.instance          = 'prod'

config.Data.outputPrimaryDataset = 'WRToTauNu_13TeV'
config.Data.splitting            = 'EventBased'
config.Data.unitsPerJob          = 1000
config.Data.totalUnits           = 50000
config.Data.outLFNDirBase        = '/store/user/ezhemchu/crab/WRToTauNu_13TeV'
config.Data.publication          = False

config.JobType.pluginName        = 'PrivateMC'
config.JobType.psetName          = 'WRToTauNu_13TeV_GEN-SIM.py'
config.JobType.maxJobRuntimeMin  = 1200
config.JobType.inputFiles        = [
    'WRToTauNu_13TeV_GEN-SIM-RAW.py',
    'WRToTauNu_13TeV_GEN-SIM-RAW-RECO.py',
    'fwjr-merge'
]
config.JobType.scriptExe         = 'WRToTauNu_13TeV.sh'
config.JobType.outputFiles       = [ 
    'WRToTauNu_13TeV_GEN-SIM-RAW.root',
    'WRToTauNu_13TeV_GEN-SIM-RAW-RECO.root'
]

config.JobType.allowUndistributedCMSSW = True

config.Site.storageSite          = 'T3_RU_MEPhI'
