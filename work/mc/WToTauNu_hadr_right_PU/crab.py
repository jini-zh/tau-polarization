from CRABClient.UserUtilities import config

config = config()

config.General.requestName     = 'WToTauNu_hadr_right_PU'
config.General.workArea        = 'crab'
config.General.transferOutputs = True
config.General.transferLogs    = False

config.Data.outputPrimaryDataset = 'WToTauNu_hadr_right_PU'
config.Data.splitting            = 'EventBased'
config.Data.unitsPerJob          = 100
config.Data.totalUnits           = 100000
config.Data.outLFNDirBase        = '/store/user/ezhemchu/crab/WToTauNu_hadr_right_PU'
config.Data.publication          = False

config.JobType.pluginName        = 'PrivateMC'
config.JobType.psetName          = 'WToTauNu_hadr_right_PU_GEN-SIM-RAW.py'
config.JobType.inputFiles        = [
    'WToTauNu_hadr_right_PU_GEN-SIM-RAW.py',
    'WToTauNu_hadr_right_PU_GEN-SIM-RAW-RECO.py',
    'WToTauNu_hadr_right_PU_MINIAODSIM.py',
    'fwjr-merge'
]
config.JobType.scriptExe         = 'crab.sh'
config.JobType.outputFiles       = [
  'WToTauNu_hadr_right_PU_MINIAODSIM.root',
]
config.JobType.allowUndistributedCMSSW = True
config.JobType.disableAutomaticOutputCollection = True
config.JobType.maxMemoryMB       = 3000

config.Site.storageSite          = 'T2_RU_JINR'
