from CRABClient.UserUtilities import config

config = config()

config.General.requestName     = 'TTbarLepton_TauRight'
config.General.workArea        = 'crab'
config.General.transferOutputs = True
config.General.transferLogs    = False

config.Data.outputPrimaryDataset = 'TTbarLepton_TauRight'
config.Data.splitting            = 'EventBased'
config.Data.unitsPerJob          = 1000
config.Data.totalUnits           = 100000
config.Data.outLFNDirBase        = '/store/user/ezhemchu/crab/TTbarLepton_TauRight'
config.Data.publication          = False

config.JobType.pluginName        = 'PrivateMC'
config.JobType.psetName          = 'TTbarLepton_TauRight_GEN-SIM-RAW.py'
config.JobType.inputFiles        = [
    'TTbarLepton_TauRight_GEN-SIM-RAW.py',
    'TTbarLepton_TauRight_GEN-SIM-RAW-RECO.py',
    'TTbarLepton_TauRight_MINIAODSIM.py',
    'fwjr-merge'
]
config.JobType.scriptExe         = 'crab.sh'
config.JobType.outputFiles       = [
  'TTbarLepton_TauRight_MINIAODSIM.root',
]
config.JobType.allowUndistributedCMSSW = True
config.JobType.disableAutomaticOutputCollection = True
config.JobType.maxMemoryMB       = 3000

config.Site.storageSite          = 'T2_RU_JINR'
