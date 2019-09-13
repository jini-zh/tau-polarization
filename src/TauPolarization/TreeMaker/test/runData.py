import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMaker")


process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.GlobalTag.globaltag='80X_dataRun2_Prompt_v14'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	'root://cms-xrd-global.cern.ch//store/data/Run2016H/Tau/AOD/PromptReco-v2/000/281/207/00000/CA27ACAD-6282-E611-BF27-FA163ED9B340.root'
      #'root://cms-xrd-global.cern.ch//store/data/Run2016D/Tau/AOD/23Sep2016-v1/100000/000A4C44-C391-E611-B70D-001E67E6F4CC.root'
  #  'root://cms-xrd-global.cern.ch//store/data/Run2016B/JetHT/AOD/22Feb2017_ver2-v1/00000/00335C56-32FA-E611-B230-0025905B85F6.root'
	#'root://cms-xrd-global.cern.ch//store/data/Run2016G/Tau/AOD/PromptReco-v1/000/278/820/00000/004A9416-2364-E611-9E82-02163E011C51.root'
	#'root://srmcms.pic.es//store/data/Run2016G/Tau/AOD/PromptReco-v1/000/278/820/00000/004A9416-2364-E611-9E82-02163E011C51.root'
	#'root://cms-xrd-global.cern.ch//store/data/Run2016D/Tau/AOD/PromptReco-v2/000/276/315/00000/1448B5D3-FC44-E611-85BF-02163E014522.root'
    )
)

process.treemaker = cms.EDAnalyzer('TreeMaker'
)

process.load('Tau.TreeMaker.treeMaker_Data-MET_cfi')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('test_data_aod.root')
)

process.p = cms.Path(process.treemaker)
