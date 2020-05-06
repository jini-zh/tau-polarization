import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMakerMiniAOD")

process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

# UL2017 GlobalTag
#process.GlobalTag.globaltag='106X_dataRun2_v20'
# 31Mar2018 GlobalTag for 2017 Tau data
#process.GlobalTag.globaltag='94X_dataRun2_ReReco_EOY17_v6'
# JetHT UL2017 Data GlobalTag
process.GlobalTag.globaltag='106X_dataRun2_v20'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'/store/data/Run2017D/Tau/MINIAOD/09Aug2019_UL2017-v1/240000/54FFF36B-09F5-0844-9970-9FD18F028F6F.root'
     # '/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/30000/A69A903E-4EF2-CE4D-AD5E-7FCAEE9F8326.root',
	#'/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/60000/C88B9135-9EF6-8443-83B9-618381FAB602.root',
	#'/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/60000/C7F2018B-B41E-4F44-BB1E-18F84BC5E825.root',
	#'/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/60000/CA47CC2E-702D-5647-9A8C-C3796BAF29E6.root',
	#'/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/60000/D0C81207-F01D-A044-99C3-24BB77219FF5.root',
      #'/store/data/Run2017F/Tau/MINIAOD/31Mar2018-v1/80000/FAB1EC5A-7737-E811-AA41-008CFAC94018.root'
	'/store/data/Run2017B/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/B6EC130B-9D55-8140-AC2E-154EB59E3BBC.root'
    )
)

#process.treemaker = cms.EDAnalyzer('TreeMakerMiniAOD')

process.TreeMakerMiniAOD = cms.EDAnalyzer("TreeMakerMiniAOD",
	monitoring = cms.bool(False),
	isMC = cms.bool(False),
	looseTauID = cms.bool(True),
	DeepTau = cms.bool(True),
	tauPtMin = cms.double(20),
	piPtMin = cms.double(0),
	tauEtaMax = cms.double(2.3),
	tauDzMax = cms.double(0.2),
        METcut = cms.double(40),
	null = cms.double(-5),
	Triggers = cms.vstring(""),
	tauCollection = cms.string("slimmedTaus"),
	muonCollection = cms.string("slimmedMuons"),
	electronCollection = cms.string("slimmedElectrons"),
	jetCollection = cms.string("slimmedJets"),
	PuppijetCollection = cms.string("slimmedJetsPuppi"),
	metCollection = cms.string("slimmedMETs"),
	vertexCollection = cms.string("offlineSlimmedPrimaryVertices"),
	genParticleCollection = cms.string("genParticles"),
	trackCollection = cms.string("isolatedTracks"), # generalTracks for AOD							
)

#process.load('Tau.TreeMaker.treeMaker_Data-MET_cfi')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('test_data_miniaod_JetHT.root')
)

process.p = cms.Path(process.TreeMakerMiniAOD)
