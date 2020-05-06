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
#process.GlobalTag.globaltag='106X_dataRun2_v20'
# MC run2 GlobalTag
#from Configuration.AlCa.autoCond import autoCond
#process.GlobalTag.globaltag=autoCond['run2_mc']
process.GlobalTag.globaltag='106X_mc2017_realistic_v6'


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#'/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/30000/A69A903E-4EF2-CE4D-AD5E-7FCAEE9F8326.root',
	#'/store/data/Run2017B/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/B6EC130B-9D55-8140-AC2E-154EB59E3BBC.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/100000/E059A7E6-BA71-DB42-BF0A-2381E40C21AB.root'
	'/store/mc/RunIISummer19UL17MiniAOD/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/FB807F94-A777-A64B-A9C1-416A53B40464.root'
    )
)

#process.treemaker = cms.EDAnalyzer('TreeMakerMiniAOD')

process.TreeMakerMiniAOD = cms.EDAnalyzer("TreeMakerMiniAOD",
	monitoring = cms.bool(True),
	isMC = cms.bool(True),
	TauSpinnerOn = cms.bool(False),
	looseTauID = cms.bool(True),
	DeepTau = cms.bool(True),
	tauPtMin = cms.double(20),
	piPtMin = cms.double(0),
	tauEtaMax = cms.double(2.3),
	tauDzMax = cms.double(0.2),
        METcut = cms.double(50),
	null = cms.double(-5),
	Triggers = cms.vstring(""),
	tauCollection = cms.string("slimmedTaus"),
	muonCollection = cms.string("slimmedMuons"),
	electronCollection = cms.string("slimmedElectrons"),
	jetCollection = cms.string("slimmedJets"),
	PuppijetCollection = cms.string("slimmedJetsPuppi"),
	metCollection = cms.string("slimmedMETs"),
	vertexCollection = cms.string("offlineSlimmedPrimaryVertices"),
	#genParticleCollection = cms.string("genParticles"),
	#genParticleCollection = cms.string("packedGenParticles"),
	genParticleCollection = cms.string("prunedGenParticles"),
	trackCollection = cms.string("isolatedTracks"), # generalTracks for AOD							
)

#process.load('Tau.TreeMaker.treeMaker_Data-MET_cfi')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('test_MC_miniaod_WJetsToLNu.root')
)

process.p = cms.Path(process.TreeMakerMiniAOD)
