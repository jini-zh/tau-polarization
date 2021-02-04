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
# 31Mar2018 GlobalTag for 2017 data
#process.GlobalTag.globaltag='94X_dataRun2_v10'
process.GlobalTag.globaltag='94X_mc2017_realistic_v14'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'/store/data/Run2017D/Tau/MINIAOD/09Aug2019_UL2017-v1/240000/54FFF36B-09F5-0844-9970-9FD18F028F6F.root'
      #'/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/30000/A69A903E-4EF2-CE4D-AD5E-7FCAEE9F8326.root'
      #'/store/data/Run2017F/Tau/MINIAOD/31Mar2018-v1/80000/FAB1EC5A-7737-E811-AA41-008CFAC94018.root'
      #'/store/data/Run2016D/Tau/MINIAOD/17Jul2018-v1/00000/4C16E40F-898B-E811-BD35-A4BF0112BE0A.root'
      #'/store/mc/RunIIFall17MiniAODv2/TTToSemiLeptonic_TuneCP5up_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/9814626D-8843-E811-A239-0CC47A1DF81C.root'
	#'/store/mc/RunIIFall17MiniAODv2/WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/60000/7A96A80F-E9C2-E811-AB1C-0CC47A7C357A.root',
	'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/100000/186120D6-C6A7-E811-9271-FA163E809085.root'
    )
)

#process.treemaker = cms.EDAnalyzer('TreeMakerMiniAOD')

process.TreeMakerMiniAOD = cms.EDAnalyzer("TreeMakerMiniAOD",
	monitoring = cms.bool(False),
	isMC = cms.bool(True),
	TauSpinnerOn = cms.bool(False),
	looseTauID = cms.bool(True),
	DeepTau = cms.bool(False),
	useHLT = cms.bool(False),
	tauPtMin = cms.double(20),
	piPtMin = cms.double(0),
	tauEtaMax = cms.double(2.3),
	tauDzMax = cms.double(0.2),
        METcut = cms.double(0),
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
	genParticleCollection = cms.string("prunedGenParticles"),
	trackCollection = cms.string("isolatedTracks"), # generalTracks for AOD							
)

#process.load('Tau.TreeMaker.treeMaker_Data-MET_cfi')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('test_MC94X_miniaod.root')
)

process.p = cms.Path(process.TreeMakerMiniAOD)
