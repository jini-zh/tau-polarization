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
#process.GlobalTag.globaltag='106X_mc2017_realistic_v6'
process.GlobalTag.globaltag='94X_mc2017_realistic_v14'


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#'/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/30000/A69A903E-4EF2-CE4D-AD5E-7FCAEE9F8326.root',
	#'/store/data/Run2017B/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/B6EC130B-9D55-8140-AC2E-154EB59E3BBC.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/100000/E059A7E6-BA71-DB42-BF0A-2381E40C21AB.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/FB807F94-A777-A64B-A9C1-416A53B40464.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/F27778B6-F8BD-2049-A473-5CF83F145DEE.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/FC2D7CE4-5223-3942-8295-C8E7E8ED1E1D.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/WZ_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/ECC01154-AC39-584E-836B-4D5B7946508D.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/B9059574-2ACF-A34B-8BBE-08473A6FB46C.root'
	#'/store/mc/RunIIFall17MiniAODv2/WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/60000/7A96A80F-E9C2-E811-AB1C-0CC47A7C357A.root'
	'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/100000/186120D6-C6A7-E811-9271-FA163E809085.root'	
    )
)

#process.treemaker = cms.EDAnalyzer('TreeMakerMiniAOD')

# Rerun tau ID 

updatedTauName = "slimmedTausNewID" #name of pat::Tau collection with new tau-Ids
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, debug = False,
                    updatedTauName = updatedTauName,
                    toKeep = ["dR0p32017v2"#"newDM2017v2"#"2017v2",#"deepTau2017v2", #deepTau TauIDs
                               ])
tauIdEmbedder.runTauID()

# TreeMakder for miniAOD

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
	#tauCollection = cms.string("slimmedTaus"),
	tauCollection = cms.string("slimmedTausNewID"),
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
  fileName = cms.string('test_MC_94X_new_ID_v1.root')
)

process.p = cms.Path(process.rerunMvaIsolationSequence * getattr(process,updatedTauName) * process.TreeMakerMiniAOD)
