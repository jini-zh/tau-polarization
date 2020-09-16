import FWCore.ParameterSet.Config as cms

process = cms.Process("TTbarTauLepton")

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


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#'/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/30000/A69A903E-4EF2-CE4D-AD5E-7FCAEE9F8326.root',
	#'/store/data/Run2017B/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/B6EC130B-9D55-8140-AC2E-154EB59E3BBC.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/100000/E059A7E6-BA71-DB42-BF0A-2381E40C21AB.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/FB807F94-A777-A64B-A9C1-416A53B40464.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/F27778B6-F8BD-2049-A473-5CF83F145DEE.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/F7134BB1-4D96-EB48-9D94-665CA6A02505.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/270000/84EDFFEE-EA6A-EE42-8149-7CBEB1A3F520.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/FC2D7CE4-5223-3942-8295-C8E7E8ED1E1D.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/WZ_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/ECC01154-AC39-584E-836B-4D5B7946508D.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/70000/B9059574-2ACF-A34B-8BBE-08473A6FB46C.root'
	#'/store/mc/RunIIFall17MiniAODv2/WWToLNuQQ_NNPDF31_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/70000/C8362E91-2087-E811-A442-FA163E0726ED.root'
	#'/store/mc/RunIIFall17MiniAODv2/WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/60000/7A96A80F-E9C2-E811-AB1C-0CC47A7C357A.root'
	#'/store/mc/RunIISummer19UL17MiniAOD/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/60000/FC36084C-F598-9B40-9040-F09C968CDC9D.root'
	'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/100000/186120D6-C6A7-E811-9271-FA163E809085.root'
	#'/store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/40000/F296B7FF-629F-E811-86FE-F4E9D4AF0AF0.root'
    )
)

#process.treemaker = cms.EDAnalyzer('TTbarTauLepton')

# TreeMakder for miniAOD

process.TTbarTauLepton = cms.EDAnalyzer("TTbarTauLepton",
	monitoring        = cms.bool(True),
        monitoringHLT     = cms.bool(False),
        monitoringTau     = cms.bool(True),
        monitoringGen     = cms.bool(False),
        monitoringJets    = cms.bool(False),
        monitoringBJets   = cms.bool(False),
        monitoringLeptons = cms.bool(False),
	monitoringMET     = cms.bool(False),
	isMC = cms.bool(True),
	TauSpinnerOn = cms.bool(False),
	looseTauID = cms.bool(True),
	DeepTau = cms.bool(False),
	tauPtMin = cms.double(20),
	BJetPtMin = cms.double(30),
	MuElePtMin = cms.double(20),
	EtaMax = cms.double(2.4),
	piPtMin = cms.double(0),
	tauEtaMax = cms.double(2.3),
	tauDzMax = cms.double(0.2),
	METcut = cms.double(0),
	null = cms.double(-5),
	tauCollection = cms.string("slimmedTaus"),
	muonCollection = cms.string("slimmedMuons"),
	electronCollection = cms.string("slimmedElectrons"),
	jetCollection = cms.string("slimmedJets"),
	PuppijetCollection = cms.string("slimmedJetsPuppi"),
	metCollection = cms.string("slimmedMETs"),
	# Puppi MET
	PuppimetCollection = cms.string("slimmedMETsPuppi"),
	vertexCollection = cms.string("offlineSlimmedPrimaryVertices"),
	SVCollection     = cms.string("slimmedSecondaryVertices"),
	#genParticleCollection = cms.string("genParticles"),
	#genParticleCollection = cms.string("packedGenParticles"),
	genParticleCollection = cms.string("prunedGenParticles"),
	trackCollection = cms.string("isolatedTracks"), # generalTracks for AOD
	PackedCandidateCollection = cms.string("packedPFCandidates"),
	GenEventInfo = cms.string("generator"),
	Triggers = cms.vstring(""),
	TauTriggers = cms.vstring(
		"HLT_DoubleLooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
		"HLT_DoubleLooseChargedIsoPFTau35_Trk1_eta2p1_Reg_v",
		"HLT_DoubleLooseChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
		"HLT_DoubleLooseChargedIsoPFTau40_Trk1_eta2p1_Reg_v",
		"HLT_DoubleMediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
		"HLT_DoubleMediumChargedIsoPFTau35_Trk1_eta2p1_Reg_v",
		"HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
		"HLT_DoubleMediumChargedIsoPFTau40_Trk1_eta2p1_Reg_v",
		"HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
		"HLT_DoubleTightChargedIsoPFTau35_Trk1_eta2p1_Reg_v",
		"HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
		"HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v",
		"HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v",
		"HLT_IsoMu27_LooseChargedIsoPFTau20_SingleL1_v",
		"HLT_IsoMu27_MediumChargedIsoPFTau20_SingleL1_v",
		"HLT_IsoMu27_TightChargedIsoPFTau20_SingleL1_v",
		"HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v",
		"HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v",
		"HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v",
		"HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_v",
		"HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120_v",
		"HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_v",
		"HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v",
		"HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v",
		"HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v",
		"HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg_v",
		"HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg_v",
		),
	JetHTTriggers = cms.vstring(
		"HLT_AK8PFHT750_TrimMass50_v",
		"HLT_AK8PFHT800_TrimMass50_v",
		"HLT_AK8PFHT850_TrimMass50_v",
		"HLT_AK8PFHT900_TrimMass50_v",
		"HLT_AK8PFJet140_v",
		"HLT_AK8PFJet200_v",
		"HLT_AK8PFJet260_v",
		"HLT_AK8PFJet320_v",
		"HLT_AK8PFJet360_TrimMass30_v",
		"HLT_AK8PFJet380_TrimMass30_v",
		"HLT_AK8PFJet400_TrimMass30_v",
		"HLT_AK8PFJet400_v",
		"HLT_AK8PFJet40_v",
		"HLT_AK8PFJet420_TrimMass30_v",
		"HLT_AK8PFJet450_v",
		"HLT_AK8PFJet500_v",
		"HLT_AK8PFJet550_v",
		"HLT_AK8PFJet60_v",
		"HLT_AK8PFJet80_v",
		"HLT_AK8PFJetFwd140_v",
		"HLT_AK8PFJetFwd200_v",
		"HLT_AK8PFJetFwd260_v",
		"HLT_AK8PFJetFwd320_v",
		"HLT_AK8PFJetFwd400_v",
		"HLT_AK8PFJetFwd40_v",
		"HLT_AK8PFJetFwd450_v",
		"HLT_AK8PFJetFwd500_v",
		"HLT_AK8PFJetFwd60_v",
		"HLT_AK8PFJetFwd80_v",
		"HLT_CaloJet500_NoJetID_v",
		"HLT_CaloJet550_NoJetID_v",
		"HLT_DiPFJetAve100_HFJEC_v",
		"HLT_DiPFJetAve140_v",
		"HLT_DiPFJetAve160_HFJEC_v",
		"HLT_DiPFJetAve200_v",
		"HLT_DiPFJetAve220_HFJEC_v",
		"HLT_DiPFJetAve260_v",
		"HLT_DiPFJetAve300_HFJEC_v",
		"HLT_DiPFJetAve320_v",
		"HLT_DiPFJetAve400_v",
		"HLT_DiPFJetAve40_v",
		"HLT_DiPFJetAve500_v",
		"HLT_DiPFJetAve60_HFJEC_v",
		"HLT_DiPFJetAve60_v",
		"HLT_DiPFJetAve80_HFJEC_v",
		"HLT_DiPFJetAve80_v",
		"HLT_PFHT1050_v",
		"HLT_PFHT180_v",
		"HLT_PFHT250_v",
		"HLT_PFHT350MinPFJet15_v",
		"HLT_PFHT350_v",
		"HLT_PFHT370_v",
		"HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v",
		"HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v",
		"HLT_PFHT380_SixPFJet32_v",
		"HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v",
		"HLT_PFHT430_SixPFJet40_v",
		"HLT_PFHT430_v",
		"HLT_PFHT510_v",
		"HLT_PFHT590_v",
		"HLT_PFHT680_v",
		"HLT_PFHT780_v",
		"HLT_PFHT890_v",
		"HLT_PFJet140_v",
		"HLT_PFJet200_v",
		"HLT_PFJet260_v",
		"HLT_PFJet320_v",
		"HLT_PFJet400_v",
		"HLT_PFJet40_v",
		"HLT_PFJet450_v",
		"HLT_PFJet500_v",
		"HLT_PFJet550_v",
		"HLT_PFJet60_v",
		"HLT_PFJet80_v",
		"HLT_PFJetFwd140_v",
		"HLT_PFJetFwd200_v",
		"HLT_PFJetFwd260_v",
		"HLT_PFJetFwd320_v",
		"HLT_PFJetFwd400_v",
		"HLT_PFJetFwd40_v",
		"HLT_PFJetFwd450_v",
		"HLT_PFJetFwd500_v",
		"HLT_PFJetFwd60_v",
		"HLT_PFJetFwd80_v",
		"HLT_QuadPFJet103_88_75_15_v",
		"HLT_QuadPFJet105_88_76_15_v",
		"HLT_QuadPFJet111_90_80_15_v",
		"HLT_QuadPFJet98_83_71_15_v",
		),
	METTriggers = cms.vstring(
		"HLT_CaloMET100_HBHECleaned_v",
		"HLT_CaloMET100_NotCleaned_v",
		"HLT_CaloMET110_NotCleaned_v",
		"HLT_CaloMET250_HBHECleaned_v",
		"HLT_CaloMET250_NotCleaned_v",
		"HLT_CaloMET300_HBHECleaned_v",
		"HLT_CaloMET350_HBHECleaned_v",
		"HLT_CaloMET70_HBHECleaned_v",
		"HLT_CaloMET80_HBHECleaned_v",
		"HLT_CaloMET80_NotCleaned_v",
		"HLT_CaloMET90_HBHECleaned_v",
		"HLT_CaloMET90_NotCleaned_v",
		"HLT_CaloMHT90_v",
		"HLT_DiJet110_35_Mjj650_PFMET110_v",
		"HLT_DiJet110_35_Mjj650_PFMET120_v",
		"HLT_DiJet110_35_Mjj650_PFMET130_v",
		"HLT_L1ETMHadSeeds_v",
		"HLT_MET105_IsoTrk50_v",
		"HLT_MET120_IsoTrk50_v",
		"HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v",
		"HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
		"HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_v",
		"HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_v",
		"HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v",
		"HLT_PFMET100_PFMHT100_IDTight_PFHT60_v",
		"HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v",
		"HLT_PFMET110_PFMHT110_IDTight_v",
		"HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v",
		"HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
		"HLT_PFMET120_PFMHT120_IDTight_v",
		"HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v",
		"HLT_PFMET130_PFMHT130_IDTight_v",
		"HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v",
		"HLT_PFMET140_PFMHT140_IDTight_v",
		"HLT_PFMET200_HBHECleaned_v",
		"HLT_PFMET200_HBHE_BeamHaloCleaned_v",
		"HLT_PFMET200_NotCleaned_v",
		"HLT_PFMET250_HBHECleaned_v",
		"HLT_PFMET300_HBHECleaned_v",
		"HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v",
		"HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v",
		"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",
		"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
		"HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v",
		"HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v",
		"HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v",
		"HLT_PFMETTypeOne110_PFMHT110_IDTight_v",
		"HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v",
		"HLT_PFMETTypeOne120_PFMHT120_IDTight_v",
		"HLT_PFMETTypeOne130_PFMHT130_IDTight_v",
		"HLT_PFMETTypeOne140_PFMHT140_IDTight_v",
		"HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v",
		"HLT_TripleJet110_35_35_Mjj650_PFMET110_v",
		"HLT_TripleJet110_35_35_Mjj650_PFMET120_v",
		"HLT_TripleJet110_35_35_Mjj650_PFMET130_v",
		),
							
)

#process.load('Tau.TreeMaker.treeMaker_Data-MET_cfi')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('TTbar_MC_v1.root')
)

process.p = cms.Path(process.TTbarTauLepton)
