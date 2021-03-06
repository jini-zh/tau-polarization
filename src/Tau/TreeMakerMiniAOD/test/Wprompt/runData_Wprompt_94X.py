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

# Single Electron 2016 17 Jul 2018 data
process.GlobalTag.globaltag='94X_dataRun2_v10'



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/30000/A69A903E-4EF2-CE4D-AD5E-7FCAEE9F8326.root',
    #'/store/data/Run2017B/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/B6EC130B-9D55-8140-AC2E-154EB59E3BBC.root'
    #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/100000/186120D6-C6A7-E811-9271-FA163E809085.root'
    #'/store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/40000/F296B7FF-629F-E811-86FE-F4E9D4AF0AF0.root'
    'root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleMuon/MINIAOD/31Mar2018-v1/90000/E40B3455-5939-E811-94C5-0CC47A7C34C8.root'
    #'root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/30000/E2A2C81B-0638-E811-89CE-008CFAC93EA8.root'
    #'root://cms-xrd-global.cern.ch//store/data/Run2017D/Tau/MINIAOD/31Mar2018-v1/00000/12F092A0-3F37-E811-AD27-7845C4F92C96.root'
    #'root://cms-xrd-global.cern.ch//store/data/Run2016C/SingleMuon/MINIAOD/17Jul2018-v1/20000/065D2BF3-9198-E811-844D-90E2BAC9B7A8.root'
    #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/100000/186120D6-C6A7-E811-9271-FA163E809085.root'
    ),
    dropDescendantsOfDroppedBranches=cms.untracked.bool(False),
)

#process.treemaker = cms.EDAnalyzer('TTbarTauLepton')

# MVA MET

#process.load("RecoJets.JetProducers.ak4PFJets_cfi")
#process.task.add(process.ak4PFJets)
#process.ak4PFJets.src = cms.InputTag("packedPFCandidates")
#process.ak4PFJets.doAreaFastjet = cms.bool(True)

#from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3

#process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
#process.task.add(process.pfMVAMEtTask)
#process.MVAMET = process.pfMVAMEtTask
#process.pfMVAMEt.srcLeptons = cms.VInputTag("slimmedElectrons")
#process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
#process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")

# TreeMakder for miniAOD

process.TTbarTauLepton = cms.EDAnalyzer("TTbarTauLepton",
    monitoring        = cms.bool(True),
    monitoringHLT     = cms.bool(False),
    monitoringTau     = cms.bool(False),
    monitoringGen     = cms.bool(False),
    monitoringJets    = cms.bool(True),
    monitoringBJets   = cms.bool(True),
    monitoringLeptons = cms.bool(True),
    monitoringMET     = cms.bool(False),
    isMC = cms.bool(False),
    #fullMC = cms.bool(False),
    #useHLT = cms.bool(False),
    #useTargetHLT = cms.bool(False),
    TauSpinnerOn = cms.bool(False),
    looseTauID = cms.bool(True),
    DeepTau = cms.bool(False),
    tauPtMin = cms.double(20),
    BJetPtMin = cms.double(30),
    MuElePtMin = cms.double(20),
    EtaMax = cms.double(2.4),
    piPtMin = cms.double(0),
    tauEtaMax = cms.double(2.4),
    tauDzMax = cms.double(0.2),
    METcut = cms.double(0),
    null = cms.double(-10),
    # Parameters added at January 2021
    JetEtaMax = cms.double(2.4),
    LeptonRequired = cms.bool(False),
    #NrequiredJets = cms.int32(2),
    #NrequiredBJets = cms.int32(1),
    #NrequiredLeptons = cms.int32(1),
    #requiredLeptonPDGID = cms.int32(13),
    UseTau = cms.bool(False),
    ####
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
    # for 2017 data
    Triggerobjects = cms.InputTag('slimmedPatTrigger', '', 'PAT'),
    prescales      = cms.InputTag('patTrigger', '', 'PAT'),
    prescalesL1min = cms.InputTag('patTrigger', 'l1min', 'PAT'),
    prescalesL1max = cms.InputTag('patTrigger', 'l1max', 'PAT'),
    # for 2016 data
    #Triggerobjects = cms.InputTag('slimmedPatTrigger', '', 'DQM'),
    #prescales      = cms.InputTag('patTrigger', '', 'DQM'),
    #prescalesL1min = cms.InputTag('patTrigger', 'l1min', 'DQM'),
    #prescalesL1max = cms.InputTag('patTrigger', 'l1max', 'DQM'),
    Triggers = cms.vstring(""),
    TriggerTarget1 = cms.vstring(
        "HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
    ),
    TriggerTarget2 = cms.vstring(
        "HLT_PFMET120_PFMHT120_IDTight_v",
    ),
    TriggerTarget3 = cms.vstring(
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
    ),
    # MET 2017 part 1
    Triggers1 = cms.vstring(
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
    ),
    # MET 2017 part 2
    Triggers2 = cms.vstring(
        "HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v",
        "HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v",
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
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",
        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
        "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v",
        "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v",
        "HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v",
        "HLT_PFMETTypeOne110_PFMHT110_IDTight_v",
        "HLT_PFMETTypeOne120_PFMHT120_IDTight_HFCleaned_v",
        "HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v",
        "HLT_PFMETTypeOne120_PFMHT120_IDTight_v",
        "HLT_PFMETTypeOne130_PFMHT130_IDTight_v",
        "HLT_PFMETTypeOne140_PFMHT140_IDTight_v",
        "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v",
        "HLT_TripleJet110_35_35_Mjj650_PFMET110_v",
        "HLT_TripleJet110_35_35_Mjj650_PFMET120_v",
        "HLT_TripleJet110_35_35_Mjj650_PFMET130_v",
    ),
    # Tau 2017 part 
    Triggers3 = cms.vstring(
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
        "HLT_MediumChargedIsoPFTau100HighPtRelaxedIso_Trk50_eta2p1_1pr_v",
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
    # SingleElectron 2017 part 2
    Triggers4 = cms.vstring(
        "",
    ),
    Triggers5 = cms.vstring(
        "",
    ),     
)

#process.load('Tau.TreeMaker.treeMaker_Data-MET_cfi')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('Data_SingleMuon.root')
)

process.p = cms.Path(process.TTbarTauLepton)
