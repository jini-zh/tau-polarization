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



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/data/Run2017C/Tau/MINIAOD/09Aug2019_UL2017-v1/30000/A69A903E-4EF2-CE4D-AD5E-7FCAEE9F8326.root',
    #'/store/data/Run2017B/JetHT/MINIAOD/09Aug2019_UL2017-v1/130000/B6EC130B-9D55-8140-AC2E-154EB59E3BBC.root'
    #'/store/mc/RunIISummer19UL17MiniAOD/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/100000/E059A7E6-BA71-DB42-BF0A-2381E40C21AB.root'
    #'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/100000/186120D6-C6A7-E811-9271-FA163E809085.root'
    '/store/mc/RunIIFall17MiniAODv2/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/40000/F296B7FF-629F-E811-86FE-F4E9D4AF0AF0.root'
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
    monitoringJets    = cms.bool(False),
    monitoringBJets   = cms.bool(False),
    monitoringLeptons = cms.bool(False),
    monitoringMET     = cms.bool(False),
    isMC = cms.bool(True),
    fullMC = cms.bool(True),
    useHLT = cms.bool(False),
    useTargetHLT = cms.bool(False),
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
    NrequiredJets = cms.int32(0),
    NrequiredBJets = cms.int32(-1),
    NrequiredLeptons = cms.int32(-1),
    requiredLeptonPDGID = cms.int32(0),
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
    # 2017 data triggers
    TriggerTarget = cms.vstring(
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v",
    ),
    TriggerTarget2 = cms.vstring(
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1_v",
    ),
    TriggersTarget = cms.vstring(
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v",
        "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1_v",
    ),
    TriggersTargetPrescaled = cms.vstring(
        "emptyHLT",
    ),
    TriggersUnprescaled = cms.vstring(
        "HLT_Ele115_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele135_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele145_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele200_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele250_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele27_WPTight_Gsf_v",
        "HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v",
        "HLT_Ele300_CaloIdVT_GsfTrkIdT_v",
        "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v",
        "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v",
        "HLT_Ele32_WPTight_Gsf_v",
        "HLT_Ele35_WPTight_Gsf_L1EGMT_v",
        "HLT_Ele35_WPTight_Gsf_v",
        "HLT_Ele38_WPTight_Gsf_v",
        "HLT_Ele40_WPTight_Gsf_v",
        "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v",

        ),
    TriggersPrescaled = cms.vstring(
        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",
        "HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v",
        "HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v",
        "HLT_Ele15_IsoVVVL_PFHT450_v",
        "HLT_Ele15_IsoVVVL_PFHT600_v",
        "HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",
        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v",
        "HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v",
        "HLT_Ele50_IsoVVVL_PFHT450_v",
        "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",
        "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v",
        ),     
)

#process.load('Tau.TreeMaker.treeMaker_Data-MET_cfi')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('Wprompt_PU.root')
)

process.p = cms.Path(process.TTbarTauLepton)
