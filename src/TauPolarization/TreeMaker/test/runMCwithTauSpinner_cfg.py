import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMaker")


#process.load("Configuration.StandardSequences.Services_cff")
#process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("GeneratorInterface.TauolaInterface.TauSpinner_cfi")

#process.GlobalTag.globaltag='80X_dataRun2_Prompt_v14'
#process.GlobalTag.globaltag='80X_dataRun2_2016SeptRepro_v3'
#process.GlobalTag.globaltag='80X_dataRun2_2016LegacyRepro_v4'
#process.GlobalTag.globaltag='106X_dataRun2_v20'
#process.GlobalTag.globaltag='94X_dataRun2_ReReco17_forValidation'
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag=autoCond['run2_mc']

# For TauSpinner
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   TauSpinnerReco = cms.PSet(
    initialSeed = cms.untracked.uint32(123456789),
    engineName = cms.untracked.string('HepJamesRandom')
    )
                                                   )
process.randomEngineStateProducer = cms.EDProducer("RandomEngineStateProducer")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#'root://cms-xrd-global.cern.ch//store/data/Run2016H/Tau/AOD/PromptReco-v2/000/281/207/00000/CA27ACAD-6282-E611-BF27-FA163ED9B340.root'
      #'root://cms-xrd-global.cern.ch//store/data/Run2016D/Tau/AOD/23Sep2016-v1/100000/000A4C44-C391-E611-B70D-001E67E6F4CC.root'
  #  'root://cms-xrd-global.cern.ch//store/data/Run2016B/JetHT/AOD/22Feb2017_ver2-v1/00000/00335C56-32FA-E611-B230-0025905B85F6.root'
	#'root://cms-xrd-global.cern.ch//store/data/Run2016G/Tau/AOD/PromptReco-v1/000/278/820/00000/004A9416-2364-E611-9E82-02163E011C51.root'
	#'root://srmcms.pic.es//store/data/Run2016G/Tau/AOD/PromptReco-v1/000/278/820/00000/004A9416-2364-E611-9E82-02163E011C51.root'
	#set=/Tau*/Run201*/AOD'/store/data/Run2016D/Tau/AOD/PromptReco-v2/000/276/315/00000/1448B5D3-FC44-E611-85BF-02163E014522.root'
        #'root://cms-xrd-global.cern.ch/store/data/Run2016H/Tau/MINIAOD/18Apr2017-v1/10000/0EA13355-B650-E711-AA53-D4856445C5F8.root'
        #'/store/data/Run2016D/Tau/MINIAOD/23Sep2016-v1/90000/F037EAD8-D691-E611-B319-0CC47A4C8E86.root'

        # Data which works well with 80X_dataRun2_2016LegacyRepro_v4 gloabal tag
        # '/store/data/Run2016H/Tau/AOD/07Aug17-v1/90001/98DBABD5-E39C-E711-8FD2-0242AC130002.root'
        # Data with 106X_dataRun2_v20 global tag (some problems now)
        #'root://cms-xrd-global.cern.ch/store/data/Run2017D/Tau/AOD/09Aug2019_UL2017-v1/60000/21459DD8-DD28-984D-93E4-9757385D236F.root'
        #'/store/data/Run2017D/Tau/AOD/09Aug2019_UL2017-v1/60000/1A854349-E00B-4845-A597-5E8EFC258DE9.root'
        #'root://cms-xrd-global.cern.ch//store/data/Run2017D/Tau/AOD/17Nov2017-v1/00000/4C06B693-D2F2-E711-AF43-A4BF011253C0.root'
        #'root://cms-xrd-global.cern.ch//store/data/Run2017D/Tau/AOD/17Nov2017-v1/00000/7261CC72-75ED-E711-995D-1CC1DE1CF622.root'
	#'root://cms-xrd-global.cern.ch//store/data/Run2017D/Tau/AOD/17Nov2017-v1/710000/0A46559A-AAF4-E711-91BE-0242AC1C0500.root'
        #'/store/data/Run2017C/Tau/AOD/17Nov2017-v1/70001/EE79F9C6-F8E5-E711-BC9D-0CC47A5FBDC1.root'
	#'/store/mc/RunIIFall17DRPremix/ZHToTauTau_M125_13TeV_powheg_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/D6CF698A-5427-E811-8463-3417EBE64888.root'
	#'file:/afs/cern.ch/user/a/aoskin/TempFiles/TauSpinerRecoTest_1.root'
	'file:/afs/cern.ch/user/a/aoskin/TempFiles/WToTauNu_13TeV_GEN-SIM-RAW-RECO_1.root'

    )
)

process.treemaker = cms.EDAnalyzer("TreeMaker",
                                                        isMC = cms.bool(True),
							printOutput = cms.bool(True),
							monitoring = cms.bool(False),
							tauPtMin = cms.double(20),
                                                        piPtMin = cms.double(0),
                                                        tauEtaMax = cms.double(2.3),
                                                        tauDzMax = cms.double(0.2),
                                                        null = cms.double(-5),
                                                        Triggers = cms.vstring(""),
                                                        tauCollection = cms.string("hpsPFTauProducer"),
                                                        muonCollection = cms.string("muons"),
                                                        electronCollection = cms.string("gedGsfElectrons"),
                                                        jetCollection = cms.string("ak4PFJets"),
                                                        metCollection = cms.string("pfMet"),
                                                        vertexCollection = cms.string("offlinePrimaryVertices"),
                                                        genParticleCollection = cms.string("genParticles"),
                                                        trackCollection = cms.string("generalTracks"),
                                                        absIsoDiscriminator = cms.string("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits"),
                                                        looseCombinedIsoDiscriminator = cms.string("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),
                                                        mediumCombinedIsoDiscriminator = cms.string("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits"),
                                                        tightCombinedIsoDiscriminator = cms.string("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits"),
                                                        looseMvaIsoDiscriminator = cms.string("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT"),
                                                        mediumMvaIsoDiscriminator = cms.string("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT"),
                                                        tightMvaIsoDiscriminator = cms.string("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT"),
                                                        looseMuonRejectionDiscriminator = cms.string("hpsPFTauDiscriminationByLooseMuonRejection3"),
                                                        tightMuonRejectionDiscriminator = cms.string("hpsPFTauDiscriminationByTightMuonRejection3"),
                                                        looseElectronRejectionDiscriminator = cms.string("hpsPFTauDiscriminationByMVA6VLooseElectronRejection"),
                                                        tightElectronRejectionDiscriminator = cms.string("hpsPFTauDiscriminationByMVA6VTightElectronRejection"),
							WTisValidCollection = cms.InputTag("TauSpinnerReco", "TauSpinnerWTisValid"),
							WTCollection        = cms.InputTag("TauSpinnerReco", "TauSpinnerWT"),
							WTFlipCollection    = cms.InputTag("TauSpinnerReco", "TauSpinnerWTFlip"),
							WThminusCollection  = cms.InputTag("TauSpinnerReco", "TauSpinnerWThminus"),
							WThplusCollection   = cms.InputTag("TauSpinnerReco", "TauSpinnerWThplus"),
                                                        MotherCollection    = cms.InputTag("TauSpinnerReco", "TauSpinnerMother"),
)

#process.treemaker = cms.EDAnalyzer('TreeMaker'
#)
#
#process.load('Tau.TreeMaker.treeMaker_Data-MET_cfi')

process.TFileService = cms.Service("TFileService",
  fileName = cms.string(
                       'MCTauspinner_test.root'
  )
)

process.p = cms.Path(process.TauSpinnerReco + process.treemaker)
