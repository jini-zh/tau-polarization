import FWCore.ParameterSet.Config as cms
# Define the CMSSW process
process = cms.Process("TauAnalyzer")

# Load the standard set of configuration modules

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("TauPolarization.TreeMaker.MC_v2")

from TauPolarization.TreeMaker.MC_v2 import *

#process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
#from PhysicsTools.PatAlgos.tools.tauTools import * # to get the most up-to-date discriminators when using patTaus
#switchToPFTauHPS(process)

from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag=autoCond['run2_mc']

#process.load("GeneratorInterface.TauolaInterface.TauSpinner_cfi")

process.TauSpinnerReco = cms.EDProducer( "TauSpinnerCMS",
   isReco = cms.bool(True),
   isTauolaConfigured = cms.bool(False),
   isLHPDFConfigured = cms.bool(False),
   LHAPDFname = cms.untracked.string('MSTW2008nnlo90cl.LHgrid'),
   CMSEnergy = cms.double(13000.0),
   gensrc = cms.InputTag('genParticles'),
   monitoring = cms.bool(True),
   MotherPDGID = cms.untracked.int32(24),
   Ipol = cms.untracked.int32(0),
   nonSM2 = cms.untracked.int32(0),
   nonSMN = cms.untracked.int32(0),
)

# For TauSpinner
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   TauSpinnerReco = cms.PSet(
    initialSeed = cms.untracked.uint32(123456789),
    engineName = cms.untracked.string('HepJamesRandom')
    )
                                                   )
process.randomEngineStateProducer = cms.EDProducer("RandomEngineStateProducer")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")

# How many events to process

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(10000)
)
# Define the input source

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#  'root://cms-xrd-global.cern.ch//store/mc/RunIIFall15DR76/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0058EAA8-CFC6-E511-9912-B083FED12B5C.root'
#'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/00000/006006AF-7B18-E611-A59E-0019B9CAFC2D.root'
#'file:/afs/cern.ch/user/a/aoskin/CMSSW_9_4_0/src/GeneratorInterface/TauolaInterface/test/MyTauSpinnerTests/TauSpinerRecoTest_2.root'
#'file:/eos/user/a/aoskin/Tau/MC_new/JenyaProdMC.root'
#'file:/eos/user/a/aoskin/Tau/MC_new/WToTauNu_13TeV_GEN-SIM-RAW-RECO_40k.root'
#'file:/afs/cern.ch/user/a/aoskin/TempFiles/WToTauNu_13TeV_GEN-SIM-RAW-RECO_1.root'
#WRToTauNu_191023()
WToTauNu_190830()
#ZToTauTau_191020()
#WToTauNu_hadr_nopolar_191031()
#WToTauNu_hadr_right_191102()
#WToTauNu_hadr_nopolar_191117()
        ),
                            skipEvents = cms.untracked.uint32(0)
)

process.TauAnalyzer = cms.EDAnalyzer("TauAnalyzer",
   null                = cms.int32(-5),
   #genParticleCollection = cms.InputTag("genParticles"),
   WTisValidCollection = cms.InputTag("TauSpinnerReco", "TauSpinnerWTisValid"),
   WTCollection        = cms.InputTag("TauSpinnerReco", "TauSpinnerWT"),
   WTFlipCollection    = cms.InputTag("TauSpinnerReco", "TauSpinnerWTFlip"),
   WThminusCollection  = cms.InputTag("TauSpinnerReco", "TauSpinnerWThminus"),
   WThplusCollection   = cms.InputTag("TauSpinnerReco", "TauSpinnerWThplus"),
   TauPolarisation     = cms.InputTag("TauSpinnerReco", "TauPolarisation"),
   PolarimetricVector  = cms.InputTag("TauSpinnerReco", "PolarimetricVector"),
   tauPtMin = cms.double(0),
   piPtMin = cms.double(0),
   tauEtaMax = cms.double(2.3),
   tauDzMax = cms.double(0.5),
   Triggers = cms.vstring("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v1",
      "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v2",
      "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v3",
      "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v4",
      "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v5",
      "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v6",
      "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v7",
      "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v8",
      "HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90_v9"
   ),
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
)

#process.load('Tau.TauAnalyzer.TauAnalyzer_MC_cfi')

process.p = cms.Path(process.TauSpinnerReco + process.TauAnalyzer)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(
					#"/afs/cern.ch/user/a/aoskin/CMSSW_9_4_0/src/Tau/TauAnalyzer/test/TauSpinAnalyzer_MC_40k_init.root"
					"W_Left_ZhenyaMC_10k.root"
					)
                                   )
