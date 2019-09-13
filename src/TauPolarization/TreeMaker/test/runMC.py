import FWCore.ParameterSet.Config as cms
# Define the CMSSW process
process = cms.Process("TreeMaker")

# Load the standard set of configuration modules

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

#process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
#from PhysicsTools.PatAlgos.tools.tauTools import * # to get the most up-to-date discriminators when using patTaus
#switchToPFTauHPS(process)

from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag=autoCond['run2_mc']



# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")

# How many events to process

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(100)
)
# Define the input source

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#  'root://cms-xrd-global.cern.ch//store/mc/RunIIFall15DR76/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0058EAA8-CFC6-E511-9912-B083FED12B5C.root'
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/00000/006006AF-7B18-E611-A59E-0019B9CAFC2D.root'
        ),
                            skipEvents = cms.untracked.uint32(0)
)

process.treemaker = cms.EDAnalyzer("TreeMaker")

process.load('Tau.TreeMaker.treeMaker_MC-April_cfi')

process.p = cms.Path(process.treemaker)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test_MC_aod.root")
                                   )
