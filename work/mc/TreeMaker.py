import re

import XRootD.client

import FWCore.ParameterSet.Config as cms
# Define the CMSSW process
process = cms.Process("TreeMaker")

# Load the standard set of configuration modules

process.load('FWCore.MessageLogger.MessageLogger_cfi')
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



## Message Logger settings
#process.load("FWCore.MessageService.MessageLogger_cfi")

# How many events to process

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(-1)
)
# Define the input source

source = '/data/phedex/user/ezhemchu/crab/WToTauNu_13TeV/WToTauNu_13TeV/crab_WToTauNu_13TeV/190830_120718/0000'
status, dirlist = XRootD.client.FileSystem('cms-phedex.lxfarm.mephi.ru').dirlist(source)
files = [ 'root://cms-phedex.lxfarm.mephi.ru' + source + '/' + f.name for f in dirlist if re.search(r'GEN-SIM-RAW-RECO_\d+\.root$', f.name) ]

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(*files),
                            skipEvents = cms.untracked.uint32(0)
)

process.treemaker = cms.EDAnalyzer("TreeMaker")

process.load('TauPolarization.TreeMaker.treeMaker_MC-April_cfi')

process.p = cms.Path(process.treemaker)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test_MC_aod.root")
                                   )

#process.MessageLogger = cms.Service(
#    'MessageLogger',
#    destinations = cms.untracked.vstring('debug'),
#    debugModules = cms.untracked.vstring('treemaker'),
#    debug = cms.untracked.PSet(
#      threshold = cms.untracked.string('DEBUG'),
#    )
#)
