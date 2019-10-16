import re
import argparse

import XRootD.client

import FWCore.ParameterSet.Config as cms

argparser = argparse.ArgumentParser()
argparser.add_argument('pset')
argparser.add_argument('files', nargs='+')
args = argparser.parse_args()

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

files = []
for f in args.files:
  if re.match('file:', f):
    files.append(f)
    continue

  m = re.match('root://([^/]+)/(.*)/', f)
  if m:
    host = m.group(1)
    path = [ '' ]
    for directory in re.group(2).split('/'):
      if directory == '': continue
      rd = re.compile(d + '$')
      newpath = []
      for p in path:
        status, dirlist = XRootD.client.FileSystem(host).dirlist(p)
        for x in dirlist:
          if rd.match(x):
            newpath.append(x)
      if len(newpath) == 0:
        if len(path) == 1:
          raise Exception('Path root://%s/%s/%s does not exist' % host % path[0] % directory)
        else:
          msg = 'None of the pathes\n'
          for p in path:
            msg += 'root://%s/%s/%s\n' % host % p % directory
          msg += 'exist'
          raise Exception(msg)
      path = newpath
    files += path
    continue

  files.append('file:' + f)

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
