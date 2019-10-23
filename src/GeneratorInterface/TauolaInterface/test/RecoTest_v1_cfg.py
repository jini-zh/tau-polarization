import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST1")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions#Global_Tags_for_Monte_Carlo_Prod
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("GeneratorInterface.TauolaInterface.TauSpinner_cfi")

process.MessageLogger.cerr = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
    FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(0)),
    DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(-1))
    )

numberOfEvents = 10

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   TauSpinnerReco = cms.PSet(
    initialSeed = cms.untracked.uint32(123456789),
    engineName = cms.untracked.string('HepJamesRandom')
    )
                                                   )
process.randomEngineStateProducer = cms.EDProducer("RandomEngineStateProducer")

#process.GlobalTag.globaltag = 'MC_50_V13::All'
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag=autoCond['run2_mc']

# List of input files
Files = []
S1 = 'root://cms-phedex.lxfarm.mephi.ru/data/phedex/user/ezhemchu/crab/WToTauNu_13TeV/WToTauNu_13TeV/crab_WToTauNu_13TeV/190830_120718/0000/WToTauNu_13TeV_GEN-SIM-RAW-RECO_'
S3 = '.root'
for i in range(41):
  if i == 0 or i == 26 or i == 27 or i == 38:
    continue
  S2 = str(i)
  string = S1 + S2 + S3
  Files.append(string)

process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring(
#'file:/afs/cern.ch/user/i/inugent/tmp/5C3DF315-CF96-E111-9323-0025B3E05BF4.root'
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/C2622FE9-5DBD-E611-B78D-5065F38122D1.root'
#'file:/afs/cern.ch/user/a/aoskin/TempFiles/TauSpinerRecoTest_1.root'
#'file:/afs/cern.ch/user/a/aoskin/TempFiles/WToTauNu_13TeV_GEN-SIM-RAW-RECO_1.root'
#'file:/afs/cern.ch/user/a/aoskin/CMSSW_9_4_0/src/GeneratorInterface/TauolaInterface/test/MyTauSpinnerTests/Test_Sample.root'
#Files
#'/store/mc/RunIISummer16DR80Premix/ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/50000/42F04128-54E0-E611-83C0-FA163E0C309F.root'
# there are some events with taus
'/store/mc/RunIIFall17DRPremix/ZHToTauTau_M125_13TeV_powheg_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/70000/D6CF698A-5427-E811-8463-3417EBE64888.root'
)
)

process.debugOutput = cms.OutputModule("PoolOutputModule",
                                       outputCommands = cms.untracked.vstring('keep *'),
                                       fileName = cms.untracked.string(
                                       #'/eos/user/a/aoskin/Tau/MC_new/JenyaProdMC.root'
                                       'TauSpinner_PublicMC_1.root'
                                       ),
                                       )
process.out_step = cms.EndPath(process.debugOutput)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(numberOfEvents) )
process.p1 = cms.Path(process.TauSpinnerReco )
process.schedule = cms.Schedule(process.p1)
process.schedule.append(process.out_step)
