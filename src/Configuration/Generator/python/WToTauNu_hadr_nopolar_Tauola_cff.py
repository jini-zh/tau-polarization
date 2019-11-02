import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from GeneratorInterface.ExternalDecays.TauolaSettings_cff import *

generator = cms.EDFilter(
    "Pythia8GeneratorFilter",
    comEnergy = cms.double(13000.0),
    crossSection = cms.untracked.double(6.44),
    filterEfficiency = cms.untracked.double(1),
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    PythiaParameters = cms.PSet(
      pythia8CommonSettingsBlock,
      pythia8CUEP8M1SettingsBlock,
      processParameters = cms.vstring(
        'WeakSingleBoson:ffbar2W = on',
        '24:onMode = off',    # disable W decays
        '24:onIfAny = 15 16', # enable W -> tau nu decay
        '24:mMin = 50.',      # minimal invariant mass of decay products
        '15:onMode = off'     # disable tau decays
        ),
      parameterSets = cms.vstring('pythia8CommonSettings',
        'pythia8CUEP8M1Settings',
        'processParameters',
      )
    ),
    ExternalDecays = cms.PSet(
      Tauola = cms.untracked.PSet(
        TauolaNoPolar,           # make unpolarized tau
        InputCards = cms.PSet(
          pjak1 = cms.int32(0),
          pjak2 = cms.int32(0),
          mdtau = cms.int32(130) # only hadron decays for tau
        )
      ),
      parameterSets = cms.vstring('Tauola')
    )
)


ProductionFilterSequence = cms.Sequence(generator)
