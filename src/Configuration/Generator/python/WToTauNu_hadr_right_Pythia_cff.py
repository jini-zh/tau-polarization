import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *



generator = cms.EDFilter("Pythia8GeneratorFilter",
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
        'WeakSingleBoson:ffbar2W = on', # enable W production
        '24:onMode = off',              # disable W decays
        '24:onIfAny = 15 16',           # enable W -> tau nu
        '24:mMin = 50.',                # minimal invariant mass of W decay products
        '15:offIfAny = 11 13',          # disable tau -> mu nu and tau -> e nu decays
        'TauDecays:mode = 2',              # set tau polarization
        'TauDecays:tauPolarization = 1.0', #   to 1.0
        'TauDecays:tauMother = 24'         #   for taus in W decays
        ),
      parameterSets = cms.vstring('pythia8CommonSettings',
        'pythia8CUEP8M1Settings',
        'processParameters',
        )
      )
    )

ProductionFilterSequence = cms.Sequence(generator)
