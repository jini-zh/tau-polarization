# Adapted from src/Configuration/Generator/python/TTbarLepton_13TeV_Tune_CUEP8M1_cfi.py

import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

generator = cms.EDFilter(
    'Pythia8GeneratorFilter',
    pythiaHepMCVerbosity  = cms.untracked.bool(False),
    maxEventsToPrint      = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency      = cms.untracked.double(1.0),
    comEnergy             = cms.double(13000.0),
    PythiaParameters = cms.PSet(
      pythia8CommonSettingsBlock,
      pythia8CUEP8M1SettingsBlock,
      processParameters = cms.vstring(
        'Top:gg2ttbar = on ',              # gg -> t tbar
        '6:m0 = 175 ',                     # set t mass to 175 GeV
        '24:onMode = off',                 # disable W decays
        '24:onIfAny = 11 12',              # enable W -> e nu_e
        '24:onIfAny = 13 14',              # enable W -> mu nu_mu
        '24:onIfAny = 15 16',              # enable W -> tau nu_tau
        '15:offIfAny = 11 13',             # disable tau -> mu nu_mu and tau->e nu_e
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
