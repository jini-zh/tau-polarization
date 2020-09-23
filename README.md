# TauAnalysis

Current event selection:
- One or more unprescaled triggers from Tau or JetHT or MET or BTagCSV or SingleEletron or SingleMuon
- One (or more) tau-lepton (Pt > 20 GeV, |Eta| < 2.3, dZ(pv) < 0.2) satisfying at least following discriminators:
  MVA:
  * byVVLooseIsolationMVArun2v1DBoldDMwLT
  * againstElectronVLooseMVA6
  * againstMuonLoose3
  or DeepTau:
  * byVLooseDeepTau2017v2VSmu
  * byVVVLooseDeepTau2017v2VSe
  * byVVVLooseDeepTau2017v2VSjet
- At least 2 loose btagged Puppi jets (Pt > 30 GeV, |Eta| < 2.4)
- At least 1 muon or electron (Pt > 20 GeV, |Eta| < 2.4)

Output variables:
- Tau parameters
  * "tau_dm" is old decay mode (i.e. 0 for 1 track 0 pi0, 1 for 1 track n pi0, 10 and 11 for 3 tracks)
  * It is preferably to use "oldDM" ones as Raw Iso discriminators (for example "tau_IsoMVArun2v1DBoldDMwLTraw")
  * "tau_againstElectronRaw" for raw discrimination against electorns
  * Only fixed working points disriminator aganst muons
  * basic parameters of second tau
  * Leading cahrged pion parameters "PiChar_*" (and only pi charged if dm is 0 or 1)
  * "piZero_*" varibles obtained from difference of tau and PiChar four-vectors (OK for dm = 1)
  * "DiPhoton_*" variables is the parameter of photons pair with closest to Pi0 invariant mass
  * "tau_hasSV" (secondary vertex) > 0 if there is SV close to tau tracks source
    "tau_SVdR" is a distance (dR = sqrt(deltaEta^2 + deltaPhi^2)) between SV and track vertex closest to SV
    (do not use SV information if the value of "tau_SVdR" is large)
    "pv_SVdR" is distance between primary vertex and secondary vertex associated with tau tracks
    "tauPVtoSVdPhi" is azimuth angle between tau 3-vector and 3-vector which connects PV and SV of tau
  * "tau_dR2", tau_dR2 = sum(dR^2 * pT^2)/sum(pT^2), sum over tau decay products (tracks and photons)
- Parameters of generated particles
  * "dR" is a distance between reconstructed tau and generated tau lepton closest to reconstructed tau
    "leptondR", "bquark1dR", "bquark2dR" the same for lepton and 2 b-quarks
- BJets and 2 leading jets parameters
- Two leptons with highest pT (and dR(tau, lepton) > 0.3) parameters:
  * "*Iso" parameters is sum of pT of particles in cone (dR < 0.4 (?)) around lepton
    So higher "*Iso" value corresponds to less isolated particle
- Two types of MET (may be Puppi MET is better as PU subtracted one).
  * "dPhi" is asimuth angle between MET and tau
  * "dPhiPuppimetTau" is asimuth angle between Puppi MET and tau
- "n*Triggers" is a number of particular data type triggers (JetHT or MET or BTagCSV ...) in event
