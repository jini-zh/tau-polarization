// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATTauDiscriminator.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

#include "RecoTauTag/RecoTau/interface/RecoTauBuilderPlugins.h"
#include "RecoTauTag/RecoTau/interface/RecoTauCommonUtilities.h"
#include "RecoTauTag/RecoTau/interface/PFTauDecayModeTools.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Association.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"

#include <cstdio>
#include <iostream>
#include <math.h>

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Definition of class LeptonCandidate
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
class LeptonCandidate {
public :
    //--- variables
    double Pt;
    double Eta, Phi;
    int Flavor;
    int Charge;
    // Methods from pat::Muon
    float caloIso;
    float ecalIso;
    float hcalIso;
    float puppiChargedHadronIso;
    float puppiNeutralHadronIso;
    float puppiPhotonIso;
    float puppiNoLeptonsChargedHadronIso;
    float puppiNoLeptonsNeutralHadronIso;
    float puppiNoLeptonsPhotonIso;
    float trackIso;
    math::XYZTLorentzVector FourMomentum;
    double Dz;
    double tauAbsIso;

    //--- constructor & destructor
    LeptonCandidate();
    LeptonCandidate(const pat::Muon &muon, math::XYZPoint pv_position) 
        : Pt(muon.pt()), Eta(muon.eta()), Phi(muon.phi()), Flavor(13), Charge(muon.charge()), caloIso(muon.caloIso()), ecalIso(muon.ecalIso()), hcalIso(muon.hcalIso()),
        puppiChargedHadronIso(muon.puppiChargedHadronIso()), puppiNeutralHadronIso(muon.puppiNeutralHadronIso()), puppiPhotonIso(muon.puppiPhotonIso()),
        puppiNoLeptonsChargedHadronIso(muon.puppiNoLeptonsChargedHadronIso()), puppiNoLeptonsNeutralHadronIso(muon.puppiNoLeptonsNeutralHadronIso()), puppiNoLeptonsPhotonIso(muon.puppiNoLeptonsPhotonIso()),
        trackIso(muon.trackIso()), FourMomentum(muon.p4()), Dz((pv_position - muon.vertex()).R()), tauAbsIso(-10)
    {
    }
    LeptonCandidate(const pat::Electron &electron, math::XYZPoint pv_position) 
        : Pt(electron.pt()), Eta(electron.eta()), Phi(electron.phi()), Flavor(11), Charge(electron.charge()), caloIso(electron.caloIso()), ecalIso(electron.ecalIso()), hcalIso(electron.hcalIso()),
        puppiChargedHadronIso(electron.puppiChargedHadronIso()), puppiNeutralHadronIso(electron.puppiNeutralHadronIso()), puppiPhotonIso(electron.puppiPhotonIso()),
        puppiNoLeptonsChargedHadronIso(electron.puppiNoLeptonsChargedHadronIso()), puppiNoLeptonsNeutralHadronIso(electron.puppiNoLeptonsNeutralHadronIso()), puppiNoLeptonsPhotonIso(electron.puppiNoLeptonsPhotonIso()),
        trackIso(electron.trackIso()), FourMomentum(electron.p4()), Dz((pv_position - electron.vertex()).R()), tauAbsIso(-10)
    {
    } 
    LeptonCandidate(const pat::Tau &tau, math::XYZPoint pv_position) 
        : Pt(tau.pt()), Eta(tau.eta()), Phi(tau.phi()), Flavor(15), Charge(tau.charge()), caloIso(-10), ecalIso(-10), hcalIso(-10),
        puppiChargedHadronIso(-10), puppiNeutralHadronIso(-10), puppiPhotonIso(-10),
        puppiNoLeptonsChargedHadronIso(-10), puppiNoLeptonsNeutralHadronIso(-10), puppiNoLeptonsPhotonIso(-10),
        trackIso(-10), FourMomentum(tau.p4()), Dz((pv_position - tau.vertex()).R()), tauAbsIso(tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"))
    {
    }
    void Monitoring();
    //double LepdR (double eta, double phi);
    virtual ~LeptonCandidate();

    //--- functions

};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Constructor of LeptonCandidate object
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
LeptonCandidate::LeptonCandidate() {

}

//LeptonCandidate::LeptonCandidate(const pat::Muon &muon) {
//
//}

LeptonCandidate::~LeptonCandidate()
{
// Destructor
}

// Show lepton parameters

void LeptonCandidate::Monitoring() {

    std::cout << "flavor = " << Flavor << std::endl
    << "charge = " << Charge << std::endl
    << "dz     = " << Dz << std::endl
    << "pt     = " << Pt << std::endl
    << "eta    = " << Eta << std::endl
    << "phi    = " << Phi << std::endl
    << "caloIso                        = " << caloIso << std::endl
    << "ecalIso                        = " << ecalIso << std::endl
    << "hcalIso                        = " << hcalIso << std::endl
    << "puppiChargedHadronIso          = " << puppiChargedHadronIso << std::endl
    << "puppiNeutralHadronIso          = " << puppiNeutralHadronIso << std::endl
    << "puppiPhotonIso                 = " << puppiPhotonIso << std::endl
    << "puppiNoLeptonsChargedHadronIso = " << puppiNoLeptonsChargedHadronIso << std::endl
    << "puppiNoLeptonsNeutralHadronIso = " << puppiNoLeptonsNeutralHadronIso << std::endl
    << "puppiNoLeptonsPhotonIso        = " << puppiNoLeptonsPhotonIso << std::endl
    << "trackIso                       = " << trackIso << std::endl;
}

/*
double LeptonCandidate::LepdR(double eta, double phi, ) {

    << "delta(BJet1) = " << sqrt(deltaR2(BJet1_eta, BJet1_phi, Eta, Phi)) << std::endl
    << "delta(BJet2) = " << sqrt(deltaR2(BJet2_eta, BJet2_phi, Eta, Phi)) << std::endl
    << "delta(Tau)   = " << sqrt(deltaR2(tau_eta, tau_phi, Eta, Phi)) << std::endl

}
*/

// Algos for leptons sorting

void SwapLeptons(std::vector<LeptonCandidate*> items, int left, int right) {
  if (left != right) {
    LeptonCandidate* temp = items[left];
    items[left] = items[right];
    items[right] = temp;
  }
}

void SortLeptons(std::vector<LeptonCandidate>& items) {
  bool swapped;
  do {
    swapped = false;
    for (unsigned i = 1; i < items.size(); i++) {
      if (items[i-1].Pt < items[i].Pt) {
        std::swap(items[i-1], items[i]);
        swapped = true;
      }
    }
  } while (swapped != false);
}

void SortElectrons(std::vector<pat::Electron>& items) {
  bool swapped;
  do {
    swapped = false;
    for (unsigned i = 1; i < items.size(); i++) {
      if (items[i-1].pt() < items[i].pt()) {
        std::swap(items[i-1], items[i]);
        swapped = true;
      }
    }
  } while (swapped != false);
}

void SortMuons(std::vector<pat::Muon>& items) {
  bool swapped;
  do {
    swapped = false;
    for (unsigned i = 1; i < items.size(); i++) {
      if (items[i-1].pt() < items[i].pt()) {
        std::swap(items[i-1], items[i]);
        swapped = true;
      }
    }
  } while (swapped != false);
}

void SortTaus(std::vector<pat::Tau>& items) {
  bool swapped;
  do {
    swapped = false;
    for (unsigned i = 1; i < items.size(); i++) {
      if (items[i-1].pt() < items[i].pt()) {
        std::swap(items[i-1], items[i]);
        swapped = true;
      }
    }
  } while (swapped != false);
}

void printElectronCandidate(pat::Electron &electron, bool PrivateMC) {
    std::cout << "Electron parameters" << std::endl
    << "pt     = " << electron.pt() << std::endl
    << "eta    = " << electron.eta() << std::endl
    << "phi    = " << electron.phi() << std::endl
    << "energy = " << electron.energy() << std::endl
    << "caloIso                        = " << electron.caloIso() << std::endl
    << "ecalIso                        = " << electron.ecalIso() << std::endl
    << "hcalIso                        = " << electron.hcalIso() << std::endl
    << "puppiChargedHadronIso          = " << electron.puppiChargedHadronIso() << std::endl
    << "puppiNeutralHadronIso          = " << electron.puppiNeutralHadronIso() << std::endl
    << "puppiPhotonIso                 = " << electron.puppiPhotonIso() << std::endl
    << "puppiNoLeptonsChargedHadronIso = " << electron.puppiNoLeptonsChargedHadronIso() << std::endl
    << "puppiNoLeptonsNeutralHadronIso = " << electron.puppiNoLeptonsNeutralHadronIso() << std::endl
    << "puppiNoLeptonsPhotonIso        = " << electron.puppiNoLeptonsPhotonIso() << std::endl
    << "trackIso                       = " << electron.trackIso() << std::endl;
    if (PrivateMC) {
        std::cout << "cutBasedElectronID-Spring15-25ns-V1-standalone-loose  = " << electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose") << std::endl;
        std::cout << "cutBasedElectronID-Spring15-25ns-V1-standalone-medium = " << electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium") << std::endl;
        std::cout << "cutBasedElectronID-Spring15-25ns-V1-standalone-tight  = " << electron.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight") << std::endl;
        std::cout << "mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose = " << electron.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose") << std::endl;
        std::cout << "mvaEleID-Spring15-25ns-nonTrig-V1-wp80    = " << electron.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wp80") << std::endl;
        std::cout << "mvaEleID-Spring15-25ns-nonTrig-V1-wp90    = " << electron.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wp90") << std::endl;
    } else {
        std::cout << "cutBasedElectronID-Fall17-94X-V1-loose  = " << electron.electronID("cutBasedElectronID-Fall17-94X-V1-loose") << std::endl
        << "cutBasedElectronID-Fall17-94X-V1-medium = " << electron.electronID("cutBasedElectronID-Fall17-94X-V1-medium") << std::endl
        << "cutBasedElectronID-Fall17-94X-V1-tight  = " << electron.electronID("cutBasedElectronID-Fall17-94X-V1-tight") << std::endl
        << "cutBasedElectronID-Fall17-94X-V1-veto   = " << electron.electronID("cutBasedElectronID-Fall17-94X-V1-veto") << std::endl
        << "mvaEleID-Fall17-iso-V1-wp80      = " << electron.electronID("mvaEleID-Fall17-iso-V1-wp80") << std::endl
        << "mvaEleID-Fall17-iso-V1-wp90      = " << electron.electronID("mvaEleID-Fall17-iso-V1-wp90") << std::endl
        << "mvaEleID-Fall17-iso-V1-wpLoose   = " << electron.electronID("mvaEleID-Fall17-iso-V1-wpLoose") << std::endl
        //<< "mvaEleID-Fall17-iso-V2-wp80      = " << electron.electronID("mvaEleID-Fall17-iso-V2-wp80") << std::endl
        //<< "mvaEleID-Fall17-iso-V2-wp90      = " << electron.electronID("mvaEleID-Fall17-iso-V2-wp90") << std::endl
        //<< "mvaEleID-Fall17-iso-V2-wpLoose   = " << electron.electronID("mvaEleID-Fall17-iso-V2-wpLoose") << std::endl
        << "mvaEleID-Fall17-noIso-V1-wp80    = " << electron.electronID("mvaEleID-Fall17-noIso-V1-wp80") << std::endl
        << "mvaEleID-Fall17-noIso-V1-wp90    = " << electron.electronID("mvaEleID-Fall17-noIso-V1-wp90") << std::endl
        << "mvaEleID-Fall17-noIso-V1-wpLoose = " << electron.electronID("mvaEleID-Fall17-noIso-V1-wpLoose") << std::endl;
        //<< "mvaEleID-Fall17-noIso-V2-wp80    = " << electron.electronID("mvaEleID-Fall17-noIso-V2-wp80") << std::endl
        //<< "mvaEleID-Fall17-noIso-V2-wp90    = " << electron.electronID("mvaEleID-Fall17-noIso-V2-wp90") << std::endl;
        //<< "mvaEleID-Fall17-noIso-V2-wpLoose = " << electron.electronID("mvaEleID-Fall17-noIso-V2-wpLoose") << std::endl;
    }
    std::cout << "ecalTrkEnergyPostCorr    = " << electron.userFloat("ecalTrkEnergyPostCorr") << std::endl
    << "ecalTrkEnergyErrPostCorr = " << electron.userFloat("ecalTrkEnergyErrPostCorr") << std::endl;
}

void printMuonCandidate(pat::Muon &muon) {
    std::cout << "Muon parameters" << std::endl
    << "pt     = " << muon.pt() << std::endl
    << "eta    = " << muon.eta() << std::endl
    << "phi    = " << muon.phi() << std::endl
    << "energy = " << muon.energy() << std::endl
    << "caloIso                        = " << muon.caloIso() << std::endl
    << "ecalIso                        = " << muon.ecalIso() << std::endl
    << "hcalIso                        = " << muon.hcalIso() << std::endl
    << "puppiChargedHadronIso          = " << muon.puppiChargedHadronIso() << std::endl
    << "puppiNeutralHadronIso          = " << muon.puppiNeutralHadronIso() << std::endl
    << "puppiPhotonIso                 = " << muon.puppiPhotonIso() << std::endl
    << "puppiNoLeptonsChargedHadronIso = " << muon.puppiNoLeptonsChargedHadronIso() << std::endl
    << "puppiNoLeptonsNeutralHadronIso = " << muon.puppiNoLeptonsNeutralHadronIso() << std::endl
    << "puppiNoLeptonsPhotonIso        = " << muon.puppiNoLeptonsPhotonIso() << std::endl
    << "trackIso                       = " << muon.trackIso() << std::endl
    << "CutBasedIdLoose            = " << muon.passed(reco::Muon::CutBasedIdLoose) << std::endl
    << "CutBasedIdMedium           = " << muon.passed(reco::Muon::CutBasedIdMedium) << std::endl
    << "CutBasedIdTight            = " << muon.passed(reco::Muon::CutBasedIdTight) << std::endl
    << "CutBasedIdGlobalHighPt     = " << muon.passed(reco::Muon::CutBasedIdGlobalHighPt) << std::endl
    << "CutBasedIdTrkHighPt        = " << muon.passed(reco::Muon::CutBasedIdTrkHighPt) << std::endl
    << "PFIsoLoose                 = " << muon.passed(reco::Muon::PFIsoLoose) << std::endl
    << "PFIsoMedium                = " << muon.passed(reco::Muon::PFIsoMedium) << std::endl
    << "PFIsoTight                 = " << muon.passed(reco::Muon::PFIsoTight) << std::endl
    << "TkIsoLoose                 = " << muon.passed(reco::Muon::TkIsoLoose) << std::endl
    << "TkIsoTight                 = " << muon.passed(reco::Muon::TkIsoTight) << std::endl
    << "SoftCutBasedId             = " << muon.passed(reco::Muon::SoftCutBasedId) << std::endl
    << "SoftMvaId                  = " << muon.passed(reco::Muon::SoftMvaId) << std::endl
    << "MiniIsoLoose               = " << muon.passed(reco::Muon::MiniIsoLoose) << std::endl
    << "MiniIsoMedium              = " << muon.passed(reco::Muon::MiniIsoMedium) << std::endl
    << "MiniIsoTight               = " << muon.passed(reco::Muon::MiniIsoTight) << std::endl
    //<< "MultiIsoLoose              = " << muon.passed(reco::Muon::MultiIsoLoose) << std::endl
    //<< "PuppiIsoLoose              = " << muon.passed(reco::Muon::PuppiIsoLoose) << std::endl
    << "MvaLoose                   = " << muon.passed(reco::Muon::MvaLoose) << std::endl
    << "MvaMedium                  = " << muon.passed(reco::Muon::MvaMedium) << std::endl
    << "MvaTight                   = " << muon.passed(reco::Muon::MvaTight) << std::endl;
    //<< "MvaVTight                  = " << muon.passed(reco::Muon::MvaVTight) << std::endl
    //<< "MvaVVTight                 = " << muon.passed(reco::Muon::MvaVVTight) << std::endl;
}