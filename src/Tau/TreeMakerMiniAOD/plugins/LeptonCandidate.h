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

    //--- constructor & destructor
    LeptonCandidate();
    LeptonCandidate(const pat::Muon &muon, math::XYZPoint pv_position) 
        : Pt(muon.pt()), Eta(muon.eta()), Phi(muon.phi()), Flavor(13), Charge(muon.charge()), caloIso(muon.caloIso()), ecalIso(muon.ecalIso()), hcalIso(muon.hcalIso()),
        puppiChargedHadronIso(muon.puppiChargedHadronIso()), puppiNeutralHadronIso(muon.puppiNeutralHadronIso()), puppiPhotonIso(muon.puppiPhotonIso()),
        puppiNoLeptonsChargedHadronIso(muon.puppiNoLeptonsChargedHadronIso()), puppiNoLeptonsNeutralHadronIso(muon.puppiNoLeptonsNeutralHadronIso()), puppiNoLeptonsPhotonIso(muon.puppiNoLeptonsPhotonIso()),
        trackIso(muon.trackIso()), FourMomentum(muon.p4()), Dz((pv_position - muon.vertex()).R())
    {
    }
    LeptonCandidate(const pat::Electron &electron, math::XYZPoint pv_position) 
        : Pt(electron.pt()), Eta(electron.eta()), Phi(electron.phi()), Flavor(11), Charge(electron.charge()), caloIso(electron.caloIso()), ecalIso(electron.ecalIso()), hcalIso(electron.hcalIso()),
        puppiChargedHadronIso(electron.puppiChargedHadronIso()), puppiNeutralHadronIso(electron.puppiNeutralHadronIso()), puppiPhotonIso(electron.puppiPhotonIso()),
        puppiNoLeptonsChargedHadronIso(electron.puppiNoLeptonsChargedHadronIso()), puppiNoLeptonsNeutralHadronIso(electron.puppiNoLeptonsNeutralHadronIso()), puppiNoLeptonsPhotonIso(electron.puppiNoLeptonsPhotonIso()),
        trackIso(electron.trackIso()), FourMomentum(electron.p4()),Dz((pv_position - electron.vertex()).R())
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