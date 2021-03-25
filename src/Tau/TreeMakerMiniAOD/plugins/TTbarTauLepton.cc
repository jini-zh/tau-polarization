// -*- C++ -*-
//
// Package:    Tau/TTbarTauLepton
// Class:      TTbarTauLepton
//
/**\class TTbarTauLepton TTbarTauLepton.cc Tau/TTbarTauLepton/plugins/TTbarTauLepton.cc

 Description: [one line class summary]

 Implementation:
         [Notes on implementation]
*/


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

#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "PhysicsTools/Heppy/interface/TriggerBitChecker.h"

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
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

// Pi Zero libs
#include "RecoTauTag/RecoTau/interface/RecoTauPiZeroPlugins.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"
#include "RecoTauTag/RecoTau/interface/CombinatoricGenerator.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "CommonTools/CandUtils/src/AddFourMomenta.cc"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// For Btag calibration
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"

#include "Tau/TreeMakerMiniAOD/plugins/ParticleMonitor.h"
#include "Tau/TreeMakerMiniAOD/plugins/BJetCandidate.h"
#include "Tau/TreeMakerMiniAOD/plugins/LeptonCandidate.h"
#include "Tau/TreeMakerMiniAOD/plugins/GenRecoMonitor.h"
//#include "Tau/TreeMakerMiniAOD/plugins/PiZeroReconstructor.h"
//#include "Tau/TauAnalyzer/plugins/MySimpleParticle.h"

void GetLastDaughter(const reco::Candidate* &particle) {
    if (particle->numberOfDaughters() == 1) {
        particle = particle->daughter(0);
        //std::cout << "Not last particle" << std::endl;
        GetLastDaughter(particle);
    } else return;
}

void SortGenTaus(std::vector <reco::GenParticle> &items) {
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

std::string decimal_to_binary_string(UInt_t n) {
    std::string result = "";
    do {
        result = (std::to_string(n % 2)) + result;
        n = n / 2;
    } while (n > 0);
    return result;
}

//
// class declaration
//

class TTbarTauLepton : public edm::EDAnalyzer {
public:
    explicit TTbarTauLepton(const edm::ParameterSet&);
    ~TTbarTauLepton();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    
    bool AddTau            (const edm::Event&);
    void FindGenTau        (const edm::Event&);
    bool AddLepton         (const edm::Event&);
    bool AddMET            (const edm::Event&);
    bool AddVertex         (const edm::Event&);
    void CountTracks       (const edm::Event&);
    bool JetPtSum          (const edm::Event&);
    void AddWT             (const edm::Event&);
    void AddWTData         (const edm::Event&);
    void GenTauDM          (const edm::Event&);
    bool TriggerOK         (const edm::Event&);

    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    
    // ----------member data ---------------------------

    edm::EDGetTokenT<pat::TauCollection> TauCollectionToken_;
    edm::EDGetTokenT<pat::MuonCollection> MuonCollectionToken_;
    edm::EDGetTokenT<pat::ElectronCollection> ElectronCollectionToken_;
    edm::EDGetTokenT<pat::JetCollection> PuppiJetCollectionToken_;
    edm::EDGetTokenT<pat::METCollection> MetCollectionToken_;
    edm::EDGetTokenT<pat::METCollection> PuppiMetCollectionToken_;
    edm::EDGetTokenT<reco::VertexCollection> PVToken_;
    // SV coolection
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> SVToken_;
    //edm::EDGetTokenT<pat::PackedCandidateCollection> GenParticleToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> GenParticleToken_;
    // PF candidates collection
    edm::EDGetTokenT<pat::PackedCandidateCollection> PackedCandidateCollectionToken_;
    
    edm::EDGetTokenT<pat::IsolatedTrackCollection> TrackToken_;
    edm::EDGetTokenT<GenEventInfoProduct> GenEventInfoToken_;
    edm::EDGetTokenT<GenRunInfoProduct> GenRunInfoToken_;
    //edm::EDGetTokenT<reco::GenMETCollection> tok_GenMetTrue_;
    // HLT
    edm::InputTag theTriggerResultsLabel;
    edm::EDGetTokenT<edm::TriggerResults> tok_trigRes;
    std::vector<std::string>  trigNamesTarget1;
    std::vector<std::string>  trigNamesTarget2;
    std::vector<std::string>  trigNamesTarget3;
    std::vector<std::string>  trigNames1;
    std::vector<std::string>  trigNames2;
    std::vector<std::string>  trigNames3;
    std::vector<std::string>  trigNames4;
    std::vector<std::string>  trigNames5;
    // trigger prescales
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> tok_triggerObjects;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> tok_triggerPrescales;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> tok_triggerPrescalesL1min;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> tok_triggerPrescalesL1max;
    // To obtain trigger prescale
    int resultTriggerWeight;
    int triggerPrescaleHLT;
    int triggerPrescaleL1max;
    int triggerPrescaleL1min;

    // new method of trigger prescales exatraction
    //edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    //edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    //edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

    edm::EDGetTokenT<double> TauSpinnerWTToken_;
    edm::EDGetTokenT<double> TauSpinnerWTFlipToken_;
    edm::EDGetTokenT<double> TauSpinnerWThminusToken_;
    edm::EDGetTokenT<double> TauSpinnerWThplusToken_;
    edm::EDGetTokenT<bool>   TauSpinnerWTisValidToken_;
    edm::EDGetTokenT<int>    TauSpinnerMotherToken_;
    
    
    TTree * tree;
    TH1F * allTauPt;
    
    int t_Run;
    int t_Event;

    int nVtx;
    int nTrks;
    int nTau;
    int nTauC;

    // taus
    
    double tau_pt;
    double tau_eta;
    double tau_phi;
    double tau_dm;
    double tau_q;
    double tau_m;
    double tau_dz;

    // tau SV
    int tau_hasSV;
    double tau_SVdR;
    double pv_SVdR;
    double SV_Chi2;
    double SV_Chi2NDF;
    double tauPVtoSVdPhi;
    double tau_dR2; 
    
    // MVA raw
    double tau_absIso;
    double tau_againstElectronRaw;
    double tau_IsoMVArun2v1DBdR03oldDMwLTraw;
    double tau_IsoMVArun2v1DBnewDMwLTraw;
    double tau_IsoMVArun2v1DBoldDMwLTraw;
    double tau_IsoMVArun2v1PWdR03oldDMwLTraw;
    double tau_IsoMVArun2v1PWnewDMwLTraw;
    double tau_IsoMVArun2v1PWoldDMwLTraw;

    // Deep raw
    double tau_Deep2017v2ElectronRejection;
    double tau_Deep2017v2MuonRejection;
    double tau_Deep2017v2JetRejection;

    // Deep WP
    int tau_VVLooseDeepTau2017v2VSjet;
    int tau_VLooseDeepTau2017v2VSjet;
    int tau_LooseDeepTau2017v2VSjet;
    int tau_MediumDeepTau2017v2VSjet;
    int tau_TightDeepTau2017v2VSjet;
    int tau_VTightDeepTau2017v2VSjet;
    int tau_VVTightDeepTau2017v2VSjet;

    int tau_LooseDeepTau2017v2VSmu;
    int tau_MediumDeepTau2017v2VSmu;
    int tau_TightDeepTau2017v2VSmu;

    int tau_VVLooseDeepTau2017v2VSe;
    int tau_VLooseDeepTau2017v2VSe;
    int tau_LooseDeepTau2017v2VSe;
    int tau_MediumDeepTau2017v2VSe;
    int tau_TightDeepTau2017v2VSe;
    int tau_VTightDeepTau2017v2VSe;
    int tau_VVTightDeepTau2017v2VSe;

    // MVA WP
    int tau_looseCombinedIso;
    int tau_mediumCombinedIso;
    int tau_tightCombinedIso;
    int tau_VlooseMvaIso;
    int tau_looseMvaIso;
    int tau_mediumMvaIso;
    int tau_tightMvaIso;
    int tau_VtightMvaIso;
    int tau_VVtightMvaIso;
    int tau_tightMuonRejection;
    int tau_looseElectronRejection;
    int tau_mediumElectronRejection;
    int tau_tightElectronRejection;
    int tau_VtightElectronRejection;
    int decayModeFindingNewDMs;
    int decayModeFinding;

    // Tau decay products
    double piChar_pt;
    double piChar_eta;
    double piChar_phi;
    double piChar_q;
    double piChar_m;

    double piZero_pt;
    double piZero_eta;
    double piZero_phi;
    double piZero_m;

    double DiPhoton_pt;
    double DiPhoton_eta;
    double DiPhoton_phi;
    double DiPhoton_m;

    // Gamma candidates
    double TauPhoton1_pt;
    double TauPhoton1_eta;
    double TauPhoton1_phi;
    double TauPhoton1_energy;

    double TauPhoton2_pt;
    double TauPhoton2_eta;
    double TauPhoton2_phi;
    double TauPhoton2_energy;
    
    double pipiMass;
    double upsilon;

    double tau_found;
    double gentau_found;
    int    gentau_dm;
    double dR;
    double bquark1dR;
    double bquark2dR;
    double leptondR;
    double genTauFromW;
    double genTauFromWFromt;
    int genTauMother;
    int genLeptonMother;
    double W_pt, W_eta, W_phi, W_energy, W_charge;
    double genLeptonFromW;
    double genLeptonFromWFromt;
    int genb1Fromt;
    int genb2Fromt;
    int genb1Mother;
    int genb2Mother;
    int genb1_flavor;
    double genb1_pt, genb1_eta, genb1_phi, genb1_energy;
    int genb2_flavor;
    double genb2_pt, genb2_eta, genb2_phi, genb2_energy;
    int genTauMisID;

    // Generated tau parameters
    double gentau_pt, gentau_energy, gentau_eta, gentau_phi;
    double genPiChar_pt, genPiChar_energy, genPiChar_eta, genPiChar_phi;
    double genPi0_pt, genPi0_energy, genPi0_eta, genPi0_phi;
    double nuW_pt, nuW_energy, nuW_eta, nuW_phi;
    double nutau_pt, nutau_energy, nutau_eta, nutau_phi;
    int nNu;
    double SumNu_pt, SumNu_eta, SumNu_phi, SumNu_energy;
    double nunu_pt;
    double gentau_vis_pt, gentau_vis_eta, gentau_vis_phi, gentau_vis_energy;
    double genLepton_eta, genLepton_phi, genLepton_pt, genLepton_energy;
    int genLepton_mother, genLepton_flavor;

    double PuppijetPtSum20, PuppijetPtSum30, PuppijetPtSum20PV, PuppijetPtSum30PV;
    int nPuppiJets20, nPuppiJets20PV, nPuppiJets30, nPuppiJets30PV;
    int nLooseBtagedPuppiJets, nMediumBtagedPuppiJets, nTightBtagedPuppiJets;
    int nLooseBtagedPuppiJetsPV, nMediumBtagedPuppiJetsPV, nTightBtagedPuppiJetsPV;
    std::vector <BJetCandidate> looseBJetsCopy;
    // Jet separately
    double Jet1_pt;
    double Jet1_eta;
    double Jet1_phi;
    double Jet1_m;
    double Jet1_E;
    double Jet1_bprob, Jet1_bbprob;
    int Jet1_hadronFlavour;
    //
    double Jet2_pt;
    double Jet2_eta;
    double Jet2_phi;
    double Jet2_m;
    double Jet2_E;
    double Jet2_bprob, Jet2_bbprob;
    int Jet2_hadronFlavour;
    //
    double Jet3_pt;
    double Jet3_eta;
    double Jet3_phi;
    double Jet3_m;
    double Jet3_E;
    double Jet3_bprob, Jet3_bbprob;
    int Jet3_hadronFlavour;
    //
    double Jet4_pt;
    double Jet4_eta;
    double Jet4_phi;
    double Jet4_m;
    double Jet4_E;
    double Jet4_bprob, Jet4_bbprob;
    int Jet4_hadronFlavour;
    //
    double Jet5_pt;
    double Jet5_eta;
    double Jet5_phi;
    double Jet5_m;
    double Jet5_E;
    double Jet5_bprob, Jet5_bbprob;
    int Jet5_hadronFlavour;
    //
    double Jet6_pt;
    double Jet6_eta;
    double Jet6_phi;
    double Jet6_m;
    double Jet6_E;
    double Jet6_bprob, Jet6_bbprob;
    int Jet6_hadronFlavour;

    // Tau jet
    double TauJet_pt;
    double TauJet_eta;
    double TauJet_phi;
    double TauJet_E;
    double TauJet_bprob;
    int TauJet_hadronFlavour;

    // Leptons
    double lepton1_pt;
    double lepton1_eta;
    double lepton1_phi;
    double lepton1_E;
    double lepton1_dz;
    int    lepton1_flavor;
    int    lepton1_charge;
    float  lepton1_trackIso;
    float  lepton1_sumPuppiIso;
    float  lepton1_sumPuppiNoLeptonIso;
    // Jet in cone dR = 0.4
    double lepton1Jet_pt;
    double lepton1Jet_eta;
    double lepton1Jet_phi;
    double lepton1Jet_E;
    double lepton1Jet_bprob;
    int    lepton1_electron_cutBasedID_loose;
    int    lepton1_electron_cutBasedID_medium;
    int    lepton1_electron_cutBasedID_tight;
    int    lepton1_electron_mvaNoIsoID_loose;
    int    lepton1_electron_mvaNoIsoID_wp80;
    int    lepton1_electron_mvaNoIsoID_wp90;
    int    lepton1_muon_CutBasedIdLoose;
    int    lepton1_muon_CutBasedIdMedium;
    int    lepton1_muon_CutBasedIdTight;
    int    lepton1_muon_CutBasedIdGlobalHighPt;
    int    lepton1_muon_CutBasedIdTrkHighPt;
    int    lepton1_muon_PFIsoLoose;
    int    lepton1_muon_PFIsoMedium;
    int    lepton1_muon_PFIsoTight;
    int    lepton1_muon_TkIsoLoose;
    int    lepton1_muon_TkIsoTight;
    int    lepton1_muon_MvaLoose;
    int    lepton1_muon_MvaMedium;
    int    lepton1_muon_MvaTight;
    double lepton1_tau_dm;
    double lepton1_tau_m;
    double lepton1_tau_absIso;
    // Raw
    double lepton1_tau_againstElectronRaw;
    double lepton1_tau_IsoMVArun2v1DBdR03oldDMwLTraw;
    // Deep raw
    double lepton1_tau_Deep2017v2ElectronRejection;
    double lepton1_tau_Deep2017v2MuonRejection;
    double lepton1_tau_Deep2017v2JetRejection;
    // Deep WP
    int    lepton1_tau_VVLooseDeepTau2017v2VSjet;
    int    lepton1_tau_VLooseDeepTau2017v2VSjet;
    int    lepton1_tau_LooseDeepTau2017v2VSjet;
    int    lepton1_tau_MediumDeepTau2017v2VSjet;
    int    lepton1_tau_TightDeepTau2017v2VSjet;
    int    lepton1_tau_VTightDeepTau2017v2VSjet;
    int    lepton1_tau_VVTightDeepTau2017v2VSjet;
    int    lepton1_tau_LooseDeepTau2017v2VSmu;
    int    lepton1_tau_MediumDeepTau2017v2VSmu;
    int    lepton1_tau_TightDeepTau2017v2VSmu;
    int    lepton1_tau_VVLooseDeepTau2017v2VSe;
    int    lepton1_tau_VLooseDeepTau2017v2VSe;
    int    lepton1_tau_LooseDeepTau2017v2VSe;
    int    lepton1_tau_MediumDeepTau2017v2VSe;
    int    lepton1_tau_TightDeepTau2017v2VSe;
    int    lepton1_tau_VTightDeepTau2017v2VSe;
    int    lepton1_tau_VVTightDeepTau2017v2VSe;
    // MVA WP
    int    lepton1_tau_looseCombinedIso;
    int    lepton1_tau_tightCombinedIso;
    int    lepton1_tau_looseMvaIso;
    int    lepton1_tau_mediumMvaIso;
    int    lepton1_tau_tightMvaIso;
    int    lepton1_tau_VtightMvaIso;
    int    lepton1_tau_VVtightMvaIso;
    int    lepton1_tau_tightMuonRejection;
    int    lepton1_tau_looseElectronRejection;
    int    lepton1_tau_mediumElectronRejection;
    int    lepton1_tau_tightElectronRejection;
    int    lepton1_tau_VtightElectronRejection;
    // lepton 2
    double lepton2_pt;
    double lepton2_eta;
    double lepton2_phi;
    double lepton2_E;
    double lepton2_dz;
    int    lepton2_flavor;
    int    lepton2_charge;
    float  lepton2_trackIso;
    float  lepton2_sumPuppiIso;
    float  lepton2_sumPuppiNoLeptonIso;
    // Jet in cone dR = 0.4
    double lepton2Jet_pt;
    double lepton2Jet_eta;
    double lepton2Jet_phi;
    double lepton2Jet_E;
    double lepton2Jet_bprob;
    int    lepton2_electron_cutBasedID_loose;
    int    lepton2_electron_cutBasedID_medium;
    int    lepton2_electron_cutBasedID_tight;
    int    lepton2_electron_mvaNoIsoID_loose;
    int    lepton2_electron_mvaNoIsoID_wp80;
    int    lepton2_electron_mvaNoIsoID_wp90;
    int    lepton2_muon_CutBasedIdLoose;
    int    lepton2_muon_CutBasedIdMedium;
    int    lepton2_muon_CutBasedIdTight;
    int    lepton2_muon_CutBasedIdGlobalHighPt;
    int    lepton2_muon_CutBasedIdTrkHighPt;
    int    lepton2_muon_PFIsoLoose;
    int    lepton2_muon_PFIsoMedium;
    int    lepton2_muon_PFIsoTight;
    int    lepton2_muon_TkIsoLoose;
    int    lepton2_muon_TkIsoTight;
    int    lepton2_muon_MvaLoose;
    int    lepton2_muon_MvaMedium;
    int    lepton2_muon_MvaTight;
    double lepton2_tau_dm;
    double lepton2_tau_m;
    double lepton2_tau_absIso;
    // Raw
    double lepton2_tau_againstElectronRaw;
    double lepton2_tau_IsoMVArun2v1DBdR03oldDMwLTraw;
    // Deep raw
    double lepton2_tau_Deep2017v2ElectronRejection;
    double lepton2_tau_Deep2017v2MuonRejection;
    double lepton2_tau_Deep2017v2JetRejection;
    // Deep WP
    int    lepton2_tau_VVLooseDeepTau2017v2VSjet;
    int    lepton2_tau_VLooseDeepTau2017v2VSjet;
    int    lepton2_tau_LooseDeepTau2017v2VSjet;
    int    lepton2_tau_MediumDeepTau2017v2VSjet;
    int    lepton2_tau_TightDeepTau2017v2VSjet;
    int    lepton2_tau_VTightDeepTau2017v2VSjet;
    int    lepton2_tau_VVTightDeepTau2017v2VSjet;
    int    lepton2_tau_LooseDeepTau2017v2VSmu;
    int    lepton2_tau_MediumDeepTau2017v2VSmu;
    int    lepton2_tau_TightDeepTau2017v2VSmu;
    int    lepton2_tau_VVLooseDeepTau2017v2VSe;
    int    lepton2_tau_VLooseDeepTau2017v2VSe;
    int    lepton2_tau_LooseDeepTau2017v2VSe;
    int    lepton2_tau_MediumDeepTau2017v2VSe;
    int    lepton2_tau_TightDeepTau2017v2VSe;
    int    lepton2_tau_VTightDeepTau2017v2VSe;
    int    lepton2_tau_VVTightDeepTau2017v2VSe;
    // MVA WP
    int    lepton2_tau_looseCombinedIso;
    int    lepton2_tau_tightCombinedIso;
    int    lepton2_tau_looseMvaIso;
    int    lepton2_tau_mediumMvaIso;
    int    lepton2_tau_tightMvaIso;
    int    lepton2_tau_VtightMvaIso;
    int    lepton2_tau_VVtightMvaIso;
    int    lepton2_tau_tightMuonRejection;
    int    lepton2_tau_looseElectronRejection;
    int    lepton2_tau_mediumElectronRejection;
    int    lepton2_tau_tightElectronRejection;
    int    lepton2_tau_VtightElectronRejection;
    //
    int    nLeptonCandidates;
    int    VetoLeptons;
    int    VetoTaus;

    int nPhotons;
    int nGamma;

    // MET
    double met;
    double met_phi;
    double met_eta;
    double met_energy;
    double tauMET_mass;
    double met_significance;
    double met_mEtSig;
    
    double m_t;
    double dPhi;

    // Puppi MET
    double Puppimet;
    double Puppimet_phi;
    double Puppimet_eta;
    double Puppimet_energy;
    double tauPuppimet_mass;
    double Puppimet_significance;
    double Puppimet_metSig;
    double dPhiPuppimetTau;

    // Trigger matching
    int nTargetTrigger1;
    int nTargetTrigger2;
    int nTargetTrigger3;
    UInt_t TriggerBit1; // 32 bit
    UInt_t TriggerBit2;
    UInt_t TriggerBit3;
    UInt_t TriggerBit4;
    UInt_t TriggerBit5;
    
    math::XYZPoint pv_position;
    math::XYZPoint SV_position;

    double WT;
    double WTFlip;
    double WThminus;
    double WThplus;
    bool WTisValid;
    int  TauSpinnerMother;

    // Type of event
    int ttbarEvent;

    //////////////////////////////////////////////////////
    bool isMC;
    bool TauSpinnerOn;
    double tauPtMin;
    double piPtMin;
    double tauEtaMax;
    double tauDzMax;
    double BJetPtMin;
    double JetEtaMax;
    double MuElePtMin;
    double EtaMax;
    double null;
    bool monitoring;
    bool monitoringHLT;
    bool monitoringTau;
    bool monitoringGen;
    bool monitoringJets;
    bool monitoringBJets;
    bool monitoringLeptons;
    bool monitoringMET;
    bool looseTauID;
    bool DeepTau;
    bool UseTau;
    bool LeptonRequired; 

    int iT;
    int nTrigger = 0;
    int nTau1 = 0;
    int nTauJet = 0;
    int nJet = 0;
    int nLepton = 0;
    int nEvent = 0;
    int nPassed = 0;
    int nTauJetInEvent;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TTbarTauLepton::TTbarTauLepton(const edm::ParameterSet& iConfig) {
    //now do what ever initialization is needed
    iT =0;

    isMC                        = iConfig.getParameter<bool>("isMC");
    monitoring                  = iConfig.getParameter<bool>("monitoring");
    monitoringHLT               = iConfig.getParameter<bool>("monitoringHLT");
    monitoringTau               = iConfig.getParameter<bool>("monitoringTau");
    monitoringGen               = iConfig.getParameter<bool>("monitoringGen");
    monitoringJets              = iConfig.getParameter<bool>("monitoringJets");
    monitoringBJets             = iConfig.getParameter<bool>("monitoringBJets");
    monitoringLeptons           = iConfig.getParameter<bool>("monitoringLeptons");
    monitoringMET               = iConfig.getParameter<bool>("monitoringMET");
    TauSpinnerOn                = iConfig.getParameter<bool>("TauSpinnerOn");
    tauPtMin                    = iConfig.getParameter<double>("tauPtMin");
    piPtMin                     = iConfig.getParameter<double>("piPtMin");
    tauEtaMax                   = iConfig.getParameter<double>("tauEtaMax");
    tauDzMax                    = iConfig.getParameter<double>("tauDzMax");
    looseTauID                  = iConfig.getParameter<bool>("looseTauID");
    BJetPtMin                   = iConfig.getParameter<double>("BJetPtMin"); // 30
    JetEtaMax                   = iConfig.getParameter<double>("JetEtaMax");
    MuElePtMin                  = iConfig.getParameter<double>("MuElePtMin"); // 20
    EtaMax                      = iConfig.getParameter<double>("EtaMax"); // 2.4
    DeepTau                     = iConfig.getParameter<bool>("DeepTau");
    UseTau                      = iConfig.getParameter<bool>("UseTau");
    LeptonRequired              = iConfig.getParameter<bool>("LeptonRequired");
    null                        = iConfig.getParameter<double>("null");

    std::string tauCollection         = iConfig.getParameter<std::string>("tauCollection");
    std::string muonCollection        = iConfig.getParameter<std::string>("muonCollection");
    std::string electronCollection    = iConfig.getParameter<std::string>("electronCollection");
    std::string PuppijetCollection    = iConfig.getParameter<std::string>("PuppijetCollection");
    std::string metCollection         = iConfig.getParameter<std::string>("metCollection");
    std::string PuppimetCollection    = iConfig.getParameter<std::string>("PuppimetCollection");
    std::string vertexCollection      = iConfig.getParameter<std::string>("vertexCollection");
    std::string SVCollection          = iConfig.getParameter<std::string>("SVCollection");
    std::string genParticleCollection = iConfig.getParameter<std::string>("genParticleCollection");
    std::string trackCollection       = iConfig.getParameter<std::string>("trackCollection");
    theTriggerResultsLabel            = edm::InputTag("TriggerResults","","HLT");
    trigNamesTarget1                  = iConfig.getParameter<std::vector<std::string>>("TriggerTarget1");
    trigNamesTarget2                  = iConfig.getParameter<std::vector<std::string>>("TriggerTarget2");
    trigNamesTarget3                  = iConfig.getParameter<std::vector<std::string>>("TriggerTarget3");
    trigNames1                        = iConfig.getParameter<std::vector<std::string>>("Triggers1");
    trigNames2                        = iConfig.getParameter<std::vector<std::string>>("Triggers2");
    trigNames3                        = iConfig.getParameter<std::vector<std::string>>("Triggers3");
    trigNames4                        = iConfig.getParameter<std::vector<std::string>>("Triggers4");
    trigNames5                        = iConfig.getParameter<std::vector<std::string>>("Triggers5");
    std::string PackedCandidateCollection = iConfig.getParameter<std::string>("PackedCandidateCollection");
    
    TauCollectionToken_         = consumes<pat::TauCollection>(edm::InputTag(tauCollection));
    MuonCollectionToken_        = consumes<pat::MuonCollection>(edm::InputTag(muonCollection));
    ElectronCollectionToken_    = consumes<pat::ElectronCollection>(edm::InputTag(electronCollection));
    PuppiJetCollectionToken_    = consumes<pat::JetCollection>(edm::InputTag(PuppijetCollection));
    MetCollectionToken_         = consumes<pat::METCollection>(edm::InputTag(metCollection));
    PuppiMetCollectionToken_    = consumes<pat::METCollection>(edm::InputTag(PuppimetCollection));
    PVToken_                    = consumes<reco::VertexCollection>(edm::InputTag(vertexCollection));
    SVToken_                    = consumes<reco::VertexCompositePtrCandidateCollection>(edm::InputTag(SVCollection));
    GenParticleToken_           = consumes<reco::GenParticleCollection>(edm::InputTag(genParticleCollection));
    TrackToken_                 = consumes<pat::IsolatedTrackCollection>(edm::InputTag(trackCollection));
    tok_trigRes                 = consumes<edm::TriggerResults>(theTriggerResultsLabel);
    PackedCandidateCollectionToken_ = consumes<pat::PackedCandidateCollection>(edm::InputTag(PackedCandidateCollection));
    tok_triggerObjects          = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("Triggerobjects"));
    tok_triggerPrescales        = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")); 
    tok_triggerPrescalesL1min   = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescalesL1min")); 
    tok_triggerPrescalesL1max   = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescalesL1max"));
    //tok_GenMetTrue_                 = consumes<reco::GenMETCollection>( iConfig.getParameter<edm::InputTag>("genMetTrue"));

    if (isMC && TauSpinnerOn) {
        TauSpinnerWTToken_        = consumes<double>(iConfig.getParameter<edm::InputTag>("WTCollection"));
        TauSpinnerWTFlipToken_    = consumes<double>(iConfig.getParameter<edm::InputTag>("WTFlipCollection"));
        TauSpinnerWThminusToken_  = consumes<double>(iConfig.getParameter<edm::InputTag>("WThminusCollection"));
        TauSpinnerWThplusToken_   = consumes<double>(iConfig.getParameter<edm::InputTag>("WThplusCollection"));
        TauSpinnerWTisValidToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("WTisValidCollection"));
        TauSpinnerMotherToken_    = consumes<int>(iConfig.getParameter<edm::InputTag>("MotherCollection"));
    }

}


TTbarTauLepton::~TTbarTauLepton() {
     // do anything here that needs to be done at desctruction time
     // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void TTbarTauLepton::analyze(const edm::Event& event, const edm::EventSetup&) {
    t_Run   = event.id().run();
    t_Event = event.id().event();
    nEvent++;
    nTauJetInEvent = 0;

    // Operation with triggers (filtering, scaling)
    if (!TriggerOK(event)) {
        nTrigger++;
        if (monitoring) std::cout << "Trigger" << std::endl;
        return;
    }
    // Primary vertex position
    if (!AddVertex(event)) {
        if (monitoring) std::cout << "Vertex" << std::endl;
        return;
    }
    // MET parameters
    if (!AddMET(event)) {
        if (monitoring) std::cout << "MET" << std::endl;
        return;
    }
    // Find the most isolated Tau
    bool tauFound;
    if (!AddTau(event)) {
        nTau1++;
        if (monitoring) std::cout << "Tau" << std::endl;
        return;
    }
    //AddPackedCandidates(event);
    //
    CountTracks(event);
    // Find second lepton for ttbar
    if (!AddLepton(event)) {
        nLepton++;
        if (monitoring) std::cout << "Lepton" << std::endl;
        return;
    }
    // Summary of jets parameters
    if (!JetPtSum(event)) {
        nJet++;
        if (monitoring) std::cout << "Jet" << std::endl;
        return;
    }
    // Generated particles parameters
    FindGenTau(event);
    // Decay mode of generated tau
    //GenTauDM(event);
    // TauSpinner
    if (isMC && TauSpinnerOn) {
        AddWT(event);
    } else {
        AddWTData(event);
    }

    nPassed++;

    if (monitoring) {
        std::cout << "Events  = " << nEvent << std::endl;
        std::cout << "Trigger = " << nTrigger << std::endl;
        std::cout << "Tau     = " << nTau1 << std::endl;
        std::cout << "Lepton  = " << nLepton << std::endl;
        std::cout << "Jet     = " << nJet << std::endl;
        std::cout << "TauJets = " << nTauJetInEvent << std::endl;
        std::cout << "Passed  = " << nPassed << std::endl;
    }

    TLorentzVector pi0, pi1;
    pi0.SetPtEtaPhiM(piZero_pt, piZero_eta, piZero_phi, piZero_m);
    pi1.SetPtEtaPhiM(piChar_pt, piChar_eta, piChar_phi, piChar_m);

    pipiMass = (pi0 + pi1).M();
    //upsilon  = (pi1.E() - pi0.E()) / tau_pt;
    upsilon  = 2 * piChar_pt / tau_pt - 1;
    dPhi            = dphi(tau_phi, met_phi);
    dPhiPuppimetTau = dphi(tau_phi, Puppimet_phi);
    m_t             = sqrt(2 * tau_pt * met * (1 - cos(dPhi)));
        
    tree->Fill();
};


// ------------ method called once each job just before starting event loop  ------------
void TTbarTauLepton::beginJob() {
    // declaring the tree and its branches.
    edm::Service<TFileService> FS;
    tree = FS->make<TTree>("tree", "tree", 1);
    
    tree->Branch("t_Run",  &t_Run,  "t_Run/I");
    tree->Branch("t_Event",&t_Event,"t_Event/I");

    // Tau kinematic parameters
    tree->Branch("tau_pt",&tau_pt,"tau_pt/D");
    tree->Branch("tau_eta",&tau_eta,"tau_eta/D");
    tree->Branch("tau_phi",&tau_phi,"tau_phi/D");
    tree->Branch("tau_q",&tau_q,"tau_q/D");
    tree->Branch("tau_m",&tau_m,"tau_m/D");
    tree->Branch("tau_dm",&tau_dm,"tau_dm/D");
    tree->Branch("decayModeFindingNewDMs",&decayModeFindingNewDMs,"decayModeFindingNewDMs/D");
    tree->Branch("decayModeFinding",&decayModeFinding,"decayModeFinding/D");
    tree->Branch("tau_dz",&tau_dz,"tau_dz/D");

    // Raw discriminators
    tree->Branch("tau_absIso",&tau_absIso,"tau_absIso/D");
    tree->Branch("tau_againstElectronRaw",&tau_againstElectronRaw,"tau_againstElectronRaw/D");
    tree->Branch("tau_IsoMVArun2v1DBdR03oldDMwLTraw",&tau_IsoMVArun2v1DBdR03oldDMwLTraw,"tau_IsoMVArun2v1DBdR03oldDMwLTraw/D");
    tree->Branch("tau_IsoMVArun2v1DBnewDMwLTraw",&tau_IsoMVArun2v1DBnewDMwLTraw,"tau_IsoMVArun2v1DBnewDMwLTraw/D");
    tree->Branch("tau_IsoMVArun2v1DBoldDMwLTraw",&tau_IsoMVArun2v1DBoldDMwLTraw,"tau_IsoMVArun2v1DBoldDMwLTraw/D");
    tree->Branch("tau_IsoMVArun2v1PWdR03oldDMwLTraw",&tau_IsoMVArun2v1PWdR03oldDMwLTraw,"tau_IsoMVArun2v1PWdR03oldDMwLTraw/D");
    tree->Branch("tau_IsoMVArun2v1PWnewDMwLTraw",&tau_IsoMVArun2v1PWnewDMwLTraw,"tau_IsoMVArun2v1PWnewDMwLTraw/D");
    tree->Branch("tau_IsoMVArun2v1PWoldDMwLTraw",&tau_IsoMVArun2v1PWoldDMwLTraw,"tau_IsoMVArun2v1PWoldDMwLTraw/D");

    // Discriminators
    tree->Branch("tau_looseCombinedIso",&tau_looseCombinedIso,"tau_looseCombinedIso/I");
    tree->Branch("tau_mediumCombinedIso",&tau_mediumCombinedIso,"tau_mediumCombinedIso/I");
    tree->Branch("tau_tightCombinedIso",&tau_tightCombinedIso,"tau_tightCombinedIso/I");
    tree->Branch("tau_VlooseMvaIso",&tau_VlooseMvaIso,"tau_VlooseMvaIso/I");
    tree->Branch("tau_looseMvaIso",&tau_looseMvaIso,"tau_looseMvaIso/I");
    tree->Branch("tau_mediumMvaIso",&tau_mediumMvaIso,"tau_mediumMvaIso/I");
    tree->Branch("tau_tightMvaIso",&tau_tightMvaIso,"tau_tightMvaIso/I");
    tree->Branch("tau_VtightMvaIso",&tau_VtightMvaIso,"tau_VtightMvaIso/I");
    tree->Branch("tau_VVtightMvaIso",&tau_VVtightMvaIso,"tau_VVtightMvaIso/I");
    tree->Branch("tau_tightMuonRejection",&tau_tightMuonRejection,"tau_tightMuonRejection/I");
    tree->Branch("tau_looseElectronRejection",&tau_looseElectronRejection,"tau_looseElectronRejection/I");
    tree->Branch("tau_mediumElectronRejection",&tau_mediumElectronRejection,"tau_mediumElectronRejection/I");
    tree->Branch("tau_tightElectronRejection",&tau_tightElectronRejection,"tau_tightElectronRejection/I");
    tree->Branch("tau_VtightElectronRejection",&tau_VtightElectronRejection,"tau_VtightElectronRejection/I");

    // Deep 2017v2
    if (DeepTau) {
        tree->Branch("tau_Deep2017v2ElectronRejection",&tau_Deep2017v2ElectronRejection,"tau_Deep2017v2ElectronRejection/D");
        tree->Branch("tau_Deep2017v2MuonRejection",&tau_Deep2017v2MuonRejection,"tau_Deep2017v2MuonRejection/D");
        tree->Branch("tau_Deep2017v2JetRejection",&tau_Deep2017v2JetRejection,"tau_Deep2017v2JetRejection/D");

        tree->Branch("tau_VVLooseDeepTau2017v2VSjet",&tau_VVLooseDeepTau2017v2VSjet,"tau_VVLooseDeepTau2017v2VSjet/I");
        tree->Branch("tau_VLooseDeepTau2017v2VSjet",&tau_VLooseDeepTau2017v2VSjet,"tau_VLooseDeepTau2017v2VSjet/I");
        tree->Branch("tau_LooseDeepTau2017v2VSjet",&tau_LooseDeepTau2017v2VSjet,"tau_LooseDeepTau2017v2VSjet/I");
        tree->Branch("tau_MediumDeepTau2017v2VSjet",&tau_MediumDeepTau2017v2VSjet,"tau_MediumDeepTau2017v2VSjet/I");
        tree->Branch("tau_TightDeepTau2017v2VSjet",&tau_TightDeepTau2017v2VSjet,"tau_TightDeepTau2017v2VSjet/I");
        tree->Branch("tau_VTightDeepTau2017v2VSjet",&tau_VTightDeepTau2017v2VSjet,"tau_VTightDeepTau2017v2VSjet/I");
        tree->Branch("tau_VVTightDeepTau2017v2VSjet",&tau_VVTightDeepTau2017v2VSjet,"tau_VVTightDeepTau2017v2VSjet/I");
        tree->Branch("tau_LooseDeepTau2017v2VSmu",&tau_LooseDeepTau2017v2VSmu,"tau_LooseDeepTau2017v2VSmu/I");
        tree->Branch("tau_MediumDeepTau2017v2VSmu",&tau_MediumDeepTau2017v2VSmu,"tau_MediumDeepTau2017v2VSmu/I");
        tree->Branch("tau_TightDeepTau2017v2VSmu",&tau_TightDeepTau2017v2VSmu,"tau_TightDeepTau2017v2VSmu/I");

        tree->Branch("tau_VVLooseDeepTau2017v2VSe",&tau_VVLooseDeepTau2017v2VSe,"tau_VVLooseDeepTau2017v2VSe/I");
        tree->Branch("tau_VLooseDeepTau2017v2VSe",&tau_VLooseDeepTau2017v2VSe,"tau_VLooseDeepTau2017v2VSe/I");
        tree->Branch("tau_LooseDeepTau2017v2VSe",&tau_LooseDeepTau2017v2VSe,"tau_LooseDeepTau2017v2VSe/I");
        tree->Branch("tau_MediumDeepTau2017v2VSe",&tau_MediumDeepTau2017v2VSe,"tau_MediumDeepTau2017v2VSe/I");
        tree->Branch("tau_TightDeepTau2017v2VSe",&tau_TightDeepTau2017v2VSe,"tau_TightDeepTau2017v2VSe/I");
        tree->Branch("tau_VTightDeepTau2017v2VSe",&tau_VTightDeepTau2017v2VSe,"tau_VTightDeepTau2017v2VSe/I");
        tree->Branch("tau_VVTightDeepTau2017v2VSe",&tau_VVTightDeepTau2017v2VSe,"tau_VVTightDeepTau2017v2VSe/I");
    }

    // Charged Pi parameters
    tree->Branch("piChar_pt", &piChar_pt, "piChar_pt/D");
    tree->Branch("piChar_eta", &piChar_eta, "piChar_eta/D");
    tree->Branch("piChar_phi", &piChar_phi, "piChar_phi/D");
    tree->Branch("piChar_q", &piChar_q, "piChar_q/D");
    tree->Branch("piChar_m", &piChar_m, "piChar_m/D");

    // Neutral Pi parameters
    tree->Branch("piZero_pt", &piZero_pt, "piZero_pt/D");
    tree->Branch("piZero_eta", &piZero_eta, "piZero_eta/D");
    tree->Branch("piZero_phi", &piZero_phi, "piZero_phi/D");
    tree->Branch("piZero_m", &piZero_m, "piZero_m/D");

    // Combinated Photons (Pi0 candidate)
    tree->Branch("DiPhoton_pt", &DiPhoton_pt, "DiPhoton_pt/D");
    tree->Branch("DiPhoton_eta", &DiPhoton_eta, "DiPhoton_eta/D");
    tree->Branch("DiPhoton_phi", &DiPhoton_phi, "DiPhoton_phi/D");
    tree->Branch("DiPhoton_m", &DiPhoton_m, "DiPhoton_m/D");

    // Gamma candidates
    tree->Branch("TauPhoton1_pt", &TauPhoton1_pt, "TauPhoton1_pt/D");
    tree->Branch("TauPhoton1_eta", &TauPhoton1_eta, "TauPhoton1_eta/D");
    tree->Branch("TauPhoton1_phi", &TauPhoton1_phi, "TauPhoton1_phi/D");
    tree->Branch("TauPhoton1_energy", &TauPhoton1_energy, "TauPhoton1_energy/D");

    tree->Branch("TauPhoton2_pt", &TauPhoton2_pt, "TauPhoton2_pt/D");
    tree->Branch("TauPhoton2_eta", &TauPhoton2_eta, "TauPhoton2_eta/D");
    tree->Branch("TauPhoton2_phi", &TauPhoton2_phi, "TauPhoton2_phi/D");
    tree->Branch("TauPhoton2_energy", &TauPhoton2_energy, "TauPhoton2_energy/D");

    // SV parameters
    tree->Branch("tau_hasSV",&tau_hasSV,"tau_hasSV/I");
    tree->Branch("tau_SVdR", &tau_SVdR, "tau_SVdR/D");
    tree->Branch("pv_SVdR", &pv_SVdR, "pv_SVdR/D");
    tree->Branch("SV_Chi2", &SV_Chi2, "SV_Chi2/D");
    tree->Branch("SV_Chi2NDF", &SV_Chi2NDF, "SV_Chi2NDF/D");
    tree->Branch("tauPVtoSVdPhi", &tauPVtoSVdPhi, "tauPVtoSVdPhi/D");
    tree->Branch("tau_dR2", &tau_dR2, "tau_dR2/D");

    tree->Branch("pipiMass", &pipiMass, "pipiMass/D");
    tree->Branch("upsilon", &upsilon, "upsilon/D");

    // Generated particles
    if (isMC) {
        tree->Branch("tau_found",&tau_found, "tau_found/D");
        tree->Branch("gentau_found",&gentau_found, "gentau_found/D");
        tree->Branch("dR", &dR, "dR/D");
        tree->Branch("leptondR", &leptondR, "leptondR/D");
        tree->Branch("bquark1dR", &bquark1dR, "bquark1dR/D");
        tree->Branch("bquark2dR", &bquark2dR, "bquark2dR/D");
        // generated W-boson (Tau mother)
        tree->Branch("genTauFromW", &genTauFromW, "genTauFromW/D");
        tree->Branch("W_pt", &W_pt, "W_pt/D");
        tree->Branch("W_eta", &W_eta, "W_eta/D");
        tree->Branch("W_phi", &W_phi, "W_phi/D");
        tree->Branch("W_energy", &W_energy, "W_energy/D");
        tree->Branch("genTauFromWFromt", &genTauFromWFromt, "genTauFromWFromt/D");
        tree->Branch("genTauMother", &genTauMother, "genTauMother/I");
        tree->Branch("genLeptonMother", &genLeptonMother, "genLeptonMother/I");
        tree->Branch("genLeptonFromW", &genLeptonFromW, "genLeptonFromW/D");
        tree->Branch("genLepton_pt", &genLepton_pt, "genLepton_pt/D");
        tree->Branch("genLepton_eta", &genLepton_eta, "genLepton_eta/D");
        tree->Branch("genLepton_phi", &genLepton_phi, "genLepton_phi/D");
        tree->Branch("genLepton_energy", &genLepton_energy, "genLepton_energy/D");
        tree->Branch("genLepton_flavor", &genLepton_flavor, "genLepton_flavor/I");
        tree->Branch("genb1Fromt", &genb1Fromt, "genb1Fromt/I");
        tree->Branch("genb2Fromt", &genb2Fromt, "genb2Fromt/I");
        tree->Branch("genb1Mother", &genb1Mother, "genb1Mother/I");
        tree->Branch("genb2Mother", &genb2Mother, "genb2Mother/I");

        tree->Branch("genb1_pt", &genb1_pt, "genb1_pt/D");
        tree->Branch("genb1_eta", &genb1_eta, "genb1_eta/D");
        tree->Branch("genb1_phi", &genb1_phi, "genb1_phi/D");
        tree->Branch("genb1_energy", &genb1_energy, "genb1_energy/D");
        tree->Branch("genb1_flavor", &genb1_flavor, "genb1_flavor/D");
        tree->Branch("genb2_pt", &genb2_pt, "genb2_pt/D");
        tree->Branch("genb2_eta", &genb2_eta, "genb2_eta/D");
        tree->Branch("genb2_phi", &genb2_phi, "genb2_phi/D");
        tree->Branch("genb2_energy", &genb2_energy, "genb2_energy/D");
        tree->Branch("genb2_flavor", &genb2_flavor, "genb2_flavor/D");

        tree->Branch("genTauMisID", &genTauMisID, "genTauMisID/I");
        // gen Tau
        tree->Branch("gentau_dm",&gentau_dm,"gentau_dm/I");
        tree->Branch("gentau_pt",&gentau_pt,"gentau_pt/D");
        tree->Branch("gentau_eta",&gentau_eta,"gentau_eta/D");
        tree->Branch("gentau_phi",&gentau_phi,"gentau_phi/D");
        tree->Branch("gentau_energy",&gentau_energy,"gentau_energy/D");
        // PiChar
        tree->Branch("genPiChar_pt",&genPiChar_pt,"genPiChar_pt/D");
        tree->Branch("genPiChar_eta",&genPiChar_eta,"genPiChar_eta/D");
        tree->Branch("genPiChar_phi",&genPiChar_phi,"genPiChar_phi/D");
        tree->Branch("genPiChar_energy",&genPiChar_energy,"genPiChar_energy/D");
        // PiZero
        tree->Branch("genPi0_pt",&genPi0_pt,"genPi0_pt/D");
        tree->Branch("genPi0_eta",&genPi0_eta,"genPi0_eta/D");
        tree->Branch("genPi0_phi",&genPi0_phi,"genPi0_phi/D");
        tree->Branch("genPi0_energy",&genPi0_energy,"genPi0_energy/D");
        // Visible tau from generator
        tree->Branch("gentau_vis_pt",&gentau_vis_pt,"gentau_vis_pt/D");
        tree->Branch("gentau_vis_eta",&gentau_vis_eta,"gentau_vis_eta/D");
        tree->Branch("gentau_vis_phi",&gentau_vis_phi,"gentau_vis_phi/D");
        tree->Branch("gentau_vis_energy",&gentau_vis_energy,"gentau_vis_energy/D");
        // Neutrinos
        tree->Branch("nutau_pt",&nutau_pt,"nutau_pt/D");
        tree->Branch("nutau_eta",&nutau_eta,"nutau_eta/D");
        tree->Branch("nutau_phi",&nutau_phi,"nutau_phi/D");
        tree->Branch("nutau_energy",&nutau_energy,"nutau_energy/D");
        tree->Branch("nuW_pt",&nuW_pt,"nuW_pt/D");
        tree->Branch("nuW_eta",&nuW_eta,"nuW_eta/D");
        tree->Branch("nuW_phi",&nuW_phi,"nuW_phi/D");
        tree->Branch("nuW_energy",&nuW_energy,"nuW_energy/D");
        // Sum of all neutrinos in evenet
        tree->Branch("nNu",&nNu,"nNu/I");
        tree->Branch("SumNu_pt",&SumNu_pt,"SumNu_pt/D");
        tree->Branch("SumNu_eta",&SumNu_eta,"SumNu_eta/D");
        tree->Branch("SumNu_phi",&SumNu_phi,"SumNu_phi/D");
        tree->Branch("SumNu_energy",&SumNu_energy,"SumNu_energy/D");
    }
    // PuppiJets
    tree->Branch("nLooseBtagedPuppiJets", &nLooseBtagedPuppiJets, "nLooseBtagedPuppiJets/I");
    tree->Branch("PuppijetPtSum20PV", &PuppijetPtSum20PV, "PuppijetPtSum20PV/D");
    tree->Branch("PuppijetPtSum20", &PuppijetPtSum20, "PuppijetPtSum20/D");
    tree->Branch("nPuppiJets30", &nPuppiJets30, "nPuppiJets30/I");

    /*
    tree->Branch("PuppijetPtSum30", &PuppijetPtSum30, "PuppijetPtSum30/D");
    tree->Branch("PuppijetPtSum30PV", &PuppijetPtSum30PV, "PuppijetPtSum30PV/D");
    tree->Branch("nPuppiJets20", &nPuppiJets20, "nPuppiJets20/I");
    tree->Branch("nPuppiJets20PV", &nPuppiJets20PV, "nPuppiJets20PV/I");
    tree->Branch("nPuppiJets30PV", &nPuppiJets30PV, "nPuppiJets30PV/I");
    tree->Branch("nMediumBtagedPuppiJets", &nMediumBtagedPuppiJets, "nMediumBtagedPuppiJets/I");
    tree->Branch("nTightBtagedPuppiJets", &nTightBtagedPuppiJets, "nTightBtagedPuppiJets/I");
    tree->Branch("nLooseBtagedPuppiJetsPV", &nLooseBtagedPuppiJetsPV, "nLooseBtagedPuppiJetsPV/I");
    tree->Branch("nMediumBtagedPuppiJetsPV", &nMediumBtagedPuppiJetsPV, "nMediumBtagedPuppiJetsPV/I");
    tree->Branch("nTightBtagedPuppiJetsPV", &nTightBtagedPuppiJetsPV, "nTightBtagedPuppiJetsPV/I");
    */

    // Leading, subleading, bJet parameters
    tree->Branch("Jet1_pt", &Jet1_pt, "Jet1_pt/D");
    tree->Branch("Jet1_eta", &Jet1_eta, "Jet1_eta/D");
    tree->Branch("Jet1_phi", &Jet1_phi, "Jet1_phi/D");
    tree->Branch("Jet1_E", &Jet1_E, "Jet1_E/D");
    tree->Branch("Jet1_bprob", &Jet1_bprob, "Jet1_bprob/D");
    tree->Branch("Jet1_bbprob", &Jet1_bbprob, "Jet1_bbprob/D");
    tree->Branch("Jet1_hadronFlavour", &Jet1_hadronFlavour, "Jet1_hadronFlavour/I");
    //
    tree->Branch("Jet2_pt", &Jet2_pt, "Jet2_pt/D");
    tree->Branch("Jet2_eta", &Jet2_eta, "Jet2_eta/D");
    tree->Branch("Jet2_phi", &Jet2_phi, "Jet2_phi/D");
    tree->Branch("Jet2_E", &Jet2_E, "Jet2_E/D");
    tree->Branch("Jet2_bprob", &Jet2_bprob, "Jet2_bprob/D");
    tree->Branch("Jet2_bbprob", &Jet2_bbprob, "Jet2_bbprob/D");
    tree->Branch("Jet2_hadronFlavour", &Jet2_hadronFlavour, "Jet2_hadronFlavour/I");
    //
    tree->Branch("Jet3_pt", &Jet3_pt, "Jet3_pt/D");
    tree->Branch("Jet3_eta", &Jet3_eta, "Jet3_eta/D");
    tree->Branch("Jet3_phi", &Jet3_phi, "Jet3_phi/D");
    tree->Branch("Jet3_E", &Jet3_E, "Jet3_E/D");
    tree->Branch("Jet3_bprob", &Jet3_bprob, "Jet3_bprob/D");
    tree->Branch("Jet3_bbprob", &Jet3_bbprob, "Jet3_bbprob/D");
    tree->Branch("Jet3_hadronFlavour", &Jet3_hadronFlavour, "Jet3_hadronFlavour/I");
    //
    tree->Branch("Jet4_pt", &Jet4_pt, "Jet4_pt/D");
    tree->Branch("Jet4_eta", &Jet4_eta, "Jet4_eta/D");
    tree->Branch("Jet4_phi", &Jet4_phi, "Jet4_phi/D");
    tree->Branch("Jet4_E", &Jet4_E, "Jet4_E/D");
    tree->Branch("Jet4_bprob", &Jet4_bprob, "Jet4_bprob/D");
    tree->Branch("Jet4_bbprob", &Jet4_bbprob, "Jet4_bbprob/D");
    tree->Branch("Jet4_hadronFlavour", &Jet4_hadronFlavour, "Jet4_hadronFlavour/I");
    //
    tree->Branch("Jet5_pt", &Jet5_pt, "Jet5_pt/D");
    tree->Branch("Jet5_eta", &Jet5_eta, "Jet5_eta/D");
    tree->Branch("Jet5_phi", &Jet5_phi, "Jet5_phi/D");
    tree->Branch("Jet5_E", &Jet5_E, "Jet5_E/D");
    tree->Branch("Jet5_bprob", &Jet5_bprob, "Jet5_bprob/D");
    tree->Branch("Jet5_bbprob", &Jet5_bbprob, "Jet5_bbprob/D");
    tree->Branch("Jet5_hadronFlavour", &Jet5_hadronFlavour, "Jet5_hadronFlavour/I");
    //
    tree->Branch("Jet6_pt", &Jet6_pt, "Jet6_pt/D");
    tree->Branch("Jet6_eta", &Jet6_eta, "Jet6_eta/D");
    tree->Branch("Jet6_phi", &Jet6_phi, "Jet6_phi/D");
    tree->Branch("Jet6_E", &Jet6_E, "Jet6_E/D");
    tree->Branch("Jet6_bprob", &Jet6_bprob, "Jet6_bprob/D");
    tree->Branch("Jet6_bbprob", &Jet6_bbprob, "Jet6_bbprob/D");
    tree->Branch("Jet6_hadronFlavour", &Jet6_hadronFlavour, "Jet6_hadronFlavour/I");

    //Leptons
    tree->Branch("nLeptonCandidates", &nLeptonCandidates, "nLeptonCandidates/I");
    tree->Branch("lepton1_pt", &lepton1_pt, "lepton1_pt/D");
    tree->Branch("lepton1_E", &lepton1_E, "lepton1_E/D");
    tree->Branch("lepton1_eta", &lepton1_eta, "lepton1_eta/D");
    tree->Branch("lepton1_phi", &lepton1_phi, "lepton1_phi/D");
    tree->Branch("lepton1_dz", &lepton1_dz, "lepton1_dz/D");
    tree->Branch("lepton1_flavor", &lepton1_flavor, "lepton1_flavor/I");
    tree->Branch("lepton1_charge", &lepton1_charge, "lepton1_charge/I");
    tree->Branch("lepton1_trackIso", &lepton1_trackIso, "lepton1_trackIso/F");
    tree->Branch("lepton1_sumPuppiIso", &lepton1_sumPuppiIso, "lepton1_sumPuppiIso/F");
    tree->Branch("lepton1_sumPuppiNoLeptonIso", &lepton1_sumPuppiNoLeptonIso, "lepton1_sumPuppiNoLeptonIso/F");
    tree->Branch("lepton1Jet_pt", &lepton1Jet_pt, "lepton1Jet_pt/D");
    tree->Branch("lepton1Jet_eta", &lepton1Jet_eta, "lepton1Jet_eta/D");
    tree->Branch("lepton1Jet_phi", &lepton1Jet_phi, "lepton1Jet_phi/D");
    tree->Branch("lepton1Jet_E", &lepton1Jet_E, "lepton1Jet_E/D");
    tree->Branch("lepton1Jet_bprob", &lepton1Jet_bprob, "lepton1Jet_bprob/D");
    // electron
    //tree->Branch("lepton1_electron_caloIso", &lepton1_electron_caloIso, "lepton1_electron_caloIso/F");
    tree->Branch("lepton1_electron_cutBasedID_loose", &lepton1_electron_cutBasedID_loose, "lepton1_electron_cutBasedID_loose/I");
    tree->Branch("lepton1_electron_cutBasedID_medium", &lepton1_electron_cutBasedID_medium, "lepton1_electron_cutBasedID_medium/I");
    tree->Branch("lepton1_electron_cutBasedID_tight", &lepton1_electron_cutBasedID_tight, "lepton1_electron_cutBasedID_tight/I");
    tree->Branch("lepton1_electron_mvaNoIsoID_loose", &lepton1_electron_mvaNoIsoID_loose, "lepton1_electron_mvaNoIsoID_loose/I");
    tree->Branch("lepton1_electron_mvaNoIsoID_wp80", &lepton1_electron_mvaNoIsoID_wp80, "lepton1_electron_mvaNoIsoID_wp80/I");
    tree->Branch("lepton1_electron_mvaNoIsoID_wp90", &lepton1_electron_mvaNoIsoID_wp90, "electron_mvaNoIsoID_wp90/I");
    // muon
    //tree->Branch("lepton1_muon_caloIso", &lepton1_muon_caloIso, "lepton1_muon_caloIso/F");
    tree->Branch("lepton1_muon_CutBasedIdLoose", &lepton1_muon_CutBasedIdLoose, "lepton1_muon_CutBasedIdLoose/I");
    tree->Branch("lepton1_muon_CutBasedIdMedium", &lepton1_muon_CutBasedIdMedium, "lepton1_muon_CutBasedIdMedium/I");
    tree->Branch("lepton1_muon_CutBasedIdTight", &lepton1_muon_CutBasedIdTight, "lepton1_muon_CutBasedIdTight/I");
    tree->Branch("lepton1_muon_CutBasedIdGlobalHighPt", &lepton1_muon_CutBasedIdGlobalHighPt, "lepton1_muon_CutBasedIdGlobalHighPt/I");
    tree->Branch("lepton1_muon_CutBasedIdTrkHighPt", &lepton1_muon_CutBasedIdTrkHighPt, "lepton1_muon_CutBasedIdTrkHighPt/I");
    tree->Branch("lepton1_muon_PFIsoLoose", &lepton1_muon_PFIsoLoose, "lepton1_muon_PFIsoLoose/I");
    tree->Branch("lepton1_muon_PFIsoMedium", &lepton1_muon_PFIsoMedium, "lepton1_muon_PFIsoMedium/I");
    tree->Branch("lepton1_muon_PFIsoTight", &lepton1_muon_PFIsoTight, "lepton1_muon_PFIsoTight/I");
    tree->Branch("lepton1_muon_TkIsoLoose", &lepton1_muon_TkIsoLoose, "lepton1_muon_TkIsoLoose/I");
    tree->Branch("lepton1_muon_MvaLoose", &lepton1_muon_MvaLoose, "lepton1_muon_MvaLoose/I");
    tree->Branch("lepton1_muon_MvaMedium", &lepton1_muon_MvaMedium, "lepton1_muon_MvaMedium/I");
    tree->Branch("lepton1_muon_MvaTight", &lepton1_muon_MvaTight, "lepton1_muon_MvaTight/I");
    // Second tau lepton
    if (UseTau) {
        // Tau kinematic parameters
        tree->Branch("lepton1_tau_m",&lepton1_tau_m,"lepton1_tau_m/D");
        tree->Branch("lepton1_tau_dm",&lepton1_tau_dm,"lepton1_tau_dm/D");
        // Raw discriminators
        tree->Branch("lepton1_tau_absIso",&lepton1_tau_absIso,"lepton1_tau_absIso/D");
        tree->Branch("lepton1_tau_againstElectronRaw",&lepton1_tau_againstElectronRaw,"lepton1_tau_againstElectronRaw/D");
        tree->Branch("lepton1_tau_IsoMVArun2v1DBdR03oldDMwLTraw",&lepton1_tau_IsoMVArun2v1DBdR03oldDMwLTraw,"lepton1_tau_IsoMVArun2v1DBdR03oldDMwLTraw/D");
        /*
        tree->Branch("lepton1_tau_IsoMVArun2v1DBnewDMwLTraw",&lepton1_tau_IsoMVArun2v1DBnewDMwLTraw,"lepton1_tau_IsoMVArun2v1DBnewDMwLTraw/D");
        tree->Branch("lepton1_tau_IsoMVArun2v1DBoldDMwLTraw",&lepton1_tau_IsoMVArun2v1DBoldDMwLTraw,"lepton1_tau_IsoMVArun2v1DBoldDMwLTraw/D");
        tree->Branch("lepton1_tau_IsoMVArun2v1PWdR03oldDMwLTraw",&lepton1_tau_IsoMVArun2v1PWdR03oldDMwLTraw,"lepton1_tau_IsoMVArun2v1PWdR03oldDMwLTraw/D");
        tree->Branch("lepton1_tau_IsoMVArun2v1PWnewDMwLTraw",&lepton1_tau_IsoMVArun2v1PWnewDMwLTraw,"lepton1_tau_IsoMVArun2v1PWnewDMwLTraw/D");
        tree->Branch("lepton1_tau_IsoMVArun2v1PWoldDMwLTraw",&lepton1_tau_IsoMVArun2v1PWoldDMwLTraw,"lepton1_tau_IsoMVArun2v1PWoldDMwLTraw/D");
        */

        // Discriminators
        tree->Branch("lepton1_tau_looseCombinedIso",&lepton1_tau_looseCombinedIso,"lepton1_tau_looseCombinedIso/I");
        //tree->Branch("lepton1_tau_mediumCombinedIso",&lepton1_tau_mediumCombinedIso,"lepton1_tau_mediumCombinedIso/I");
        tree->Branch("lepton1_tau_tightCombinedIso",&lepton1_tau_tightCombinedIso,"lepton1_tau_tightCombinedIso/I");
        tree->Branch("lepton1_tau_looseMvaIso",&lepton1_tau_looseMvaIso,"lepton1_tau_looseMvaIso/I");
        tree->Branch("lepton1_tau_mediumMvaIso",&lepton1_tau_mediumMvaIso,"lepton1_tau_mediumMvaIso/I");
        tree->Branch("lepton1_tau_tightMvaIso",&lepton1_tau_tightMvaIso,"lepton1_tau_tightMvaIso/I");
        tree->Branch("lepton1_tau_VtightMvaIso",&lepton1_tau_VtightMvaIso,"lepton1_tau_VtightMvaIso/I");
        tree->Branch("lepton1_tau_VVtightMvaIso",&lepton1_tau_VVtightMvaIso,"lepton1_tau_VVtightMvaIso/I");
        tree->Branch("lepton1_tau_tightMuonRejection",&lepton1_tau_tightMuonRejection,"lepton1_tau_tightMuonRejection/I");
        tree->Branch("lepton1_tau_looseElectronRejection",&lepton1_tau_looseElectronRejection,"lepton1_tau_looseElectronRejection/I");
        tree->Branch("lepton1_tau_mediumElectronRejection",&lepton1_tau_mediumElectronRejection,"lepton1_tau_mediumElectronRejection/I");
        tree->Branch("lepton1_tau_tightElectronRejection",&lepton1_tau_tightElectronRejection,"lepton1_tau_tightElectronRejection/I");
        tree->Branch("lepton1_tau_VtightElectronRejection",&lepton1_tau_VtightElectronRejection,"lepton1_tau_VtightElectronRejection/I");

        // Deep 2017v2
        if (DeepTau) {
            tree->Branch("lepton1_tau_Deep2017v2ElectronRejection",&lepton1_tau_Deep2017v2ElectronRejection,"lepton1_tau_Deep2017v2ElectronRejection/D");
            tree->Branch("lepton1_tau_Deep2017v2MuonRejection",&lepton1_tau_Deep2017v2MuonRejection,"lepton1_tau_Deep2017v2MuonRejection/D");
            tree->Branch("lepton1_tau_Deep2017v2JetRejection",&lepton1_tau_Deep2017v2JetRejection,"lepton1_tau_Deep2017v2JetRejection/D");

            tree->Branch("lepton1_tau_VVLooseDeepTau2017v2VSjet",&lepton1_tau_VVLooseDeepTau2017v2VSjet,"lepton1_tau_VVLooseDeepTau2017v2VSjet/I");
            tree->Branch("lepton1_tau_VLooseDeepTau2017v2VSjet",&lepton1_tau_VLooseDeepTau2017v2VSjet,"lepton1_tau_VLooseDeepTau2017v2VSjet/I");
            tree->Branch("lepton1_tau_LooseDeepTau2017v2VSjet",&lepton1_tau_LooseDeepTau2017v2VSjet,"lepton1_tau_LooseDeepTau2017v2VSjet/I");
            tree->Branch("lepton1_tau_MediumDeepTau2017v2VSjet",&lepton1_tau_MediumDeepTau2017v2VSjet,"lepton1_tau_MediumDeepTau2017v2VSjet/I");
            tree->Branch("lepton1_tau_TightDeepTau2017v2VSjet",&lepton1_tau_TightDeepTau2017v2VSjet,"lepton1_tau_TightDeepTau2017v2VSjet/I");
            tree->Branch("lepton1_tau_VTightDeepTau2017v2VSjet",&lepton1_tau_VTightDeepTau2017v2VSjet,"lepton1_tau_VTightDeepTau2017v2VSjet/I");
            tree->Branch("lepton1_tau_VVTightDeepTau2017v2VSjet",&lepton1_tau_VVTightDeepTau2017v2VSjet,"lepton1_tau_VVTightDeepTau2017v2VSjet/I");
            tree->Branch("lepton1_tau_LooseDeepTau2017v2VSmu",&lepton1_tau_LooseDeepTau2017v2VSmu,"lepton1_tau_LooseDeepTau2017v2VSmu/I");
            tree->Branch("lepton1_tau_MediumDeepTau2017v2VSmu",&lepton1_tau_MediumDeepTau2017v2VSmu,"lepton1_tau_MediumDeepTau2017v2VSmu/I");
            tree->Branch("lepton1_tau_TightDeepTau2017v2VSmu",&lepton1_tau_TightDeepTau2017v2VSmu,"lepton1_tau_TightDeepTau2017v2VSmu/I");

            tree->Branch("lepton1_tau_VVLooseDeepTau2017v2VSe",&lepton1_tau_VVLooseDeepTau2017v2VSe,"lepton1_tau_VVLooseDeepTau2017v2VSe/I");
            tree->Branch("lepton1_tau_VLooseDeepTau2017v2VSe",&lepton1_tau_VLooseDeepTau2017v2VSe,"lepton1_tau_VLooseDeepTau2017v2VSe/I");
            tree->Branch("lepton1_tau_LooseDeepTau2017v2VSe",&lepton1_tau_LooseDeepTau2017v2VSe,"lepton1_tau_LooseDeepTau2017v2VSe/I");
            tree->Branch("lepton1_tau_MediumDeepTau2017v2VSe",&lepton1_tau_MediumDeepTau2017v2VSe,"lepton1_tau_MediumDeepTau2017v2VSe/I");
            tree->Branch("lepton1_tau_TightDeepTau2017v2VSe",&lepton1_tau_TightDeepTau2017v2VSe,"lepton1_tau_TightDeepTau2017v2VSe/I");
            tree->Branch("lepton1_tau_VTightDeepTau2017v2VSe",&lepton1_tau_VTightDeepTau2017v2VSe,"lepton1_tau_VTightDeepTau2017v2VSe/I");
            tree->Branch("lepton1_tau_VVTightDeepTau2017v2VSe",&lepton1_tau_VVTightDeepTau2017v2VSe,"lepton1_tau_VVTightDeepTau2017v2VSe/I");
        }
    }
    tree->Branch("lepton2_pt", &lepton2_pt, "lepton2_pt/D");
    tree->Branch("lepton2_E", &lepton2_E, "lepton2_E/D");
    tree->Branch("lepton2_eta", &lepton2_eta, "lepton2_eta/D");
    tree->Branch("lepton2_phi", &lepton2_phi, "lepton2_phi/D");
    tree->Branch("lepton2_dz", &lepton2_dz, "lepton2_dz/D");
    tree->Branch("lepton2_flavor", &lepton2_flavor, "lepton2_flavor/I");
    tree->Branch("lepton2_charge", &lepton2_charge, "lepton2_charge/I");
    tree->Branch("lepton2_trackIso", &lepton2_trackIso, "lepton2_trackIso/F");
    tree->Branch("lepton2_sumPuppiIso", &lepton2_sumPuppiIso, "lepton2_sumPuppiIso/F");
    tree->Branch("lepton2_sumPuppiNoLeptonIso", &lepton2_sumPuppiNoLeptonIso, "lepton2_sumPuppiNoLeptonIso/F");
    tree->Branch("lepton2Jet_pt", &lepton2Jet_pt, "lepton2Jet_pt/D");
    tree->Branch("lepton2Jet_eta", &lepton2Jet_eta, "lepton2Jet_eta/D");
    tree->Branch("lepton2Jet_phi", &lepton2Jet_phi, "lepton2Jet_phi/D");
    tree->Branch("lepton2Jet_E", &lepton2Jet_E, "lepton2Jet_E/D");
    tree->Branch("lepton2Jet_bprob", &lepton2Jet_bprob, "lepton2Jet_bprob/D");
    // electron
    //tree->Branch("lepton2_electron_caloIso", &lepton2_electron_caloIso, "lepton2_electron_caloIso/F");
    tree->Branch("lepton2_electron_cutBasedID_loose", &lepton2_electron_cutBasedID_loose, "lepton2_electron_cutBasedID_loose/I");
    tree->Branch("lepton2_electron_cutBasedID_medium", &lepton2_electron_cutBasedID_medium, "lepton2_electron_cutBasedID_medium/I");
    tree->Branch("lepton2_electron_cutBasedID_tight", &lepton2_electron_cutBasedID_tight, "lepton2_electron_cutBasedID_tight/I");
    tree->Branch("lepton2_electron_mvaNoIsoID_loose", &lepton2_electron_mvaNoIsoID_loose, "lepton2_electron_mvaNoIsoID_loose/I");
    tree->Branch("lepton2_electron_mvaNoIsoID_wp80", &lepton2_electron_mvaNoIsoID_wp80, "lepton2_electron_mvaNoIsoID_wp80/I");
    tree->Branch("lepton2_electron_mvaNoIsoID_wp90", &lepton2_electron_mvaNoIsoID_wp90, "electron_mvaNoIsoID_wp90/I");
    // muon
    //tree->Branch("lepton2_muon_caloIso", &lepton2_muon_caloIso, "lepton2_muon_caloIso/F");
    tree->Branch("lepton2_muon_CutBasedIdLoose", &lepton2_muon_CutBasedIdLoose, "lepton2_muon_CutBasedIdLoose/I");
    tree->Branch("lepton2_muon_CutBasedIdMedium", &lepton2_muon_CutBasedIdMedium, "lepton2_muon_CutBasedIdMedium/I");
    tree->Branch("lepton2_muon_CutBasedIdTight", &lepton2_muon_CutBasedIdTight, "lepton2_muon_CutBasedIdTight/I");
    tree->Branch("lepton2_muon_CutBasedIdGlobalHighPt", &lepton2_muon_CutBasedIdGlobalHighPt, "lepton2_muon_CutBasedIdGlobalHighPt/I");
    tree->Branch("lepton2_muon_CutBasedIdTrkHighPt", &lepton2_muon_CutBasedIdTrkHighPt, "lepton2_muon_CutBasedIdTrkHighPt/I");
    tree->Branch("lepton2_muon_PFIsoLoose", &lepton2_muon_PFIsoLoose, "lepton2_muon_PFIsoLoose/I");
    tree->Branch("lepton2_muon_PFIsoMedium", &lepton2_muon_PFIsoMedium, "lepton2_muon_PFIsoMedium/I");
    tree->Branch("lepton2_muon_PFIsoTight", &lepton2_muon_PFIsoTight, "lepton2_muon_PFIsoTight/I");
    tree->Branch("lepton2_muon_TkIsoLoose", &lepton2_muon_TkIsoLoose, "lepton2_muon_TkIsoLoose/I");
    tree->Branch("lepton2_muon_MvaLoose", &lepton2_muon_MvaLoose, "lepton2_muon_MvaLoose/I");
    tree->Branch("lepton2_muon_MvaMedium", &lepton2_muon_MvaMedium, "lepton2_muon_MvaMedium/I");
    tree->Branch("lepton2_muon_MvaTight", &lepton2_muon_MvaTight, "lepton2_muon_MvaTight/I");
    // Second tau lepton
    if (UseTau) {
        // Tau kinematic parameters
        tree->Branch("lepton2_tau_m",&lepton2_tau_m,"lepton2_tau_m/D");
        tree->Branch("lepton2_tau_dm",&lepton2_tau_dm,"lepton2_tau_dm/D");
        // Raw discriminators
        tree->Branch("lepton2_tau_absIso",&lepton2_tau_absIso,"lepton2_tau_absIso/D");
        tree->Branch("lepton2_tau_againstElectronRaw",&lepton2_tau_againstElectronRaw,"lepton2_tau_againstElectronRaw/D");
        tree->Branch("lepton2_tau_IsoMVArun2v1DBdR03oldDMwLTraw",&lepton2_tau_IsoMVArun2v1DBdR03oldDMwLTraw,"lepton2_tau_IsoMVArun2v1DBdR03oldDMwLTraw/D");
        /*
        tree->Branch("lepton2_tau_IsoMVArun2v1DBnewDMwLTraw",&lepton2_tau_IsoMVArun2v1DBnewDMwLTraw,"lepton2_tau_IsoMVArun2v1DBnewDMwLTraw/D");
        tree->Branch("lepton2_tau_IsoMVArun2v1DBoldDMwLTraw",&lepton2_tau_IsoMVArun2v1DBoldDMwLTraw,"lepton2_tau_IsoMVArun2v1DBoldDMwLTraw/D");
        tree->Branch("lepton2_tau_IsoMVArun2v1PWdR03oldDMwLTraw",&lepton2_tau_IsoMVArun2v1PWdR03oldDMwLTraw,"lepton2_tau_IsoMVArun2v1PWdR03oldDMwLTraw/D");
        tree->Branch("lepton2_tau_IsoMVArun2v1PWnewDMwLTraw",&lepton2_tau_IsoMVArun2v1PWnewDMwLTraw,"lepton2_tau_IsoMVArun2v1PWnewDMwLTraw/D");
        tree->Branch("lepton2_tau_IsoMVArun2v1PWoldDMwLTraw",&lepton2_tau_IsoMVArun2v1PWoldDMwLTraw,"lepton2_tau_IsoMVArun2v1PWoldDMwLTraw/D");
        */

        // Discriminators
        tree->Branch("lepton2_tau_looseCombinedIso",&lepton2_tau_looseCombinedIso,"lepton2_tau_looseCombinedIso/I");
        //tree->Branch("lepton2_tau_mediumCombinedIso",&lepton2_tau_mediumCombinedIso,"lepton2_tau_mediumCombinedIso/I");
        tree->Branch("lepton2_tau_tightCombinedIso",&lepton2_tau_tightCombinedIso,"lepton2_tau_tightCombinedIso/I");
        tree->Branch("lepton2_tau_looseMvaIso",&lepton2_tau_looseMvaIso,"lepton2_tau_looseMvaIso/I");
        tree->Branch("lepton2_tau_mediumMvaIso",&lepton2_tau_mediumMvaIso,"lepton2_tau_mediumMvaIso/I");
        tree->Branch("lepton2_tau_tightMvaIso",&lepton2_tau_tightMvaIso,"lepton2_tau_tightMvaIso/I");
        tree->Branch("lepton2_tau_VtightMvaIso",&lepton2_tau_VtightMvaIso,"lepton2_tau_VtightMvaIso/I");
        tree->Branch("lepton2_tau_VVtightMvaIso",&lepton2_tau_VVtightMvaIso,"lepton2_tau_VVtightMvaIso/I");
        tree->Branch("lepton2_tau_tightMuonRejection",&lepton2_tau_tightMuonRejection,"lepton2_tau_tightMuonRejection/I");
        tree->Branch("lepton2_tau_looseElectronRejection",&lepton2_tau_looseElectronRejection,"lepton2_tau_looseElectronRejection/I");
        tree->Branch("lepton2_tau_mediumElectronRejection",&lepton2_tau_mediumElectronRejection,"lepton2_tau_mediumElectronRejection/I");
        tree->Branch("lepton2_tau_tightElectronRejection",&lepton2_tau_tightElectronRejection,"lepton2_tau_tightElectronRejection/I");
        tree->Branch("lepton2_tau_VtightElectronRejection",&lepton2_tau_VtightElectronRejection,"lepton2_tau_VtightElectronRejection/I");

        // Deep 2017v2
        if (DeepTau) {
            tree->Branch("lepton2_tau_Deep2017v2ElectronRejection",&lepton2_tau_Deep2017v2ElectronRejection,"lepton2_tau_Deep2017v2ElectronRejection/D");
            tree->Branch("lepton2_tau_Deep2017v2MuonRejection",&lepton2_tau_Deep2017v2MuonRejection,"lepton2_tau_Deep2017v2MuonRejection/D");
            tree->Branch("lepton2_tau_Deep2017v2JetRejection",&lepton2_tau_Deep2017v2JetRejection,"lepton2_tau_Deep2017v2JetRejection/D");

            tree->Branch("lepton2_tau_VVLooseDeepTau2017v2VSjet",&lepton2_tau_VVLooseDeepTau2017v2VSjet,"lepton2_tau_VVLooseDeepTau2017v2VSjet/I");
            tree->Branch("lepton2_tau_VLooseDeepTau2017v2VSjet",&lepton2_tau_VLooseDeepTau2017v2VSjet,"lepton2_tau_VLooseDeepTau2017v2VSjet/I");
            tree->Branch("lepton2_tau_LooseDeepTau2017v2VSjet",&lepton2_tau_LooseDeepTau2017v2VSjet,"lepton2_tau_LooseDeepTau2017v2VSjet/I");
            tree->Branch("lepton2_tau_MediumDeepTau2017v2VSjet",&lepton2_tau_MediumDeepTau2017v2VSjet,"lepton2_tau_MediumDeepTau2017v2VSjet/I");
            tree->Branch("lepton2_tau_TightDeepTau2017v2VSjet",&lepton2_tau_TightDeepTau2017v2VSjet,"lepton2_tau_TightDeepTau2017v2VSjet/I");
            tree->Branch("lepton2_tau_VTightDeepTau2017v2VSjet",&lepton2_tau_VTightDeepTau2017v2VSjet,"lepton2_tau_VTightDeepTau2017v2VSjet/I");
            tree->Branch("lepton2_tau_VVTightDeepTau2017v2VSjet",&lepton2_tau_VVTightDeepTau2017v2VSjet,"lepton2_tau_VVTightDeepTau2017v2VSjet/I");
            tree->Branch("lepton2_tau_LooseDeepTau2017v2VSmu",&lepton2_tau_LooseDeepTau2017v2VSmu,"lepton2_tau_LooseDeepTau2017v2VSmu/I");
            tree->Branch("lepton2_tau_MediumDeepTau2017v2VSmu",&lepton2_tau_MediumDeepTau2017v2VSmu,"lepton2_tau_MediumDeepTau2017v2VSmu/I");
            tree->Branch("lepton2_tau_TightDeepTau2017v2VSmu",&lepton2_tau_TightDeepTau2017v2VSmu,"lepton2_tau_TightDeepTau2017v2VSmu/I");

            tree->Branch("lepton2_tau_VVLooseDeepTau2017v2VSe",&lepton2_tau_VVLooseDeepTau2017v2VSe,"lepton2_tau_VVLooseDeepTau2017v2VSe/I");
            tree->Branch("lepton2_tau_VLooseDeepTau2017v2VSe",&lepton2_tau_VLooseDeepTau2017v2VSe,"lepton2_tau_VLooseDeepTau2017v2VSe/I");
            tree->Branch("lepton2_tau_LooseDeepTau2017v2VSe",&lepton2_tau_LooseDeepTau2017v2VSe,"lepton2_tau_LooseDeepTau2017v2VSe/I");
            tree->Branch("lepton2_tau_MediumDeepTau2017v2VSe",&lepton2_tau_MediumDeepTau2017v2VSe,"lepton2_tau_MediumDeepTau2017v2VSe/I");
            tree->Branch("lepton2_tau_TightDeepTau2017v2VSe",&lepton2_tau_TightDeepTau2017v2VSe,"lepton2_tau_TightDeepTau2017v2VSe/I");
            tree->Branch("lepton2_tau_VTightDeepTau2017v2VSe",&lepton2_tau_VTightDeepTau2017v2VSe,"lepton2_tau_VTightDeepTau2017v2VSe/I");
            tree->Branch("lepton2_tau_VVTightDeepTau2017v2VSe",&lepton2_tau_VVTightDeepTau2017v2VSe,"lepton2_tau_VVTightDeepTau2017v2VSe/I");
        }
    }
    /*
    tree->Branch("VetoLeptons", &VetoLeptons, "VetoLeptons/I");
    tree->Branch("VetoTaus", &VetoTaus, "VetoTaus/I");
    */
    // Tau Jet
    tree->Branch("TauJet_pt", &TauJet_pt, "TauJet_pt/D");
    tree->Branch("TauJet_eta", &TauJet_eta, "TauJet_eta/D");
    tree->Branch("TauJet_phi", &TauJet_phi, "TauJet_phi/D");
    tree->Branch("TauJet_E", &TauJet_E, "TauJet_E/D");
    tree->Branch("TauJet_bprob", &TauJet_bprob, "TauJet_bprob/D");
    tree->Branch("TauJet_hadronFlavour", &TauJet_hadronFlavour, "TauJet_hadronFlavour/I");

    // MET
    tree->Branch("met", &met, "met/D");
    tree->Branch("met_phi", &met_phi, "met_phi/D");
    tree->Branch("met_eta", &met_eta, "met_eta/D");
    tree->Branch("met_significance", &met_significance, "met_significance/D");
    tree->Branch("met_mEtSig", &met_mEtSig, "met_mEtSig/D");
    tree->Branch("met_energy", &met_energy, "met_energy/D");
    /*
    tree->Branch("tauMET_mass", &tauMET_mass, "tauMET_mass/D");
    tree->Branch("m_t", &m_t, "m_t/D");
    tree->Branch("dPhi", &dPhi, "dPhi/D");
    */

    // Puppi MET
    tree->Branch("Puppimet", &Puppimet, "Puppimet/D");
    tree->Branch("Puppimet_phi", &Puppimet_phi, "Puppimet_phi/D");
    tree->Branch("Puppimet_eta", &Puppimet_eta, "Puppimet_eta/D");
    tree->Branch("Puppimet_significance", &Puppimet_significance, "Puppimet_significance/D");
    tree->Branch("Puppimet_metSig", &Puppimet_metSig, "Puppimet_metSig/D");
    tree->Branch("Puppimet_energy", &Puppimet_energy, "Puppimet_energy/D");
    /*
    tree->Branch("tauPuppimet_mass", &tauPuppimet_mass, "tauPuppimet_mass/D");
    tree->Branch("dPhiPuppimetTau", &dPhiPuppimetTau, "dPhiPuppimetTau/D");
    */
    //tree->Branch("m_t", &m_t, "m_t/D");
    
    // Vertices, tracks, tau candidates
    tree->Branch("nVtx",&nVtx,"nVtx/I");
    tree->Branch("nTrks",&nTrks,"nTrks/I");
    tree->Branch("nTau",&nTau,"nTau/I");
    tree->Branch("nTauC",&nTauC,"nTauC/I");

    tree->Branch("nGamma",&nGamma,"nGamma/I");
    tree->Branch("nPhotons",&nPhotons,"nPhotons/I");

    // TauSpinner
    /*
    tree->Branch("WT", &WT, "WT/D");
    tree->Branch("WTFlip", &WTFlip, "WTFlip/D");
    tree->Branch("WThminus", &WThminus, "WThminus/D");
    tree->Branch("WThplus", &WThplus, "WThplus/D");
    tree->Branch("TauSpinnerMother", &TauSpinnerMother, "TauSpinnerMother/I");
    */
    tree->Branch("nTargetTrigger1", &nTargetTrigger1, "nTargetTrigger1/I");
    tree->Branch("nTargetTrigger2", &nTargetTrigger2, "nTargetTrigger2/I");
    tree->Branch("nTargetTrigger3", &nTargetTrigger3, "nTargetTrigger3/I");
    tree->Branch("TriggerBit1", &TriggerBit1, "TriggerBit1/i");
    tree->Branch("TriggerBit2", &TriggerBit2, "TriggerBit2/i");
    tree->Branch("TriggerBit3", &TriggerBit3, "TriggerBit3/i");
    tree->Branch("TriggerBit4", &TriggerBit4, "TriggerBit4/i");
    tree->Branch("TriggerBit5", &TriggerBit5, "TriggerBit5/i");
    //
    tree->Branch("resultTriggerWeight",&resultTriggerWeight,"resultTriggerWeight/I");
    tree->Branch("triggerPrescaleHLT",&triggerPrescaleHLT,"triggerPrescaleHLT/I");
    tree->Branch("triggerPrescaleL1min",&triggerPrescaleL1min,"triggerPrescaleL1min/I");
    tree->Branch("triggerPrescaleL1max",&triggerPrescaleL1max,"triggerPrescaleL1max/I");
    //tree->Branch("WTisValid", &WTisValid, "WTisValid/D")
    // add more branches
    
    allTauPt = FS->make<TH1F>("allTauPt","allTauPt",300,0,300);
}

// ------------ method called once each job just after ending the event loop  ------------
void TTbarTauLepton::endJob() {
    if (monitoring){
        std::cout << "Events  = " << nEvent << std::endl;
        std::cout << "Trigger = " << nTrigger << std::endl;
        std::cout << "Tau     = " << nTau1 << std::endl;
        std::cout << "Lepton  = " << nLepton << std::endl;
        std::cout << "Jet     = " << nJet << std::endl;
        std::cout << "TauJet  = " << nTauJet << std::endl;
        std::cout << "Passed  = " << nPassed << std::endl;
    }
}

// ------------ method called when starting to processes a run  ------------

void TTbarTauLepton::beginRun(const edm::Run& Run, const edm::EventSetup& Setup) {
}

// ------------ method called when ending the processing of a run  ------------

void
TTbarTauLepton::endRun(edm::Run const& Run, edm::EventSetup const&)
{

}

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TTbarTauLepton::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TTbarTauLepton::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

//---------------------------------TRIGGER-------------------------------------------------
bool TTbarTauLepton::TriggerOK (const edm::Event& iEvent) {

    if (monitoringHLT) std::cout << std::endl << "!!! Triggers !!!" << std::endl;

    resultTriggerWeight = null;
    triggerPrescaleHLT = null;
    triggerPrescaleL1max = null;
    triggerPrescaleL1min = null;

    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(tok_triggerObjects, triggerObjects);
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    iEvent.getByToken(tok_triggerPrescales, triggerPrescales);
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1min;
    iEvent.getByToken(tok_triggerPrescalesL1min, triggerPrescalesL1min);
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesL1max;
    iEvent.getByToken(tok_triggerPrescalesL1max, triggerPrescalesL1max);
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByToken(tok_trigRes, triggerResults);

    std::vector<std::string> trigNameVec;
    std::vector<bool> trigPassVec;
    std::vector<int> trigPsVec;
    std::vector<int> trigL1minPsVec;
    std::vector<int> trigL1maxPsVec;
    std::vector<int> trigPrescaleVec;
    int finalPrescale = 0;

    bool triggerOK = false;
    nTargetTrigger1 = 0;
    nTargetTrigger2 = 0;
    nTargetTrigger3 = 0;
    TriggerBit1 = 0;
    TriggerBit2 = 0;
    TriggerBit3 = 0;
    TriggerBit4 = 0;
    TriggerBit5 = 0;
    /////////////////////////////TriggerResults////////////////////////////////////
    if (triggerResults.isValid()) {
        const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
        const std::vector<std::string> & triggerNames_ = triggerNames.triggerNames();
        for ( unsigned int iHLT=0; iHLT<triggerResults->size(); iHLT++ ) {
            int hlt    = triggerResults->accept(iHLT);
            //
            const std::string& trigName = triggerNames.triggerName(iHLT);
            int ps        = triggerPrescales->getPrescaleForIndex(iHLT);
            int psL1min   = triggerPrescalesL1min->getPrescaleForIndex(iHLT);
            int psL1max   = triggerPrescalesL1max->getPrescaleForIndex(iHLT);
            bool pass     = triggerResults->accept(iHLT);
            //std::cout << "# " << triggerNames_[iHLT] << std::endl;
            bool printed = false;
            if (monitoringHLT) {
                //std::cout << "Single Muon part 1" << std::endl;
                for ( unsigned int i=0; i<trigNames1.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames1[i].c_str())!= std::string::npos ) {
                        std::cout << iHLT << " (" << triggerResults->accept(iHLT) << ") " << trigNames1[i] <<
                        std::endl << "[" << i << "]" << " Single Muon part 1 " << std::endl;
                        printed = true;
                    }
                }
                //std::cout << "Single Muon part 2" << std::endl;
                for ( unsigned int i=0; i<trigNames2.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames2[i].c_str())!= std::string::npos ) {
                        std::cout << iHLT << " (" << triggerResults->accept(iHLT) << ") " << trigNames2[i] <<
                        std::endl << "[" << i << "]" << " Single Muon part 2 " << std::endl;
                        printed = true;
                    }
                }
                //std::cout << "Single Electron part 1" << std::endl;
                for ( unsigned int i=0; i<trigNames3.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames3[i].c_str())!= std::string::npos ) {
                        std::cout << iHLT << " (" << triggerResults->accept(iHLT) << ") " << trigNames3[i] <<
                        std::endl << "[" << i << "]" << " Single Electron part 1 " << std::endl;
                        printed = true;
                    }
                }
                //std::cout << "Single Electron part 2" << std::endl;
                for ( unsigned int i=0; i<trigNames4.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames4[i].c_str())!= std::string::npos ) {
                        std::cout << iHLT << " (" << triggerResults->accept(iHLT) << ") " << trigNames4[i] <<
                        std::endl << "[" << i << "]" << " Single Electron part 2 " << std::endl;
                        printed = true;
                    }
                }
                //std::cout << "Tau" << std::endl;
                for ( unsigned int i=0; i<trigNames5.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames5[i].c_str())!= std::string::npos ) {
                        std::cout << iHLT << " (" << triggerResults->accept(iHLT) << ") " << trigNames5[i] <<
                        std::endl << "[" << i << "]" << " Tau " << std::endl;
                        printed = true;
                    }
                }
                if (!printed) {
                    std::cout << iHLT << " (" << triggerResults->accept(iHLT) << ") " << triggerNames_[iHLT] << std::endl;
                }
            }
            if (hlt > 0 && monitoringHLT) {
                /*
                std::cout << "Name    = " << trigName << std::endl;
                std::cout << "ps      = " << ps << std::endl;
                std::cout << "psL1min = " << psL1min << std::endl;
                std::cout << "psL1max = " << psL1max << std::endl;
                */
            }
            if ( hlt > 0 ) {
                // Write prescale weights to corresponding vectors
                trigPsVec.push_back(ps);
                trigL1minPsVec.push_back(psL1min);
                trigL1maxPsVec.push_back(psL1max);
                trigNameVec.push_back(trigName);
                trigPassVec.push_back(pass);
                trigPrescaleVec.push_back(ps * psL1min);
                if (monitoringHLT) std::cout << triggerNames_[iHLT] << std::endl;
                for ( unsigned int i=0; i<trigNamesTarget1.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNamesTarget1[i].c_str())!= std::string::npos ) {
                        nTargetTrigger1++;
                        if (monitoringHLT) std::cout << "Target Trigger" << std::endl;
                    }
                }
                for ( unsigned int i=0; i<trigNamesTarget2.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNamesTarget2[i].c_str())!= std::string::npos ) {
                        nTargetTrigger2++;
                        if (monitoringHLT) std::cout << "Target Trigger 2" << std::endl;
                    }
                }
                for ( unsigned int i=0; i<trigNamesTarget3.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNamesTarget3[i].c_str())!= std::string::npos ) {
                        nTargetTrigger3++;
                        if (monitoringHLT) std::cout << "Target Trigger 3" << triggerNames_[iHLT] << std::endl;
                    }
                }
                for ( unsigned int i=0; i<trigNames1.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames1[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit1 = TriggerBit1 + TMath::Power(2, i_signed);
                        if (monitoringHLT) std::cout << "Trigger" << triggerNames_[iHLT] << std::endl;
                    }
                }
                for ( unsigned int i=0; i<trigNames2.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames2[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit2 = TriggerBit2 + TMath::Power(2, i_signed);
                        if (monitoringHLT) std::cout << "Trigger " << triggerNames_[iHLT] << std::endl;
                    }
                }
                for ( unsigned int i=0; i<trigNames3.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames3[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit3 = TriggerBit3 + TMath::Power(2, i_signed);
                        if (monitoringHLT) std::cout << "Trigger " << triggerNames_[iHLT] << std::endl;
                    }
                }
                for ( unsigned int i=0; i<trigNames4.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames4[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit4 = TriggerBit4 + TMath::Power(2, i_signed);
                        if (monitoringHLT) std::cout << "Trigger " << triggerNames_[iHLT] << std::endl;
                    }
                }
                for ( unsigned int i=0; i<trigNames5.size(); ++i ) {
                    if ( triggerNames_[iHLT].find(trigNames5[i].c_str())!= std::string::npos ) {
                        int i_signed = i;
                        TriggerBit5 = TriggerBit5 + TMath::Power(2, i_signed);
                        if (monitoringHLT) std::cout << "Trigger " << triggerNames_[iHLT] << std::endl;
                    }
                }
            }
        }
        /*
        for ( unsigned int iHLT=0; iHLT<triggerResults->size(); iHLT++ ) {
            int hlt    = triggerResults->accept(iHLT);
            if (nJetHTTriggers < 1) {
                if (hlt > 0) std::cout << triggerNames_[iHLT] << std::endl;
            }
        }
        */
    }

    int requiredTrigPrescale = 9999999;
    int requiredTrigIndex = -1;
    std::string requiredTrigName;
    // Loop over prescale vectors
    for(unsigned int l = 0; l < trigPrescaleVec.size(); l++) {
        if (!trigPassVec[l]) continue;
        if (trigPrescaleVec[l] < requiredTrigPrescale) {
            requiredTrigPrescale = trigPrescaleVec[l];
            requiredTrigIndex = l;
            requiredTrigName = trigNameVec[l];
        }
        /*
        if (monitoringHLT) {
            std::cout << "###################" << std::endl;
            std::cout << "current Index    = " << l << std::endl;
            std::cout << "current Name     = " << trigNameVec[l] << std::endl;
            std::cout << "current Prescale = " << trigPrescaleVec[l] << std::endl;
        }
        */
    }
    if (monitoringHLT) {
        std::cout << "Single Muon 1 (" << TriggerBit1 << ")" << std::endl <<
        decimal_to_binary_string(TriggerBit1) << std::endl;
        std::cout << "Single Muon 2 (" << TriggerBit2 << ")" << std::endl <<
        decimal_to_binary_string(TriggerBit2) << std::endl;
        std::cout << "Single Electron 1 (" << TriggerBit3 << ")" << std::endl <<
        decimal_to_binary_string(TriggerBit3) << std::endl;
        std::cout << "Single Electron 2 (" << TriggerBit4 << ")" << std::endl <<
        decimal_to_binary_string(TriggerBit4) << std::endl;
        std::cout << "Tau (" << TriggerBit5 << ")" << std::endl <<
        decimal_to_binary_string(TriggerBit5) << std::endl;
    }

    if (requiredTrigIndex > -1) {
        resultTriggerWeight = requiredTrigPrescale;
        triggerPrescaleHLT = trigPsVec[requiredTrigIndex];
        triggerPrescaleL1max = trigL1maxPsVec[requiredTrigIndex];
        triggerPrescaleL1min = trigL1minPsVec[requiredTrigIndex];
    }
    if (monitoringHLT) {
        std::cout << "resultTriggerWeight = " << resultTriggerWeight << std::endl
        << "triggerPrescaleHLT = " << triggerPrescaleHLT << std::endl
        << "triggerPrescaleL1max = " << triggerPrescaleL1max << std::endl
        << "triggerPrescaleL1min = " << triggerPrescaleL1min << std::endl;
    }
    return true;
};

//-----------------------------------------------------------------------------------------

bool TTbarTauLepton::AddTau(const edm::Event& event) {

    if (monitoringTau) std::cout << std::endl << "!!! Tau !!!" << std::endl;

    tau_found = 0;
    nTauC     = 0;
    math::XYZTLorentzVector DiPhoton_p4;
    int nPhotons_temp = 0;

    edm::Handle<pat::TauCollection> taus;
    event.getByToken(TauCollectionToken_, taus);
    if (monitoringTau) {
        std::cout << "Tau collection is valid = " << taus.isValid() << std::endl;
        std::cout << "Tau collection size     = " << taus->size() << std::endl;
    }
    if (!taus.isValid() || taus->size() == 0) return false;

    edm::Handle<reco::VertexCompositePtrCandidateCollection> SecondaryVertices;
    event.getByToken(SVToken_, SecondaryVertices);

    // search for the tau candidate with the minimum isolation and the
    // maximum transverse momentum
    size_t index  = taus->size();;
    //int index_2 = null;
    //std::map <int, std::pair <double, double>> TauIdMap;

    ///////////////////////////////////////////////////////////////////////////////////
    std::map <int, bool> TauIndexMap;
    for (size_t i = 0; i < taus->size(); ++i) {
#define cut(condition) if (!(condition)) continue;
        TauIndexMap.insert(std::pair<int, bool>(i, false));
        auto& tau = (*taus)[i];
        if (monitoringTau) {
            std::cout << "######### Tau candidate number " << i << " ################" << std::endl;
            std::cout << "pt = " << tau.pt() << ", eta = " << tau.eta() << ", "
            << "IsoIdMVA = " << tau.tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw") << ", "
            << "AbsIso   = " << tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") << ", "
            << "Jet = " << tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT") << ", "
            << "Ele = " << tau.tauID("againstElectronVLooseMVA6") << ", "
            << "Mu = " << tau.tauID("againstMuonLoose3") << ", "
            << "dZ = " << (pv_position - tau.vertex()).R() << std::endl;
        }
        cut(tau.pt() > tauPtMin); // tau transverse momentum
        if (monitoringTau) std::cout << "Pt thr passed" << std::endl;
        cut(TMath::Abs(tau.eta()) < tauEtaMax); // tau pseudorapidity
        if (monitoringTau) std::cout << "|Eta| thr passed" << std::endl;
        if (looseTauID && !DeepTau) {
            cut(tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT")); // at least VVloose Iso MVA
            if (monitoringTau) std::cout << "VLoose jet thr passed" << std::endl;
            cut(tau.tauID("againstElectronVLooseMVA6")); // at least loose Iso MVA
            if (monitoringTau) std::cout << "VLoose ele thr passed" << std::endl;
            cut(tau.tauID("againstMuonLoose3")); // at least loose Iso MVA
            if (monitoringTau) std::cout << "Loose mu thr passed" << std::endl;
            //cut(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")); // at least loose Iso MVA
        } else if (DeepTau) {
            cut(tau.tauID("byVLooseDeepTau2017v2VSmu"));
            cut(tau.tauID("byVVVLooseDeepTau2017v2VSe"));
            cut(tau.tauID("byVVVLooseDeepTau2017v2VSjet"));
        }

        ++nTauC;
        allTauPt->Fill(tau.pt());
        //std::pair <int, double> PairAbsIsoPt (tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"), tau.pt());
        //TauIdMap.insert(std::pair<int, std::pair<double, double>> (i, PairAbsIsoPt));

        cut((pv_position - tau.vertex()).R() < tauDzMax); // tau vertex displacement
        if (monitoringTau) std::cout << "Vertex displacement thr passed" << std::endl;
        if (monitoringTau) std::cout << "PV passed" << std::endl;
        if (monitoringTau) std::cout << "candidate passed cuts" << std::endl;
        TauIndexMap.erase(i);
        TauIndexMap.insert(std::pair<int, bool>(i, true));
        if (monitoringTau) std::cout << "Index map pair inserted" << std::endl;

        if (index < taus->size()) {
            /*
            double discriminator = (*taus)[i].tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw");
            double maxDiscriminator = (*taus)[index].tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw");
            cut(
                  discriminator > maxDiscriminator
                  || discriminator == maxDiscriminator && tau.pt() > (*taus)[index].pt()
            )
            */
            double iso = (*taus)[i].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
            double minIso = (*taus)[index].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
            cut(
                  iso < minIso
                  || iso == minIso && tau.pt() > (*taus)[index].pt()
            )
        }

        index = i;
#undef cut
    };
    ///////////////////////////////////////////////////////////////////////////////////


    if (monitoringTau) {
        std::cout << "Tau index (size) = " << index << "(" << taus->size() << ")" << std::endl;
        //for (auto iter = TauIndexMap.begin(); iter != TauIndexMap.end(); ++iter) {
        //  std::cout << "tau index = " << iter->first << ", value = " << iter->second << std::endl;
        //}
        for (unsigned l = 0; l < TauIndexMap.size(); l++) {
            std::cout << "tau index = " << l << ", value = " << TauIndexMap[l] << std::endl;
        }
    }

    // I think it's not nesessary because there may be only one tau in collection
    if (index == taus->size()) return false;

    auto& tau = (*taus)[index];
    if (monitoringTau) TauMonitor (tau, DeepTau, pv_position);

    const reco::CandidatePtr leadPiCh = tau.leadChargedHadrCand();
    reco::CandidatePtrVector VectorSignalCands = tau.signalCands();
    reco::CandidatePtrVector VectorPiCh = tau.signalChargedHadrCands();

    reco::CandidatePtrVector VectorSignalCands_Photons;
    VectorSignalCands_Photons.clear();
    //std::vector<reco::PFCandidate> newPhotons;
    //std::vector<reco::PFCandidatePtr> newPhotonsPtrs;
    //newPhotons.clear();
    // Fill vector of Photons for chosen tau lepton
    // Also calcolate dR^2 variable
    double tau_dR2_numerator = 0;
    double tau_dR2_denominator = 0;
    if (VectorSignalCands.size() > 0) {
        for (unsigned l = 0; l < VectorSignalCands.size(); l++) {
            tau_dR2_numerator += sqr(VectorSignalCands[l]->pt()) * (sqr(dphi(tau.phi(), VectorSignalCands[l]->phi())) + sqr(tau.eta() - VectorSignalCands[l]->eta()));
            tau_dR2_denominator += sqr(VectorSignalCands[l]->pt());
            if (VectorSignalCands[l]->pdgId() == 22) {
                VectorSignalCands_Photons.push_back(VectorSignalCands[l]);
            }
        }
    }
    if (tau_dR2_denominator > 0) tau_dR2 = tau_dR2_numerator/tau_dR2_denominator;
    else tau_dR2 = null;

    // Calculate DiPhoton kinematic
    math::XYZTLorentzVector Photons_p4;
    double InvMass_temp = null;

    using namespace reco;
    using namespace tau;

    std::vector<reco::RecoTauPiZero> CombPi0Vector;
    CombPi0Vector.clear();
    AddFourMomenta p4Builder;
    unsigned choose = 2;
    unsigned maxInputGammas = 10;
    double maxMass = -1;
    double minMass = 0;

    if (monitoringTau) std::cout << "--- Photons parameters:" << std::endl;
    if (VectorSignalCands_Photons.size() > 0) {
        for (unsigned l = 0; l < VectorSignalCands_Photons.size(); l++) {
            Photons_p4 += VectorSignalCands_Photons[l]->p4();
            for (unsigned n = l; n < VectorSignalCands_Photons.size(); n++) {
                if (n == l) continue;
                // Pi zero candidate
                const Candidate::LorentzVector totalP4;
                reco::RecoTauPiZero piZero(0, totalP4, Candidate::Point(0, 0, 0), 111, 10001, true, RecoTauPiZero::kCombinatoric);
                piZero.addDaughter(VectorSignalCands_Photons[l]);
                piZero.addDaughter(VectorSignalCands_Photons[n]);
                p4Builder.set(piZero);
                if (piZero.daughterPtr(0).isNonnull()) piZero.setVertex(piZero.daughterPtr(0)->vertex());
                if ((maxMass < 0 || piZero.mass() < maxMass) && piZero.mass() > minMass) CombPi0Vector.push_back(piZero);
                //
                double InvMass_pi0 = (VectorSignalCands_Photons[l]->p4() + VectorSignalCands_Photons[n]->p4()).M();
                //if (monitoringTau) std::cout << "mass Photon [" << l <<"][" << n << "] = " << InvMass_pi0 << std::endl;
                if (l == 0 && n == 1) {
                    InvMass_temp = InvMass_pi0;
                    DiPhoton_p4 = VectorSignalCands_Photons[l]->p4() + VectorSignalCands_Photons[n]->p4();
                } else if (abs(InvMass_pi0 - 0.135) < abs(InvMass_temp - 0.135)) {
                    InvMass_temp = InvMass_pi0;
                    DiPhoton_p4 = VectorSignalCands_Photons[l]->p4() + VectorSignalCands_Photons[n]->p4();
                }
            }
        }
    }

    nPhotons_temp = VectorSignalCands_Photons.size();

    if (monitoringTau) {
        std::cout << "Vector of Combined Pi zeros:" << std::endl;
        for (unsigned n = 0; n < CombPi0Vector.size(); n++) {
            std::cout << "Pi0 [" << n << "]: Pt = " << CombPi0Vector[n].pt() << ", mass = " << CombPi0Vector[n].mass() << std::endl;
        }
    }

    // ----------------------------------------------------------------------------------------------

    // Add second tau candidate if there is
    /*
    if (nTauC > 1) {
        for (auto MapId_iter = TauIdMap.begin(); MapId_iter !=  TauIdMap.end(); ++MapId_iter) {
            if ((*MapId_iter).second.first < tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") ||
                (*MapId_iter).second.first == tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") && (*MapId_iter).second.second < tau.pt()) {
                index_2 = (*MapId_iter).first;
            }
        }
    } else {
        index_2 = null;
    }
    */

    //reco::CandidateCollection Tau_daughters = tau.daughters;


    tau_pt                     = tau.pt();
    tau_eta                    = tau.eta();
    tau_phi                    = tau.phi();
    tau_dm                     = tau.decayMode();
    tau_dz                     = (pv_position - tau.vertex()).R();
    tau_q                      = tau.charge();
    tau_m                      = tau.mass();
    
    //if (monitoring) std::cout << "Selected tau discriminator (" << index << ") = " << tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") << std::endl;
    tau_absIso                  = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    tau_looseCombinedIso        = tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
    tau_mediumCombinedIso       = tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
    tau_tightCombinedIso        = tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
    tau_VlooseMvaIso            = tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
    tau_looseMvaIso             = tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT"); // there are more possible discriminators 
    tau_mediumMvaIso            = tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
    tau_tightMvaIso             = tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
    tau_VtightMvaIso            = tau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
    tau_VVtightMvaIso           = tau.tauID("byVVTightIsolationMVArun2v1DBoldDMwLT");
    // New Raw discriminators
    tau_againstElectronRaw            = tau.tauID("againstElectronMVA6Raw");
    tau_IsoMVArun2v1DBdR03oldDMwLTraw = tau.tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw");
    tau_IsoMVArun2v1DBnewDMwLTraw     = tau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
    tau_IsoMVArun2v1DBoldDMwLTraw     = tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
    tau_IsoMVArun2v1PWdR03oldDMwLTraw = tau.tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw");
    tau_IsoMVArun2v1PWnewDMwLTraw     = tau.tauID("byIsolationMVArun2v1PWnewDMwLTraw");
    tau_IsoMVArun2v1PWoldDMwLTraw     = tau.tauID("byIsolationMVArun2v1PWoldDMwLTraw");
    //
    tau_tightMuonRejection      = tau.tauID("againstMuonTight3");
    tau_looseElectronRejection  = tau.tauID("againstElectronLooseMVA6");
    tau_mediumElectronRejection = tau.tauID("againstElectronMediumMVA6");
    tau_tightElectronRejection  = tau.tauID("againstElectronTightMVA6");
    tau_VtightElectronRejection = tau.tauID("againstElectronVTightMVA6");
    decayModeFindingNewDMs      = tau.tauID("decayModeFindingNewDMs");
    decayModeFinding            = tau.tauID("decayModeFinding");
    if (DeepTau) {
        tau_Deep2017v2ElectronRejection = tau.tauID("byDeepTau2017v2VSeraw");
        tau_Deep2017v2MuonRejection     = tau.tauID("byDeepTau2017v2VSmuraw");
        tau_Deep2017v2JetRejection      = tau.tauID("byDeepTau2017v2VSjetraw");
        
        tau_VVLooseDeepTau2017v2VSjet   = tau.tauID("byVVLooseDeepTau2017v2VSjet");
        tau_VLooseDeepTau2017v2VSjet    = tau.tauID("byVLooseDeepTau2017v2VSjet");
        tau_LooseDeepTau2017v2VSjet     = tau.tauID("byLooseDeepTau2017v2VSjet");
        tau_MediumDeepTau2017v2VSjet    = tau.tauID("byMediumDeepTau2017v2VSjet");
        tau_TightDeepTau2017v2VSjet     = tau.tauID("byTightDeepTau2017v2VSjet");
        tau_VTightDeepTau2017v2VSjet    = tau.tauID("byVTightDeepTau2017v2VSjet");
        tau_VVTightDeepTau2017v2VSjet   = tau.tauID("byVVTightDeepTau2017v2VSjet");

        tau_LooseDeepTau2017v2VSmu      = tau.tauID("byLooseDeepTau2017v2VSmu");
        tau_MediumDeepTau2017v2VSmu     = tau.tauID("byMediumDeepTau2017v2VSmu");
        tau_TightDeepTau2017v2VSmu      = tau.tauID("byTightDeepTau2017v2VSmu");

        tau_VVLooseDeepTau2017v2VSe     = tau.tauID("byVVLooseDeepTau2017v2VSe");
        tau_VLooseDeepTau2017v2VSe      = tau.tauID("byVLooseDeepTau2017v2VSe");
        tau_LooseDeepTau2017v2VSe       = tau.tauID("byLooseDeepTau2017v2VSe");
        tau_MediumDeepTau2017v2VSe      = tau.tauID("byMediumDeepTau2017v2VSe");
        tau_TightDeepTau2017v2VSe       = tau.tauID("byTightDeepTau2017v2VSe");
        tau_VTightDeepTau2017v2VSe      = tau.tauID("byVTightDeepTau2017v2VSe");
        tau_VVTightDeepTau2017v2VSe     = tau.tauID("byVVTightDeepTau2017v2VSe");
    } else {
        tau_Deep2017v2ElectronRejection = null;
        tau_Deep2017v2MuonRejection     = null;
        tau_Deep2017v2JetRejection      = null;
    }

    piChar_pt = leadPiCh->pt();
    piChar_eta = leadPiCh->eta();
    piChar_phi = leadPiCh->phi();
    piChar_m = leadPiCh->mass();
    piChar_q = leadPiCh->charge();
    
    math::XYZTLorentzVector tau_p4 = tau.p4();
    math::XYZTLorentzVector piChar_p4 = leadPiCh->p4();

    piZero_pt  = (tau_p4 - piChar_p4).pt();
    piZero_eta = (tau_p4 - piChar_p4).eta();
    piZero_phi = (tau_p4 - piChar_p4).phi();
    piZero_m   = (tau_p4 - piChar_p4).M();

    if (tau.signalGammaCands().size() > 1) {
        DiPhoton_pt  = DiPhoton_p4.pt();
        DiPhoton_eta = DiPhoton_p4.eta();
        DiPhoton_phi = DiPhoton_p4.phi();
        DiPhoton_m   = DiPhoton_p4.M();
        TauPhoton1_pt      = tau.signalGammaCands()[0]->pt();
        TauPhoton1_eta     = tau.signalGammaCands()[0]->eta();
        TauPhoton1_phi     = tau.signalGammaCands()[0]->phi();
        TauPhoton1_energy  = tau.signalGammaCands()[0]->energy();
        TauPhoton2_pt      = tau.signalGammaCands()[1]->pt();
        TauPhoton2_eta     = tau.signalGammaCands()[1]->eta();
        TauPhoton2_phi     = tau.signalGammaCands()[1]->phi();
        TauPhoton2_energy  = tau.signalGammaCands()[1]->energy();
    } else {
        DiPhoton_pt  = null;
        DiPhoton_eta = null;
        DiPhoton_phi = null;
        DiPhoton_m   = null;
        if (tau.signalGammaCands().size() == 1) {
            TauPhoton1_pt      = tau.signalGammaCands()[0]->pt();
            TauPhoton1_eta     = tau.signalGammaCands()[0]->eta();
            TauPhoton1_phi     = tau.signalGammaCands()[0]->phi();
            TauPhoton1_energy  = tau.signalGammaCands()[0]->energy();
            TauPhoton2_pt      = null;
            TauPhoton2_eta     = null;
            TauPhoton2_phi     = null;
            TauPhoton2_energy  = null;
        } else {
            TauPhoton1_pt      = null;
            TauPhoton1_eta     = null;
            TauPhoton1_phi     = null;
            TauPhoton1_energy  = null;
            TauPhoton2_pt      = null;
            TauPhoton2_eta     = null;
            TauPhoton2_phi     = null;
            TauPhoton2_energy  = null;
        }
    }


    // SV analysis
    ///const reco::VertexCompositePtrCandidate *TauSVCandidate = nullptr;
    tau_hasSV     = null;
    tau_SVdR      = null;
    SV_Chi2       = null;
    SV_Chi2NDF    = null;
    tauPVtoSVdPhi = null;

    if (VectorPiCh.size() > 0) {
        for (unsigned n = 0; n < VectorPiCh.size(); n++) {
           if (monitoringTau) std::cout << "Charged Candidate " << n << ": pdgID = " << VectorPiCh[n]->pdgId() << ", Pt = " << VectorPiCh[n]->pt() << std::endl;
            //std::cout << "Vertex = (" << VectorPiCh[n]->vertex().x() << ", " << VectorPiCh[n]->vertex().y() << ", " << VectorPiCh[n]->vertex().z() << ")" << std::endl;
            for (unsigned nSV = 0; nSV < SecondaryVertices->size(); nSV++) {
                double tau_SVdRtemp = ((*SecondaryVertices)[nSV].position() - VectorPiCh[n]->vertex()).R();
                if (tau_SVdRtemp < abs(tau_SVdR)) {
                    ///TauSVCandidate = &(*SecondaryVertices)[nSV];
                    tau_SVdR    = tau_SVdRtemp;
                    SV_position = (*SecondaryVertices)[nSV].position();
                    SV_Chi2     = (*SecondaryVertices)[nSV].vertexChi2();
                    SV_Chi2NDF  = (*SecondaryVertices)[nSV].vertexNormalizedChi2();
                    tau_hasSV   = 1; 
                    if (monitoringTau) {
                        std::cout << "deltaR = " << ((*SecondaryVertices)[nSV].position() - VectorPiCh[n]->vertex()).R() << std::endl; 
                        std::cout << "vertex Chi2 = " << (*SecondaryVertices)[nSV].vertexChi2() << std::endl;
                        std::cout << "vertex Chi2/NDF = " << (*SecondaryVertices)[nSV].vertexNormalizedChi2() << std::endl;
                    }
                }
                //std::cout << "SV [" << nSV << "] position = " << "(" << (*SecondaryVertices)[nSV].position().x() << ", " << (*SecondaryVertices)[nSV].position().y() << ", " << (*SecondaryVertices)[nSV].position().z() << ")" << std::endl;
                //std::cout << "deltaR = " << ((*SecondaryVertices)[nSV].position() - VectorPiCh[n]->vertex()).R() << std::endl; 
            }
        }
    }

    if (tau.hasSecondaryVertex() && tau_hasSV > 0) tau_hasSV = 2;
    if (tau.hasSecondaryVertex() && tau_hasSV < 0) tau_hasSV = -1;
    if (tau_hasSV > 0) pv_SVdR = (pv_position - SV_position).R();

    ROOT::Math::Cartesian3D <double> Tau_vector;
    Tau_vector.SetXYZ(tau.px(), tau.py(), tau.pz());
    ROOT::Math::Cartesian3D <double> PVtoSV_vector;
    PVtoSV_vector.SetXYZ((pv_position - SV_position).x(), (pv_position - SV_position).y(), (pv_position - SV_position).z());
    if (tau_hasSV > 0) tauPVtoSVdPhi = dphi(Tau_vector.Phi(), PVtoSV_vector.Phi());

    if (monitoringTau) {
        std::cout << "Final SV:" << std::endl;
        std::cout << "deltaR          = " << pv_SVdR << std::endl; 
        std::cout << "vertex Chi2     = " << SV_Chi2 << std::endl;
        std::cout << "vertex Chi2/NDF = " << SV_Chi2NDF << std::endl;
        std::cout << "dPhi            = " << tauPVtoSVdPhi << std::endl;
    }

    nPhotons  = nPhotons_temp;
    nGamma    = tau.signalGammaCands().size();
    nTau      = taus->size();
    tau_found = 1;

    VectorSignalCands.clear();
    VectorPiCh.clear();
    VectorSignalCands_Photons.clear();

    return true;
};

void TTbarTauLepton::FindGenTau(const edm::Event& event) {

    if (monitoringGen) std::cout << std::endl << "!!! Generated particles !!!" << std::endl;

    gentau_found  = 0;
    genTauFromW   = null;
    genTauFromWFromt = null;
    genLeptonFromW   = null;
    genLeptonFromWFromt = null;
    genb1Fromt = null;
    genb2Fromt = null;
    genTauMother  = 0;
    genLeptonMother = 0;
    genb1Mother = 0;
    genb2Mother = 0;
    dR            = null;
    bquark1dR     = null;
    bquark2dR     = null;
    leptondR      = null;
    // W, b, t
    W_pt          = null;
    W_eta         = null;
    W_phi         = null;
    W_energy      = null;
    W_charge      = null;
    // Tau
    gentau_pt     = null;
    gentau_eta    = null;
    gentau_phi    = null;
    gentau_dm     = null;
    gentau_energy = null;
    // lepton1
    genLepton_pt     = null;
    genLepton_eta    = null;
    genLepton_phi    = null;
    genLepton_energy = null;
    genLepton_flavor = null;
    // Gen b-jets
    genb1_flavor = null;
    genb1_pt = null;
    genb1_eta = null;
    genb1_phi = null;
    genb1_energy = null;
    genb2_flavor = null;
    genb2_pt = null;
    genb2_eta = null;
    genb2_phi = null;
    genb2_energy = null;
    // Neutrinos
    nutau_pt      = null;
    nuW_pt        = null;
    nunu_pt       = null;
    nuW_energy    = null;
    nuW_eta       = null;
    nuW_phi       = null;
    nutau_energy  = null;
    nutau_eta     = null;
    nutau_phi     = null;
    nNu           = 0;
    SumNu_pt      = null;
    SumNu_eta     = null;
    SumNu_phi     = null;
    SumNu_energy  = null;
    // PiChar and PiZero
    genPiChar_pt     = null;
    genPiChar_energy = null;
    genPiChar_eta    = null;
    genPiChar_phi    = null;
    genPi0_pt        = null;
    genPi0_energy    = null;
    genPi0_eta       = null;
    genPi0_phi       = null;
    gentau_vis_pt     = null;
    gentau_vis_eta    = null;
    gentau_vis_phi    = null;
    gentau_vis_energy = null;

    if (!isMC) {
        if (monitoring) std::cout << "This is not MC file" << std::endl; 
        return;
    }

    //edm::Handle<pat::PackedCandidateCollection> genParticles;
    edm::Handle<reco::GenParticleCollection> genParticles;
    event.getByToken(GenParticleToken_, genParticles);
    if (!genParticles.isValid()) {
        std::cout << "gen particles collection is not valid" << std::endl;
        return;
    }

    if (monitoringGen) std::cout << "Searching for tau lepton among Gen particles" << std::endl;

    const int pdg_tau      = 15;
    const int pdg_pi0      = 111;
    const int pdg_pi1      = 211;
    const int pdg_rho_plus = 213;
    const int pdg_W        = 24;
    const int pdg_nu_tau   = 16;
    const int pdg_electron = 11;
    const int pdg_nu_ele   = 12;
    const int pdg_mu       = 13;
    const int pdg_nu_mu    = 14;
    const int pdg_bquark   = 5;
    const int pdg_tquark   = 6;

    //const reco::Candidate* tau    = nullptr;
    const reco::GenParticle* tau     = nullptr;
    const reco::GenParticle* lepton  = nullptr;
    const reco::GenParticle* bquark1 = nullptr;
    const reco::GenParticle* bquark2 = nullptr;
    const reco::Candidate* nu_W      = nullptr;
    const reco::Candidate* nu_tau    = nullptr;
    const reco::Candidate* pi0       = nullptr;
    const reco::Candidate* pi1       = nullptr;
    const reco::Candidate* W_lep     = nullptr;
    const reco::Candidate* tquark_tau = nullptr;
    const reco::Candidate* tquark_lep = nullptr;
    const reco::GenParticle* tauMisID = nullptr;
    // test
    //double dRmin = null;
    double dRmin       = 100;
    double leptondRmin = 100;
    double b1dRmin     = 100;
    double b2dRmin     = 100;
    double dRminMisID  = 100;

    math::XYZTLorentzVector TauVisP4;
    math::XYZTLorentzVector SumNuP4;

    std::vector <reco::GenParticle> GenTauCandidates;

    // Look for tau among generated particles
    for (auto& particle: *genParticles) {
#define cut(condition) if (!(condition)) continue;
        // look for the tau -> pi+ pi0 neutrino decay most oriented towards
        // the reconstructed tau (if present)
        cut(abs(particle.pdgId()) == pdg_tau);
        if (monitoringGen) {
            std::cout << "### GenTau Pt  = " << particle.pt() << std::endl;
            std::cout << "### GenTau Eta = " << particle.eta() << std::endl;
            std::cout << "### dR = " << TMath::Sqrt( sqr(dphi(particle.phi(), tau_phi)) + sqr(particle.eta() - tau_eta) ) << std::endl;
        }
        // If tau radiate photon
        if (particle.numberOfDaughters() >= 2 && (abs(particle.daughter(0)->pdgId()) == 15 || abs(particle.daughter(1)->pdgId()) == 15)) continue;

        for (unsigned i = 0; i < particle.numberOfDaughters(); ++i) {
            const reco::Candidate* daughter = particle.daughter(i);
            int id = abs(daughter->pdgId());
            if (id != pdg_nu_tau) {
                TauVisP4 += daughter->p4();
            }
            if (id == pdg_pi0)
                pi0 = daughter;
            else if (id == pdg_pi1)
                pi1 = daughter;
            else if (id == pdg_nu_tau)
                nu_tau = daughter;
        };

        // Do we need these cuts in gen?
        /*
        cut(particle.pt() > tauPtMin);
        cut(TMath::Abs(particle.eta()) < tauEtaMax);
        cut((pv_position - particle.vertex()).R() < tauDzMax);
        */

        double dR_ = null;
        if (tau_found > 0) {
            dR_ = TMath::Sqrt( sqr(dphi(particle.phi(), tau_phi)) + sqr(particle.eta() - tau_eta) );
            if (!tau || dR_ < dRmin) {
                tau = &particle;
                dRmin = dR_;
                if (monitoringGen) {
                    std::cout << "---------- Tau" << std::endl;
                    GenRecoMonitor *TauGenReco = new GenRecoMonitor(particle, pdg_tau, tau_pt, tau_eta, tau_phi);
                    TauGenReco->PrintComp(true, true);
                    delete TauGenReco;
                    std::cout << "decay mode = " << GenTauDecayMode(particle) << std::endl;
                }
            };
        }
#undef cut
        //if (monitoringGen) GenTauMonitor(particle);
    };

    //if (!tau) return;
    if (!tau && monitoringGen) std::cout << "No tau leptons among gen particles" << std::endl;

    if (monitoringGen) std::cout << "Investigate the source of tau (mother)" << std::endl;

    if (tau) {
        gentau_pt     = tau->pt();
        gentau_energy = tau->energy();
        gentau_eta    = tau->eta();
        gentau_phi    = tau->phi();
        gentau_dm     = GenTauDecayMode(*tau);
        dR           = dRmin;
        gentau_found = 1;
        genTauFromW  = 0;
        genTauFromWFromt = 0;
        // Look for W which is mother particle for tau
        for (auto p = tau->mother(); p; p = p->mother()) {
            genTauMother = p->pdgId();
            if (abs(p->pdgId()) == pdg_W) {
                if (monitoringGen) {
                    std::cout << "---------- W-boson (tau mother)" << std::endl;
                    GenRecoMonitor *WGen = new GenRecoMonitor(*p);
                    WGen->PrintGen(true, true);
                    delete WGen;
                }
                genTauFromW = 1;
                W_pt        = p->pt();
                W_eta       = p->eta();
                W_phi       = p->phi();
                W_energy    = p->energy();
                W_charge    = p->charge();
                for (unsigned l = 0; l < p->numberOfDaughters(); l++) {
                    const reco::Candidate* Wdaughter = p->daughter(l);
                    //if (monitoringGen) std::cout << "W_daughter(" << l << ") = " << Wdaughter->pdgId() << std::endl;
                    if (abs(Wdaughter->pdgId()) == pdg_nu_tau) {
                        nu_W = Wdaughter;
                    } else continue;
                }
                for (auto p1 = p->mother(); p1; p1 = p1->mother()) {
                    if (abs(p1->pdgId()) == pdg_tquark) {
                        if (monitoringGen) {
                            std::cout << "---------- t-quark 1" << std::endl;
                            GenRecoMonitor *tGen = new GenRecoMonitor(*p1);
                            tGen->PrintGen(true, true);
                            delete tGen;
                        }
                        genTauFromWFromt = 1;
                        tquark_tau = p1;
                        break;
                    }
                }
                break;
            };
        };
    }

    // Search for secondlepton (mu, ele or tau) and reconstruct the decay chain in reverse direction
    genLeptonFromW = 0;
    genLeptonFromWFromt = 0;
    genb1Fromt = 0;
    genb2Fromt = 0;
    // Look for second lepton among generated particles
    for (auto& particle: *genParticles) {
        //cut(abs(particle.pdgId()) == pdg_mu || abs(particle.pdgId()) == pdg_electron);
        if (abs(particle.pdgId()) == pdg_mu || abs(particle.pdgId()) == pdg_electron || (abs(particle.pdgId()) == pdg_tau && UseTau)) {
            //double lepdR1 = sqrt(deltaR2(particle.eta(), particle.phi(), lepton1_eta, lepton1_phi));
            //double lepdR2 = sqrt(deltaR2(particle.eta(), particle.phi(), lepton2_eta, lepton2_phi));
            //double leptondR_ = TMath::Min(lepdR1, lepdR2);
            double dRtau  = sqrt(deltaR2(particle.eta(), particle.phi(), gentau_eta, gentau_phi));
            double leptondR_ = sqrt(deltaR2(particle.eta(), particle.phi(), lepton1_eta, lepton1_phi));
            if ( (!lepton || leptondR_ < leptondRmin) && dRtau > 0.4) {
                lepton = &particle;
                leptondRmin = leptondR_;
                if (monitoringGen) {
                    std::cout << "---------- Lepton" << std::endl;
                    GenRecoMonitor *Lepton1Gen = new GenRecoMonitor(particle, lepton1_flavor, lepton1_pt, lepton1_eta, lepton1_phi);
                    Lepton1Gen->PrintComp(true, true);
                    delete Lepton1Gen;
                }
            }
        } else continue;
        //if (monitoringGen) GenTauMonitor(particle);
    };

    if (lepton) {
        leptondR = leptondRmin;
        genLepton_pt     = lepton->pt();
        genLepton_eta    = lepton->eta();
        genLepton_phi    = lepton->phi();
        genLepton_energy = lepton->energy();
        genLepton_flavor = lepton->pdgId();
        // Look for W which is mother particle for lepton
        for (auto p = lepton->mother(); p; p = p->mother()) {
            genLeptonMother = p->pdgId();
            //if (monitoringGen) std::cout << "gnetau mother pdg ID = " << genTauMother << std::endl;
            if (abs(p->pdgId()) == pdg_W) {
                if (monitoringGen) {
                    std::cout << "---------- W-boson (lepton mother)" << std::endl;
                    GenRecoMonitor *WGen = new GenRecoMonitor(*p);
                    WGen->PrintGen(true, true);
                    delete WGen;
                }
                genLeptonFromW = 1;
                for (auto p1 = p->mother(); p1; p1 = p1->mother()) {
                    if (abs(p1->pdgId()) == pdg_tquark) {
                        if (monitoringGen) {
                            std::cout << "---------- t-quark 2" << std::endl;
                            GenRecoMonitor *tGen = new GenRecoMonitor(*p1);
                            tGen->PrintGen(true, true);
                            delete tGen;
                        }
                        genLeptonFromWFromt = 1;
                        tquark_lep = p1;
                        //FirstTCharge = p1->charge();
                        break;
                    }
                }
                break;
            } else if (abs(p->pdgId()) != abs(lepton->pdgId())) {
                break;
            }
        };
    }

    // Investigate the gen source of b-quarks
    // if it has been found
    for (auto& particle: *genParticles) {
        if (abs(particle.pdgId()) == pdg_bquark) {
            double bquark1dR_ = sqrt(deltaR2(particle.eta(), particle.phi(), Jet1_eta, Jet1_phi));
            if ( !bquark1 || bquark1dR_ < b1dRmin) {
                bquark1 = &particle;
                b1dRmin = bquark1dR_;
                bquark1dR = b1dRmin;
                if (monitoringGen) {
                    std::cout << "---------- b-quark 1" << std::endl;
                    std::cout << "dR(calc) = " << bquark1dR << std::endl;
                    GenRecoMonitor *b1Gen = new GenRecoMonitor(particle, pdg_bquark, Jet1_pt, Jet1_eta, Jet1_phi);
                    b1Gen->PrintComp(true, true);
                    delete b1Gen;
                }
            }
        } else continue;
    }

    for (auto& particle: *genParticles) {
        if (abs(particle.pdgId()) == pdg_bquark) {
            double bquark2dR_ = 100;
            double BJet2_pt = null;
            double BJet2_eta = null;
            double BJet2_phi = null;
            for (unsigned ijet = 0; ijet < looseBJetsCopy.size(); ijet++) {
                bquark2dR_ = sqrt(deltaR2(particle.eta(), particle.phi(), looseBJetsCopy[ijet].Eta, looseBJetsCopy[ijet].Phi));
                if (!bquark2 || bquark2dR_ < b2dRmin) {
                    bquark2 = &particle;
                    b2dRmin = bquark2dR_;
                    BJet2_pt = looseBJetsCopy[ijet].Pt;
                    BJet2_eta = looseBJetsCopy[ijet].Eta;
                    BJet2_phi = looseBJetsCopy[ijet].Phi;
                }
                //looseBJets[ijet].Monitoring();
            }
            bquark2dR = b2dRmin;
            if (monitoringGen) {
                std::cout << "---------- b-quark 2" << std::endl;
                std::cout << "dR(calc) = " << bquark2dR << std::endl;
                GenRecoMonitor *b2Gen = new GenRecoMonitor(particle, pdg_bquark, BJet2_pt, BJet2_eta, BJet2_phi);
                b2Gen->PrintComp(true, true);
                delete b2Gen;
            }
            //
            /*
            if (BJet2_pt > 0) {
                bquark2dR_ = sqrt(deltaR2(particle.eta(), particle.phi(), BJet2_eta, BJet2_phi));
            } else if ((Jet1_bprob + Jet1_bbprob) >= (Jet2_bprob + Jet2_bbprob)) {
                bquark2dR_ = sqrt(deltaR2(particle.eta(), particle.phi(), Jet1_eta, Jet1_phi));
            } else {
                bquark2dR_ = sqrt(deltaR2(particle.eta(), particle.phi(), Jet2_eta, Jet2_phi));
            }
            if ( !bquark2 || bquark2dR_ < b2dRmin) {
                bquark2 = &particle;
                b2dRmin = bquark2dR_;
                bquark2dR = b2dRmin;
                if (monitoringGen) {
                    std::cout << "---------- b-quark 2" << std::endl;
                    std::cout << "dR(calc) = " << bquark2dR << std::endl;
                    GenRecoMonitor *b2Gen = new GenRecoMonitor(particle, pdg_bquark, BJet2_pt, BJet2_eta, BJet2_phi);
                    b2Gen->PrintComp(true, true);
                    delete b2Gen;
                }
            }
            */
        } else continue;
    }
    looseBJetsCopy.clear();

    if (bquark1) {
        genb1_flavor = bquark1->pdgId();
        genb1_pt = bquark1->pt();
        genb1_eta = bquark1->eta();
        genb1_phi = bquark1->phi();
        genb1_energy = bquark1->energy();
        for (auto p = bquark1->mother(); p; p = p->mother()) {
            genb1Mother = p->pdgId();
            //if (monitoringGen) std::cout << "gnetau mother pdg ID = " << genTauMother << std::endl;
            if (abs(p->pdgId()) == pdg_tquark) {
                double t1dR = abs(null) + 1;
                double t2dR = abs(null) + 1;
                if (tquark_tau) {
                    t1dR = sqrt(deltaR2(p->eta(), p->phi(), tquark_tau->eta(), tquark_tau->phi()));
                }
                if (tquark_lep) {
                    t2dR = sqrt(deltaR2(p->eta(), p->phi(), tquark_lep->eta(), tquark_lep->phi()));
                }
                if (t1dR < abs(null) || t2dR < abs(null)) {
                    if (t1dR < t2dR) genb1Fromt = 1;
                    else genb1Fromt = 2;
                }
                break;
            } else if (abs(p->pdgId()) != abs(bquark1->pdgId())) {
                break;
            }
        }
    }

    if (bquark2) {
        genb2_flavor = bquark2->pdgId();
        genb2_pt = bquark2->pt();
        genb2_eta = bquark2->eta();
        genb2_phi = bquark2->phi();
        genb2_energy = bquark2->energy();
        for (auto p = bquark2->mother(); p; p = p->mother()) {
            genb2Mother = p->pdgId();
            if (abs(p->pdgId()) == pdg_tquark) {
                double t1dR = abs(null) + 1;
                double t2dR = abs(null) + 1;
                if (tquark_tau) {
                    t1dR = sqrt(deltaR2(p->eta(), p->phi(), tquark_tau->eta(), tquark_tau->phi()));
                }
                if (tquark_lep) {
                    t2dR = sqrt(deltaR2(p->eta(), p->phi(), tquark_lep->eta(), tquark_lep->phi()));
                }
                if (t1dR < abs(null) || t2dR < abs(null)) {
                    if (t1dR < t2dR) genb2Fromt = 1;
                    else genb2Fromt = 2;
                }
                break;
            } else if (abs(p->pdgId()) != abs(bquark2->pdgId())) {
                break;
            }
        }
    }

    for (auto& particle: *genParticles) {
        if ( (abs(particle.pdgId()) == pdg_nu_tau || abs(particle.pdgId()) == pdg_nu_mu || abs(particle.pdgId()) == pdg_nu_ele)
            && (particle.pdgId() != particle.mother()->pdgId()) ) {
            SumNuP4 += particle.p4();
            nNu++;
        }
    }

    if (monitoringGen) std::cout << "gentau Flag 3" << std::endl;

    if (pi1) {
        genPiChar_pt      = pi1->pt();
        genPiChar_energy  = pi1->energy();
        genPiChar_eta     = pi1->eta();
        genPiChar_phi     = pi1->phi();
    } else {
        genPiChar_pt     = null;
        genPiChar_energy = null;
        genPiChar_eta    = null;
        genPiChar_phi    = null;
    }

    if (pi0) {
        genPi0_pt      = pi0->pt();
        genPi0_energy  = pi0->energy();
        genPi0_eta     = pi0->eta();
        genPi0_phi     = pi0->phi();
    } else {
        genPi0_pt     = null;
        genPi0_energy = null;
        genPi0_eta    = null;
        genPi0_phi    = null;
    }

    if (tau) {
        gentau_vis_pt = TauVisP4.pt();
        gentau_vis_eta = TauVisP4.eta();
        gentau_vis_phi = TauVisP4.phi();
        gentau_vis_energy = TauVisP4.energy();
    } else {
        gentau_vis_pt = null;
        gentau_vis_eta = null;
        gentau_vis_phi = null;
        gentau_vis_energy = null;
    }

    /*
    if (lepton) {
        genLepton_flavor = lepton->pdgId();
        genLepton_mother = lepton->mother()->pdgId();
        genLepton_eta    = lepton->eta();
        genLepton_phi    = lepton->phi();
        genLepton_pt     = lepton->pt();
        genLepton_energy = lepton->energy();
    } else {
        genLepton_flavor = null;
        genLepton_mother = null;
        genLepton_eta    = null;
        genLepton_phi    = null;
        genLepton_pt     = null;
        genLepton_energy = null;
    }
    */

    if (nNu > 0) {
        SumNu_pt = SumNuP4.pt();
        SumNu_eta = SumNuP4.eta();
        SumNu_phi = SumNuP4.phi();
        SumNu_energy = SumNuP4.energy();
    } else {
        SumNu_pt = null;
        SumNu_eta = null;
        SumNu_phi = null;
        SumNu_energy = null;
    }

    if (monitoringGen) std::cout << "gentau Flag 4" << std::endl;

    //if (!nu_W) return; 

    if (monitoringGen) std::cout << "gentau Flag 5" << std::endl;

    if (nu_tau) {
        nutau_pt     = nu_tau->pt();
        nutau_energy = nu_tau->energy();
        nutau_eta    = nu_tau->eta();
        nutau_phi    = nu_tau->phi();
    }
    if (genTauFromW > 0 && nu_W) {
        nuW_pt       = nu_W->pt();
        nuW_energy   = nu_W->energy();
        nuW_eta      = nu_W->eta();
        nuW_phi      = nu_W->phi();
    }
    if (nu_W && nu_tau)  {
        double nutau_m = nu_tau->mass();
        double nuW_m = nu_W->mass();
        TLorentzVector gennutau, gennuW;
        gennutau.SetPtEtaPhiM(nutau_pt, nutau_eta, nutau_phi, nutau_m);
        gennuW.SetPtEtaPhiM(nuW_pt, nuW_eta, nuW_phi, nuW_m);
        nunu_pt = (gennutau + gennuW).Pt();
    }

    for (auto& particle: *genParticles) {
        double dR_ = null;
        if (tau_found > 0) {
            dR_ = TMath::Sqrt( sqr(dphi(particle.phi(), tau_phi)) + sqr(particle.eta() - tau_eta) );
            if (!tauMisID || dR_ < dRminMisID) {
                tauMisID = &particle;
                dRminMisID = dR_;
                /*
                if (monitoringGen) {
                    std::cout << "---------- Tau MisID" << std::endl;
                    GenRecoMonitor *TauGenReco = new GenRecoMonitor(particle, pdg_tau, tau_pt, tau_eta, tau_phi);
                    TauGenReco->PrintComp(true, true);
                    delete TauGenReco;
                    //std::cout << "decay mode = " << GenTauDecayMode(particle) << std::endl;
                }
                */
            };
        }
    }
    if (dRminMisID < 0.1) genTauMisID = tauMisID->pdgId();
    else genTauMisID = null;
};


bool TTbarTauLepton::AddMET(const edm::Event& event) {
    edm::Handle<pat::METCollection> mets;
    event.getByToken(MetCollectionToken_, mets);
    if (!mets.isValid() || !mets->size()) return false;

    auto& MET = mets->front();
    met              = MET.pt();
    met_phi          = MET.phi();
    met_eta          = MET.eta();
    met_energy       = MET.energy();
    met_significance = MET.significance();
    met_mEtSig       = MET.mEtSig();

    TLorentzVector METp4, taup4;
    METp4.SetPtEtaPhiE(met, met_eta, met_phi, met_energy);
    taup4.SetPtEtaPhiM(tau_pt, tau_eta, tau_phi, tau_m);
    tauMET_mass = (METp4 + taup4).M();

    if (monitoringMET) {
        std::cout << "---- MET -----------" << std::endl;
        METMonitor(MET);
    }
    
    edm::Handle<pat::METCollection> Puppimets;
    event.getByToken(PuppiMetCollectionToken_, Puppimets);
    //if (!mets.isValid() || !mets->size()) return false;
    auto& PuppiMET = Puppimets->front();
    Puppimet              = PuppiMET.pt();
    Puppimet_phi          = PuppiMET.phi();
    Puppimet_eta          = PuppiMET.eta();
    Puppimet_energy       = PuppiMET.energy();
    Puppimet_significance = PuppiMET.significance();
    Puppimet_metSig       = PuppiMET.mEtSig(); 

    TLorentzVector PuppiMETp4;
    PuppiMETp4.SetPtEtaPhiE(Puppimet, Puppimet_eta, Puppimet_phi, Puppimet_energy);
    tauPuppimet_mass = (PuppiMETp4 + taup4).M();

    if (monitoringMET) {
        std::cout << "---- MET Pippi -----" << std::endl;
        METMonitor(PuppiMET);
    }

    return true;
}

bool TTbarTauLepton::AddVertex(const edm::Event& event) {
    edm::Handle<reco::VertexCollection> vertices;
    event.getByToken(PVToken_, vertices);
    if (!vertices.isValid()) return false;

    nVtx = vertices->size();
    if (nVtx == 0) return false;
    pv_position = vertices->front().position();
    return true;
};

// Jets analysis

bool TTbarTauLepton::JetPtSum (const edm::Event& event) {

    if (monitoringBJets) std::cout << std::endl << "!!! Jets !!!" << std::endl;

    // Puppi Jets
    PuppijetPtSum30 = 0;
    PuppijetPtSum20 = 0;
    nPuppiJets30    = 0;
    nPuppiJets20    = 0;
    PuppijetPtSum30PV = 0;
    PuppijetPtSum20PV = 0;
    nPuppiJets30PV    = 0;
    nPuppiJets20PV    = 0;
    nLooseBtagedPuppiJets    = 0;
    nMediumBtagedPuppiJets   = 0;
    nTightBtagedPuppiJets    = 0;
    nLooseBtagedPuppiJetsPV  = 0;
    nMediumBtagedPuppiJetsPV = 0;
    nTightBtagedPuppiJetsPV  = 0;

    Jet1_pt   = null;
    Jet1_eta  = null;
    Jet1_phi  = null;
    Jet1_m    = null;
    Jet1_E = null;
    Jet1_bprob = null;
    Jet1_bbprob = null;
    Jet1_hadronFlavour = null;

    Jet2_pt   = null;
    Jet2_eta  = null;
    Jet2_phi  = null;
    Jet2_m    = null;
    Jet2_E = null;
    Jet2_bprob = null;
    Jet2_bbprob = null;
    Jet2_hadronFlavour = null;

    Jet3_pt   = null;
    Jet3_eta  = null;
    Jet3_phi  = null;
    Jet3_m    = null;
    Jet3_E = null;
    Jet3_bprob = null;
    Jet3_bbprob = null;
    Jet3_hadronFlavour = null;

    Jet4_pt   = null;
    Jet4_eta  = null;
    Jet4_phi  = null;
    Jet4_m    = null;
    Jet4_E = null;
    Jet4_bprob = null;
    Jet4_bbprob = null;
    Jet4_hadronFlavour = null;

    Jet5_pt   = null;
    Jet5_eta  = null;
    Jet5_phi  = null;
    Jet5_m    = null;
    Jet5_E = null;
    Jet5_bprob = null;
    Jet5_bbprob = null;
    Jet5_hadronFlavour = null;

    Jet6_pt   = null;
    Jet6_eta  = null;
    Jet6_phi  = null;
    Jet6_m    = null;
    Jet6_E    = null;
    Jet6_bprob = null;
    Jet6_bbprob = null;
    Jet6_hadronFlavour = null;

    lepton1Jet_pt = null;
    lepton1Jet_eta = null;
    lepton1Jet_phi = null;
    lepton1Jet_E = null;
    lepton1Jet_bprob = null;

    TauJet_pt = null;
    TauJet_eta = null;
    TauJet_phi = null;
    TauJet_E = null;
    TauJet_bprob = null;
    TauJet_hadronFlavour = null;

    double deltaRTauJetMin = 100.;
    double deltaRLepton1JetMin = 100.;

    int nLooseBJets = 0;

    // Btag calibration
    //BTagCalibration calib("csvv1", "/afs/cern.ch/user/a/aoskin/Tau_packeges/CMSSW_9_4_6/src/Tau/TreeMakerMiniAOD/test/TTBarTauLepton_only/DeepCSV_94XSF_V5_B_F.csv");
    // operating point, central sys type, other sys types
    
    // Tracks and vertices
    edm::Handle<pat::IsolatedTrackCollection> tracks;
    event.getByToken(TrackToken_, tracks);
    if (!tracks.isValid()) return false;
    edm::Handle<reco::VertexCollection> vertices;
    event.getByToken(PVToken_, vertices);
    if (!vertices.isValid() || vertices->size() == 0) return false;

    // Working points
    double WPBTag_loose   = 0.1522;
    double WPBTag_medium  = 0.4941;
    double WPBTag_tight   = 0.8001;
    std::vector <BJetCandidate> looseBJets;

    // Selection of two bJets for ttbar analysis
    // Also study of multijets and leading/subleading jets in event

    edm::Handle<pat::JetCollection> Puppijets;
    event.getByToken(PuppiJetCollectionToken_, Puppijets);
    if (!Puppijets.isValid()) return false;

    int j = 0;
    // Lopp over Jets
    for (auto& jet: *Puppijets) {
        if (TMath::Abs(jet.eta()) > JetEtaMax) continue;
        ///if (TMath::Abs(jet.eta()) > 3.) continue;
        if (jet.pt() > 20) {
            j++;
            if (monitoringJets) {
                std::cout << std::endl;
                std::cout << "PuppiJet[" << j << "] Pt = " << jet.pt() << std::endl;
                std::cout << "(eta, phi) = (" << jet.eta() << ", " << jet.phi() << ")" << std::endl;
                std::cout << "hadronFlavour = " << jet.hadronFlavour() << std::endl; 
                std::cout << "partonFlavour = " << jet.partonFlavour() << std::endl; 
                std::cout << "b_prob  = " << jet.bDiscriminator("pfDeepCSVJetTags:probb") << std::endl;
                std::cout << "bb_prob = " << jet.bDiscriminator("pfDeepCSVJetTags:probbb") << std::endl;
                std::cout << "neutralHadronEnergyFraction = " << jet.neutralHadronEnergy() / jet.energy() << std::endl; 
                std::cout << "neutralEmEnergyFraction     = " << jet.neutralEmEnergyFraction() << std::endl; 
                std::cout << "muonEnergyFraction          = " << jet.muonEnergyFraction() << std::endl; 
                std::cout << "chargedHadronEnergyFraction = " << jet.chargedHadronEnergyFraction() << std::endl; 
                std::cout << "chargedEmEnergyFraction     = " << jet.chargedEmEnergyFraction() << std::endl; 
                std::cout << "chargedHadronMultiplicity   = " << jet.chargedHadronMultiplicity() << std::endl;
                std::cout << "chargedMultiplicity         = " << jet.chargedMultiplicity() << std::endl;
                std::cout << "neutralHadronMultiplicity   = " << jet.neutralHadronMultiplicity() << std::endl; 
                std::cout << "neutralMultiplicity         = " << jet.neutralMultiplicity() << std::endl;
                std::cout << "pileup = " << jet.pileup() << std::endl;
            }
            // Search for tau jet
            double deltaRTauJet = sqrt(deltaR2(tau_eta, tau_phi, jet.eta(), jet.phi()));
            if (deltaRTauJet < 0.4 && deltaRTauJet < deltaRTauJetMin) {
                deltaRTauJetMin = deltaRTauJet;
                if (monitoringJets) std::cout << "This is Tau jet, dR = " << deltaRTauJet << std::endl;
                nTauJetInEvent++;
                TauJet_pt = jet.pt();
                TauJet_eta = jet.eta();
                TauJet_phi = jet.phi();
                TauJet_E = jet.energy();
                TauJet_hadronFlavour = jet.hadronFlavour();
                TauJet_bprob = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb");
                continue; 
            }
            // Exclude jets close to lepton 1
            if (lepton1_pt > 0) {
                double deltaRLepton1 = sqrt(deltaR2(jet.eta(), jet.phi(), lepton1_eta, lepton1_phi));
                //deltaRLepton1 = sqrt(deltaR2(jet.eta(), jet.phi(), lepton1_electron_eta, lepton1_electron_phi));
                //deltaRLepton1 = sqrt(deltaR2(jet.eta(), jet.phi(), lepton1_muon_eta, lepton1_muon_phi));
                if (deltaRLepton1 < 0.4 && deltaRLepton1 < deltaRLepton1JetMin) {
                    deltaRLepton1JetMin = deltaRLepton1;
                    if (monitoringJets) {
                        std::cout << "Jet with lepton inside cone dR < 0.4" << std::endl;
                        std::cout << "dR = " << deltaRLepton1 << std::endl;
                    }
                    lepton1Jet_pt = jet.pt();
                    lepton1Jet_eta = jet.eta();
                    lepton1Jet_phi = jet.phi();
                    lepton1Jet_E = jet.energy();
                    lepton1Jet_bprob = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb");
                    continue; ///
                }
            }
            reco::TrackRefVector JetTracksRef = jet.associatedTracks();
            std::map <int, int> VertexTracksMap;
            std::map <int, std::pair<int, double>> VertexTracksPtMap;
            int NtracksFromPV = 0;
            bool JetFromPV = true;
            std::vector<int> TrackCounter;
            std::vector<double> TracksPtSum;
            for (unsigned count = 0; count < vertices->size(); count++) {
                TrackCounter.push_back(0);
                TracksPtSum.push_back(0);
            }
            // Loop over tracks from Jet
            for (reco::TrackRefVector::const_iterator i_trk = JetTracksRef.begin(); i_trk != JetTracksRef.end(); i_trk++) {
                // Loop over vertices
                for (unsigned nvtx = 0; nvtx < vertices->size(); nvtx++) {
                    // number of tracks from this vertex in jet
                    VertexTracksMap.insert(std::pair<int, int> (nvtx, TrackCounter.at(nvtx)));
                    std::pair <int, double> Pair1 (TrackCounter.at(nvtx), TracksPtSum.at(nvtx));
                    VertexTracksPtMap.insert(std::pair<int, std::pair<int, double>> (nvtx, Pair1));
                    reco::Vertex::trackRef_iterator it = (*vertices)[nvtx].tracks_begin();
                    reco::Vertex::trackRef_iterator lastTrack = (*vertices)[nvtx].tracks_end();
                    int ntrk = 0;
                    // Loop over tracks from vertex
                    for (; it != lastTrack; ++it, ++ntrk) {
                        if (abs((*i_trk)->eta() - (*it)->eta()) < 0.0001 && abs((*i_trk)->phi() - (*it)->phi()) < 0.0001) {
                            TrackCounter.at(nvtx)++;
                            TracksPtSum.at(nvtx) += (*i_trk)->pt();
                            std::pair <int, double> NTraksSumPt (TrackCounter.at(nvtx), TracksPtSum.at(nvtx));
                            VertexTracksMap.at(nvtx) = TrackCounter.at(nvtx);
                            VertexTracksPtMap.at(nvtx) = NTraksSumPt;
                        }
                    }
                }
            }
            if (monitoringJets) std::cout << std::endl;
            for (auto Map_iter = VertexTracksMap.begin(); Map_iter != VertexTracksMap.end(); ++Map_iter) {
                if ((*Map_iter).first == 0) NtracksFromPV = (*Map_iter).second;
                if ((*Map_iter).first > 0 && (*Map_iter).second > NtracksFromPV) JetFromPV = false;
                if ((*Map_iter).second < 1) continue;
            }
            for (auto Map_iter1 = VertexTracksPtMap.begin(); Map_iter1 != VertexTracksPtMap.end(); ++Map_iter1) {
                if ((*Map_iter1).second.first < 1) continue;
                if (monitoringJets) std::cout << "Vertex number   = " << (*Map_iter1).first << "  entries = " << (*Map_iter1).second.first << "  SumPt = " << (*Map_iter1).second.second << std::endl;
            }

            if ((jet.pt() > BJetPtMin) && (TMath::Abs(jet.eta()) < JetEtaMax)) {
                // Jets selection criteria
                if ( (jet.neutralHadronEnergy() / jet.energy() < 0.90) &&
                    (jet.neutralEmEnergyFraction() < 0.90) &&
                    (jet.muonEnergyFraction() < 0.80) &&
                    (jet.chargedHadronEnergyFraction() > 0.) &&
                    (jet.chargedMultiplicity() + jet.neutralMultiplicity() > 0.) &&
                    (jet.chargedEmEnergyFraction() < 0.80)) {
                    // Fill vector of jets
                    nLooseBJets++;
                    BJetCandidate BJet(jet, JetFromPV);
                    looseBJets.push_back(BJet);
                }
            }
            
            PuppijetPtSum20 += jet.pt();
            nPuppiJets20++;
            if (monitoringJets) std::cout << "Jet added to 20 GeV" << std::endl;
            if (JetFromPV) {
                PuppijetPtSum20PV += jet.pt();
                nPuppiJets20PV++;
                if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_loose) {
                    nLooseBtagedPuppiJetsPV++;
                    if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_medium) {
                        nMediumBtagedPuppiJetsPV++;
                        if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_tight) {
                            nTightBtagedPuppiJetsPV++;
                        }
                    }
                }
            }
            if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_loose) {
                nLooseBtagedPuppiJets++;
                if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_medium) {
                    nMediumBtagedPuppiJets++;
                    if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_tight) {
                        nTightBtagedPuppiJets++;
                    }
                }
            }

            if (jet.pt() > 30) {
                PuppijetPtSum30 += jet.pt();
                nPuppiJets30++;
                if (JetFromPV) {
                    nPuppiJets30PV++;
                    PuppijetPtSum30PV += jet.pt();
                    /*
                    if (jet.neutralHadronEnergy() / jet.energy() >= 0.99) continue; // 0.9
                    if (jet.neutralEmEnergyFraction() >= 0.99) continue; // 0.9
                    //if (jet.getPFConstituents().size() <= 1) continue; // why it was commented?
                    if (jet.muonEnergyFraction() > 0.8) continue;
                    if (jet.chargedHadronEnergyFraction() <= 0.001) continue;
                    if (jet.chargedMultiplicity() <= 0) continue;
                    if (jet.chargedEmEnergyFraction() >= 0.99) continue; // 0.8
                    */
                }
            }
        }
    }
    if (monitoringBJets) std::cout << "number of 30 GeV Jets from PV (passed selections) = " << nPuppiJets30PV
    << "(" << nPuppiJets30PV << ")" << std::endl;

    if (monitoringBJets) std::cout << std::endl << "!!! bTagged Jets !!!" << std::endl;

    SortJets(looseBJets);
    looseBJetsCopy = looseBJets;

    if (monitoringBJets) {
        std::cout << "List of loose bJets (" << looseBJets.size() << ", " << nLooseBJets << "):" << std::endl;
        for (unsigned ijet = 0; ijet < looseBJets.size(); ijet++) {
            std::cout << "bjet " << ijet << std::endl;
            looseBJets[ijet].Monitoring();
        }
    }

    if (looseBJets.size() > 5) {
        Jet1_pt     = looseBJets[0].Pt;
        Jet1_eta    = looseBJets[0].Eta;
        Jet1_phi    = looseBJets[0].Phi;
        Jet1_bprob  = looseBJets[0].Bprob;
        Jet1_bbprob = looseBJets[0].BBprob;
        Jet1_E      = looseBJets[0].FourMomentum.E();
        Jet1_hadronFlavour = looseBJets[0].HadronFlavour;

        Jet2_pt     = looseBJets[1].Pt;
        Jet2_eta    = looseBJets[1].Eta;
        Jet2_phi    = looseBJets[1].Phi;
        Jet2_bprob  = looseBJets[1].Bprob;
        Jet2_bbprob = looseBJets[1].BBprob;
        Jet2_E      = looseBJets[1].FourMomentum.E();
        Jet2_hadronFlavour = looseBJets[1].HadronFlavour;

        Jet3_pt     = looseBJets[2].Pt;
        Jet3_eta    = looseBJets[2].Eta;
        Jet3_phi    = looseBJets[2].Phi;
        Jet3_bprob  = looseBJets[2].Bprob;
        Jet3_bbprob = looseBJets[2].BBprob;
        Jet3_E      = looseBJets[2].FourMomentum.E();
        Jet3_hadronFlavour = looseBJets[2].HadronFlavour;

        Jet4_pt     = looseBJets[3].Pt;
        Jet4_eta    = looseBJets[3].Eta;
        Jet4_phi    = looseBJets[3].Phi;
        Jet4_bprob  = looseBJets[3].Bprob;
        Jet4_bbprob = looseBJets[3].BBprob;
        Jet4_E      = looseBJets[3].FourMomentum.E();
        Jet4_hadronFlavour = looseBJets[3].HadronFlavour;

        Jet5_pt     = looseBJets[4].Pt;
        Jet5_eta    = looseBJets[4].Eta;
        Jet5_phi    = looseBJets[4].Phi;
        Jet5_bprob  = looseBJets[4].Bprob;
        Jet5_bbprob = looseBJets[4].BBprob;
        Jet5_E      = looseBJets[4].FourMomentum.E();
        Jet5_hadronFlavour = looseBJets[4].HadronFlavour;

        Jet6_pt     = looseBJets[5].Pt;
        Jet6_eta    = looseBJets[5].Eta;
        Jet6_phi    = looseBJets[5].Phi;
        Jet6_bprob  = looseBJets[5].Bprob;
        Jet6_bbprob = looseBJets[5].BBprob;
        Jet6_E      = looseBJets[5].FourMomentum.E();
        Jet6_hadronFlavour = looseBJets[5].HadronFlavour;
    } else if (looseBJets.size() > 4) {
        Jet1_pt     = looseBJets[0].Pt;
        Jet1_eta    = looseBJets[0].Eta;
        Jet1_phi    = looseBJets[0].Phi;
        Jet1_bprob  = looseBJets[0].Bprob;
        Jet1_bbprob = looseBJets[0].BBprob;
        Jet1_E      = looseBJets[0].FourMomentum.E();
        Jet1_hadronFlavour = looseBJets[0].HadronFlavour;

        Jet2_pt     = looseBJets[1].Pt;
        Jet2_eta    = looseBJets[1].Eta;
        Jet2_phi    = looseBJets[1].Phi;
        Jet2_bprob  = looseBJets[1].Bprob;
        Jet2_bbprob = looseBJets[1].BBprob;
        Jet2_E      = looseBJets[1].FourMomentum.E();
        Jet2_hadronFlavour = looseBJets[1].HadronFlavour;

        Jet3_pt     = looseBJets[2].Pt;
        Jet3_eta    = looseBJets[2].Eta;
        Jet3_phi    = looseBJets[2].Phi;
        Jet3_bprob  = looseBJets[2].Bprob;
        Jet3_bbprob = looseBJets[2].BBprob;
        Jet3_E      = looseBJets[2].FourMomentum.E();
        Jet3_hadronFlavour = looseBJets[2].HadronFlavour;

        Jet4_pt     = looseBJets[3].Pt;
        Jet4_eta    = looseBJets[3].Eta;
        Jet4_phi    = looseBJets[3].Phi;
        Jet4_bprob  = looseBJets[3].Bprob;
        Jet4_bbprob = looseBJets[3].BBprob;
        Jet4_E      = looseBJets[3].FourMomentum.E();
        Jet4_hadronFlavour = looseBJets[3].HadronFlavour;

        Jet5_pt     = looseBJets[4].Pt;
        Jet5_eta    = looseBJets[4].Eta;
        Jet5_phi    = looseBJets[4].Phi;
        Jet5_bprob  = looseBJets[4].Bprob;
        Jet5_bbprob = looseBJets[4].BBprob;
        Jet5_E      = looseBJets[4].FourMomentum.E();
        Jet5_hadronFlavour = looseBJets[4].HadronFlavour;
    } else if (looseBJets.size() > 3) {
        Jet1_pt     = looseBJets[0].Pt;
        Jet1_eta    = looseBJets[0].Eta;
        Jet1_phi    = looseBJets[0].Phi;
        Jet1_bprob  = looseBJets[0].Bprob;
        Jet1_bbprob = looseBJets[0].BBprob;
        Jet1_E      = looseBJets[0].FourMomentum.E();
        Jet1_hadronFlavour = looseBJets[0].HadronFlavour;

        Jet2_pt     = looseBJets[1].Pt;
        Jet2_eta    = looseBJets[1].Eta;
        Jet2_phi    = looseBJets[1].Phi;
        Jet2_bprob  = looseBJets[1].Bprob;
        Jet2_bbprob = looseBJets[1].BBprob;
        Jet2_E      = looseBJets[1].FourMomentum.E();
        Jet2_hadronFlavour = looseBJets[1].HadronFlavour;

        Jet3_pt     = looseBJets[2].Pt;
        Jet3_eta    = looseBJets[2].Eta;
        Jet3_phi    = looseBJets[2].Phi;
        Jet3_bprob  = looseBJets[2].Bprob;
        Jet3_bbprob = looseBJets[2].BBprob;
        Jet3_E      = looseBJets[2].FourMomentum.E();
        Jet3_hadronFlavour = looseBJets[2].HadronFlavour;

        Jet4_pt     = looseBJets[3].Pt;
        Jet4_eta    = looseBJets[3].Eta;
        Jet4_phi    = looseBJets[3].Phi;
        Jet4_bprob  = looseBJets[3].Bprob;
        Jet4_bbprob = looseBJets[3].BBprob;
        Jet4_E      = looseBJets[3].FourMomentum.E();
        Jet4_hadronFlavour = looseBJets[3].HadronFlavour;
    } else if (looseBJets.size() > 2) {
        Jet1_pt     = looseBJets[0].Pt;
        Jet1_eta    = looseBJets[0].Eta;
        Jet1_phi    = looseBJets[0].Phi;
        Jet1_bprob  = looseBJets[0].Bprob;
        Jet1_bbprob = looseBJets[0].BBprob;
        Jet1_E      = looseBJets[0].FourMomentum.E();
        Jet1_hadronFlavour = looseBJets[0].HadronFlavour;

        Jet2_pt     = looseBJets[1].Pt;
        Jet2_eta    = looseBJets[1].Eta;
        Jet2_phi    = looseBJets[1].Phi;
        Jet2_bprob  = looseBJets[1].Bprob;
        Jet2_bbprob = looseBJets[1].BBprob;
        Jet2_E      = looseBJets[1].FourMomentum.E();
        Jet2_hadronFlavour = looseBJets[1].HadronFlavour;

        Jet3_pt     = looseBJets[2].Pt;
        Jet3_eta    = looseBJets[2].Eta;
        Jet3_phi    = looseBJets[2].Phi;
        Jet3_bprob  = looseBJets[2].Bprob;
        Jet3_bbprob = looseBJets[2].BBprob;
        Jet3_E      = looseBJets[2].FourMomentum.E();
        Jet3_hadronFlavour = looseBJets[2].HadronFlavour;
    } else if (looseBJets.size() > 1) {
        Jet1_pt     = looseBJets[0].Pt;
        Jet1_eta    = looseBJets[0].Eta;
        Jet1_phi    = looseBJets[0].Phi;
        Jet1_bprob  = looseBJets[0].Bprob;
        Jet1_bbprob = looseBJets[0].BBprob;
        Jet1_E      = looseBJets[0].FourMomentum.E();
        Jet1_hadronFlavour = looseBJets[0].HadronFlavour;

        Jet2_pt     = looseBJets[1].Pt;
        Jet2_eta    = looseBJets[1].Eta;
        Jet2_phi    = looseBJets[1].Phi;
        Jet2_bprob  = looseBJets[1].Bprob;
        Jet2_bbprob = looseBJets[1].BBprob;
        Jet2_E      = looseBJets[1].FourMomentum.E();
        Jet2_hadronFlavour = looseBJets[1].HadronFlavour;
    } else if (looseBJets.size() > 0) {
        Jet1_pt     = looseBJets[0].Pt;
        Jet1_eta    = looseBJets[0].Eta;
        Jet1_phi    = looseBJets[0].Phi;
        Jet1_bprob  = looseBJets[0].Bprob;
        Jet1_bbprob = looseBJets[0].BBprob;
        Jet1_E      = looseBJets[0].FourMomentum.E();
        Jet1_hadronFlavour = looseBJets[0].HadronFlavour;
    }

    if (monitoringJets) {
        std::cout << "b-Jets found " << looseBJets.size() << std::endl;
    }

    looseBJets.clear();
    return true;
};

bool TTbarTauLepton::AddLepton (const edm::Event& event) {

    if (monitoringLeptons) std::cout << std::endl << "!!! Leptons !!!" << std::endl;

    nLeptonCandidates = 0;
    VetoLeptons = 0;
    VetoTaus = 0;
    bool LeptonFound = false;

    lepton1_pt       = null;
    lepton1_eta      = null;
    lepton1_phi      = null;
    lepton1_dz       = null;
    lepton1_flavor   = null;
    lepton1_charge   = null;
    lepton1_E        = null;
    lepton1_trackIso = null;
    lepton1_sumPuppiIso = null;
    lepton1_sumPuppiNoLeptonIso = null;

    lepton1_electron_cutBasedID_loose    = null;
    lepton1_electron_cutBasedID_medium   = null;
    lepton1_electron_cutBasedID_tight    = null;
    lepton1_electron_mvaNoIsoID_loose    = null;
    lepton1_electron_mvaNoIsoID_wp80     = null;
    lepton1_electron_mvaNoIsoID_wp90     = null;

    lepton1_muon_CutBasedIdLoose = null;
    lepton1_muon_CutBasedIdMedium = null;
    lepton1_muon_CutBasedIdTight = null;
    lepton1_muon_CutBasedIdGlobalHighPt = null;
    lepton1_muon_CutBasedIdTrkHighPt = null;
    lepton1_muon_PFIsoLoose = null;
    lepton1_muon_PFIsoMedium = null;
    lepton1_muon_PFIsoTight = null;
    lepton1_muon_TkIsoLoose = null;
    lepton1_muon_TkIsoTight = null;
    lepton1_muon_MvaLoose = null;
    lepton1_muon_MvaMedium = null;
    lepton1_muon_MvaTight = null;

    lepton1_tau_dm  = null;
    lepton1_tau_m   = null;
    lepton1_tau_absIso           = null;
    lepton1_tau_looseCombinedIso = null;
    lepton1_tau_tightCombinedIso = null;
    lepton1_tau_looseMvaIso      = null;
    lepton1_tau_mediumMvaIso     = null;
    lepton1_tau_tightMvaIso      = null;
    lepton1_tau_VtightMvaIso     = null;
    lepton1_tau_VVtightMvaIso    = null;
    // New Raw discriminators
    lepton1_tau_againstElectronRaw            = null;
    lepton1_tau_IsoMVArun2v1DBdR03oldDMwLTraw = null;
    //
    lepton1_tau_tightMuonRejection      = null;
    lepton1_tau_looseElectronRejection  = null;
    lepton1_tau_mediumElectronRejection = null;
    lepton1_tau_tightElectronRejection  = null;
    lepton1_tau_VtightElectronRejection = null;

    // lepton 2 parameters

    lepton2_pt       = null;
    lepton2_eta      = null;
    lepton2_phi      = null;
    lepton2_dz       = null;
    lepton2_flavor   = null;
    lepton2_charge   = null;
    lepton2_E        = null;
    lepton2_trackIso = null;
    lepton2_sumPuppiIso = null;
    lepton2_sumPuppiNoLeptonIso = null;

    lepton2_electron_cutBasedID_loose    = null;
    lepton2_electron_cutBasedID_medium   = null;
    lepton2_electron_cutBasedID_tight    = null;
    lepton2_electron_mvaNoIsoID_loose    = null;
    lepton2_electron_mvaNoIsoID_wp80     = null;
    lepton2_electron_mvaNoIsoID_wp90     = null;

    lepton2_muon_CutBasedIdLoose = null;
    lepton2_muon_CutBasedIdMedium = null;
    lepton2_muon_CutBasedIdTight = null;
    lepton2_muon_CutBasedIdGlobalHighPt = null;
    lepton2_muon_CutBasedIdTrkHighPt = null;
    lepton2_muon_PFIsoLoose = null;
    lepton2_muon_PFIsoMedium = null;
    lepton2_muon_PFIsoTight = null;
    lepton2_muon_TkIsoLoose = null;
    lepton2_muon_TkIsoTight = null;
    lepton2_muon_MvaLoose = null;
    lepton2_muon_MvaMedium = null;
    lepton2_muon_MvaTight = null;

    lepton2_tau_dm  = null;
    lepton2_tau_m   = null;
    lepton2_tau_absIso           = null;
    lepton2_tau_looseCombinedIso = null;
    lepton2_tau_tightCombinedIso = null;
    lepton2_tau_looseMvaIso      = null;
    lepton2_tau_mediumMvaIso     = null;
    lepton2_tau_tightMvaIso      = null;
    lepton2_tau_VtightMvaIso     = null;
    lepton2_tau_VVtightMvaIso    = null;
    // New Raw discriminators
    lepton2_tau_againstElectronRaw            = null;
    lepton2_tau_IsoMVArun2v1DBdR03oldDMwLTraw = null;
    //
    lepton2_tau_tightMuonRejection      = null;
    lepton2_tau_looseElectronRejection  = null;
    lepton2_tau_mediumElectronRejection = null;
    lepton2_tau_tightElectronRejection  = null;
    lepton2_tau_VtightElectronRejection = null;

    edm::Handle<pat::ElectronCollection> electrons;
    event.getByToken(ElectronCollectionToken_, electrons);
    //if (!electrons.isValid()) return false;
    if (monitoringLeptons) {
        if (electrons.isValid()) {
            std::cout << "Number of electrons in collection = " << electrons->size() << std::endl;
        } else {
            std::cout << "Electrons are not valid" << std::endl;
        }
    }

    edm::Handle<pat::MuonCollection> muons;
    event.getByToken(MuonCollectionToken_, muons);
    if (monitoringLeptons) {
        if (muons.isValid()) {
            std::cout << "Number of muons in collection = " << muons->size() << std::endl;
        } else {
            std::cout << "Muons are not valid" << std::endl;
        }
    }
    //if (!muons.isValid()) return false;

    edm::Handle<pat::TauCollection> taus;
    event.getByToken(TauCollectionToken_, taus);
    if (monitoringLeptons) {
        if (taus.isValid()) {
            std::cout << "Number of taus in collection = " << taus->size() << std::endl;
        } else {
            std::cout << "Taus are not valid" << std::endl;
        }
    }

    std::vector<LeptonCandidate> LepCandidates;
    int  nLeptons = 0;
    //bool LeptonCandidateFound = false;

    #define cut(condition) if (!(condition)) continue;

    // Select electrons and muons which satisfy following criteria
    // Not in bJets, pt above threshold, eta below threshold
    std::vector<pat::Electron> ElectronCandidates;
    if (electrons.isValid() && !(*electrons).empty()) {
        //
        if (monitoringLeptons) std::cout << "Isolated electrons:" << std::endl;
        for (auto& electron: *electrons) {
            nLeptons++;
            if (monitoringLeptons) {
                std::cout << "Electron " << nLeptons << std::endl;
                std::cout << "Pt = " << electron.pt() << std::endl;
                std::cout << "Eta = " << electron.eta() << std::endl;
            }
            cut(electron.pt() > MuElePtMin);
            cut(abs(electron.eta()) < EtaMax);
            if (monitoringLeptons) {
                std::cout << "relIso = " << (electron.puppiChargedHadronIso() + electron.puppiNeutralHadronIso() + electron.puppiPhotonIso()) / electron.pt() << std::endl;
            }
            reco::TrackRef trackRef;
            trackRef = electron.track();
            if (trackRef.isNonnull()) {
                if (monitoringLeptons) {
                    std::cout << "dXY = " << trackRef->dxy(pv_position) <<
                    ", dZ = " << trackRef->dz(pv_position) << std::endl;
                }
                if (abs(electron.eta()) < 1.3) {
                    cut(std::abs(trackRef->dxy(pv_position)) <= 0.05);
                    cut(std::abs(trackRef->dz(pv_position)) <= 0.1);
                } else {
                    cut(std::abs(trackRef->dxy(pv_position)) <= 0.1);
                    cut(std::abs(trackRef->dz(pv_position)) <= 0.2);
                }
            }
            if (monitoringLeptons) {
                std::cout << "Cut based loose ID = " << electron.electronID("cutBasedElectronID-Fall17-94X-V1-loose") << std::endl;
            }
            cut(electron.electronID("cutBasedElectronID-Fall17-94X-V1-loose")); // loose for a start (MVA id is also availible)
            //
            if (monitoringLeptons) {
                std::cout << "Passed cuts" << std::endl;
            }
            double deltaRTau  = sqrt(deltaR2(tau_eta, tau_phi, electron.eta(), electron.phi()));
            if (deltaRTau <= 0.3) {
                std::cout << "In Tau cone: dR = " << deltaRTau << std::endl;
                return false;
            }
            if (deltaRTau > 0.3) {
                bool SameLepton = false;
                if(!LepCandidates.empty()) {
                    for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
                        double LepdR = sqrt(deltaR2(electron.eta(), electron.phi(), LepCandidates[ilep].Eta, LepCandidates[ilep].Phi));
                        if (LepdR < 0.4) {
                            if (monitoringLeptons) {
                                std::cout << "Same lepton with " << LepCandidates[ilep].Flavor << std::endl;
                                std::cout << "LepdR2 = " << LepdR << std::endl;
                                std::cout << "LepdR1 = " << TMath::Power(TMath::Power(electron.eta() - LepCandidates[ilep].Eta, 2) + TMath::Power(dphi(electron.phi(), LepCandidates[ilep].Phi), 2), 0.5) << std::endl;
                            }
                            SameLepton = true;
                        }
                    }
                }
                cut(!SameLepton);
                nLeptonCandidates++;
                LeptonCandidate Lepton(electron, pv_position);
                LepCandidates.push_back(Lepton);
                ElectronCandidates.push_back(electron);
                //LeptonCandidateFound = true;
                if (monitoringLeptons) std::cout << "Added to vector" << std::endl;
            }
        }
    }

    // Muons, minor cuts commented
    // declare the vector of muons
    std::vector<pat::Muon> MuonCandidates;
    if (muons.isValid() && !(*muons).empty()) {
        if (monitoringLeptons) std::cout << "Isolated muons:" << std::endl;
        for (auto& muon: *muons) {
            nLeptons++;
            if (monitoringLeptons) {
                std::cout << "Muon " << nLeptons << std::endl;
                std::cout << "Pt = " << muon.pt() << std::endl;
                std::cout << "Eta = " << muon.eta() << std::endl;
            }
            cut(muon.pt() > MuElePtMin);
            cut(abs(muon.eta()) < EtaMax);
            if (monitoringLeptons) {
                std::cout << "relIso = " << (muon.puppiChargedHadronIso() + muon.puppiNeutralHadronIso() + muon.puppiPhotonIso()) / muon.pt() << std::endl;
            }
            // Discard events with muon with pt > 15 GeV
            cut(muon.passed(reco::Muon::CutBasedIdLoose));
            cut(muon.isPFMuon());
            cut(muon.isGlobalMuon());
            cut(muon.numberOfMatchedStations() > 1);
            // check muon dXY and dZ
            reco::TrackRef trackRef;
            trackRef = muon.innerTrack();
            if (trackRef.isNonnull()) {
                if (monitoringLeptons) {
                    std::cout << "dXY = " << trackRef->dxy(pv_position) <<
                    ", dZ = " << trackRef->dz(pv_position) << std::endl;
                }
                cut(std::abs(trackRef->dxy(pv_position)) <= 0.05);
                cut(std::abs(trackRef->dz(pv_position)) <= 0.1);
                //cut(trackRef->found() > 5);  // more than 5 hits in inner tracker
            }
            //cut(muon.passed(reco::Muon::PFIsoLoose)); // I_rel < 0.25 (particle flow)
            if (monitoringLeptons) {
                std::cout << "Passed cuts" << std::endl;
            }
            double deltaRTau  = sqrt(deltaR2(tau_eta, tau_phi, muon.eta(), muon.phi()));
            if (deltaRTau <= 0.3) {
                if (monitoringLeptons) std::cout << "In Tau cone: dR = " << deltaRTau << std::endl;
                return false;
            }
            if (deltaRTau > 0.3) {
                bool SameLepton = false;
                if(!LepCandidates.empty()) {
                    for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
                        double LepdR = sqrt(deltaR2(muon.eta(), muon.phi(), LepCandidates[ilep].Eta, LepCandidates[ilep].Phi));
                        if (LepdR < 0.4) {
                            if (monitoringLeptons) {
                                std::cout << "Same lepton with " << LepCandidates[ilep].Flavor << std::endl;
                                std::cout << "LepdR2 = " << LepdR << std::endl;
                                std::cout << "LepdR1 = " << TMath::Power(TMath::Power(muon.eta() - LepCandidates[ilep].Eta, 2) + TMath::Power(dphi(muon.phi(), LepCandidates[ilep].Phi), 2), 0.5) << std::endl;
                            }
                            SameLepton = true;
                        }
                    }
                }
                cut(!SameLepton);
                LeptonCandidate Lepton(muon, pv_position);
                nLeptonCandidates++;
                LepCandidates.push_back(Lepton);
                MuonCandidates.push_back(muon);
                //LeptonCandidateFound = true;
                if (monitoringLeptons) std::cout << "Added to vector" << std::endl;
                //delete Lepton;
            }
        }
    }

    // tau as second lepton candidate
    std::vector<pat::Tau> TauCandidates;
    if (taus.isValid() && !(*taus).empty() && UseTau) {
        //
        // declare the vector of taus
        if (monitoringLeptons) {
            std::cout << "Isolated taus:" << std::endl;
        }
        for (auto& tau: *taus) {
            nLeptons++;
            if (monitoringLeptons) {
                std::cout << "Tau " << nLeptons << std::endl;
            }
            cut(tau.pt() > MuElePtMin);
            //cut(tau.tauID("byVVLooseIsolationMVArun2v1DBoldDMwLT"));
            cut(tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"));
            cut(tau.tauID("againstElectronVLooseMVA6")); // at least loose Iso MVA
            cut(tau.tauID("againstMuonLoose3")); // at least loose Iso MVA
            cut(abs(tau.eta()) < EtaMax);
            cut((pv_position - tau.vertex()).R() < tauDzMax);
            if (monitoringLeptons) {
                std::cout << "Passed cuts" << std::endl;
            }
            double deltaRTau  = sqrt(deltaR2(tau_eta, tau_phi, tau.eta(), tau.phi()));
            cut(deltaRTau >= 0.3);
            bool SameLepton = false;
            if(!LepCandidates.empty()) {
                for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
                    double LepdR = sqrt(deltaR2(tau.eta(), tau.phi(), LepCandidates[ilep].Eta, LepCandidates[ilep].Phi));
                    if (LepdR < 0.4) {
                        if (monitoringLeptons) {
                            std::cout << "Same lepton with " << LepCandidates[ilep].Flavor << std::endl;
                            std::cout << "LepdR2 = " << LepdR << std::endl;
                            std::cout << "LepdR1 = " << TMath::Power(TMath::Power(tau.eta() - LepCandidates[ilep].Eta, 2) + TMath::Power(dphi(tau.phi(), LepCandidates[ilep].Phi), 2), 0.5) << std::endl;
                        }
                        SameLepton = true;
                        if (monitoringLeptons) {
                            reco::CandidatePtrVector VectorSignalCands = tau.signalCands();
                            if (VectorSignalCands.size() > 0) {
                                std::cout << "tau dm = " << tau.decayMode() << std::endl;
                                for (unsigned l = 0; l < VectorSignalCands.size(); l++) {
                                    std::cout << "pdgID(" << l << ") = " << VectorSignalCands[l]->pdgId() << std::endl;
                                }
                            }
                        }
                    }
                }
            }
            cut(!SameLepton);
            nLeptonCandidates++;
            LeptonCandidate Lepton(tau, pv_position);
            LepCandidates.push_back(Lepton);
            TauCandidates.push_back(tau);
            if (monitoringLeptons) std::cout << "Added to vector" << std::endl;
        }
    }

    #undef cut

    if (monitoringLeptons) {
        std::cout << "Number of lepton candiates = " << LepCandidates.size() << ", " << nLeptonCandidates << std::endl;
    }
    SortLeptons(LepCandidates);

    /*
    if (monitoringLeptons) std::cout << "List of selected leptons:" << std::endl;
    for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
        if (monitoringLeptons) {
            std::cout << "lepton " << ilep << std::endl;
            LepCandidates[ilep].Monitoring();
            std::cout << "delta(Tau) = " << sqrt(deltaR2(tau_eta, tau_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
        }
    }
    */

    if (LepCandidates.size() > 0) {
        lepton1_pt       = LepCandidates[0].Pt;
        lepton1_eta      = LepCandidates[0].Eta;
        lepton1_phi      = LepCandidates[0].Phi;
        lepton1_dz       = LepCandidates[0].Dz;
        lepton1_flavor   = LepCandidates[0].Flavor;
        lepton1_charge   = LepCandidates[0].Charge;
        lepton1_E        = LepCandidates[0].FourMomentum.E();
        lepton1_trackIso = LepCandidates[0].trackIso;
        lepton1_sumPuppiIso = LepCandidates[0].puppiChargedHadronIso + LepCandidates[0].puppiNeutralHadronIso + LepCandidates[0].puppiPhotonIso;
        lepton1_sumPuppiNoLeptonIso = LepCandidates[0].puppiNoLeptonsChargedHadronIso + LepCandidates[0].puppiNoLeptonsNeutralHadronIso + LepCandidates[0].puppiNoLeptonsPhotonIso;
        if (LepCandidates.size() > 1) {
            lepton2_pt       = LepCandidates[1].Pt;
            lepton2_eta      = LepCandidates[1].Eta;
            lepton2_phi      = LepCandidates[1].Phi;
            lepton2_dz       = LepCandidates[1].Dz;
            lepton2_flavor   = LepCandidates[1].Flavor;
            lepton2_charge   = LepCandidates[1].Charge;
            lepton2_E        = LepCandidates[1].FourMomentum.E();
            lepton2_trackIso = LepCandidates[1].trackIso;
            lepton2_sumPuppiIso = LepCandidates[1].puppiChargedHadronIso + LepCandidates[1].puppiNeutralHadronIso + LepCandidates[1].puppiPhotonIso;
            lepton2_sumPuppiNoLeptonIso = LepCandidates[1].puppiNoLeptonsChargedHadronIso + LepCandidates[1].puppiNoLeptonsNeutralHadronIso + LepCandidates[1].puppiNoLeptonsPhotonIso;
        }
    }

    LepCandidates.clear();

    if (ElectronCandidates.size() > 0) {
        SortElectrons(ElectronCandidates);
        if (monitoringLeptons) {
            printElectronCandidate(ElectronCandidates[0]);
            std::cout << "delta(Tau) = " << sqrt(deltaR2(tau_eta, tau_phi, ElectronCandidates[0].eta(), ElectronCandidates[0].phi())) << std::endl;
        }
        double lepton1dR = sqrt(deltaR2(ElectronCandidates[0].eta(), ElectronCandidates[0].phi(), lepton1_eta, lepton1_phi));
        double lepton2dR = sqrt(deltaR2(ElectronCandidates[0].eta(), ElectronCandidates[0].phi(), lepton2_eta, lepton2_phi));
        if ((lepton1dR < lepton2dR && lepton1dR < 0.1) || lepton2_pt < 0) {
            if (monitoringLeptons) std::cout << "lepton 1 is electron 1" << std::endl;
            lepton1_electron_cutBasedID_loose    = ElectronCandidates[0].electronID("cutBasedElectronID-Fall17-94X-V1-loose");
            lepton1_electron_cutBasedID_medium   = ElectronCandidates[0].electronID("cutBasedElectronID-Fall17-94X-V1-medium");
            lepton1_electron_cutBasedID_tight    = ElectronCandidates[0].electronID("cutBasedElectronID-Fall17-94X-V1-tight");
            lepton1_electron_mvaNoIsoID_loose    = ElectronCandidates[0].electronID("mvaEleID-Fall17-noIso-V1-wpLoose");
            lepton1_electron_mvaNoIsoID_wp80     = ElectronCandidates[0].electronID("mvaEleID-Fall17-noIso-V1-wp80");
            lepton1_electron_mvaNoIsoID_wp90     = ElectronCandidates[0].electronID("mvaEleID-Fall17-noIso-V1-wp90");
        } else if (lepton2dR < 0.1) {
            if (monitoringLeptons) std::cout << "lepton 2 is electron 1" << std::endl;
            lepton2_electron_cutBasedID_loose    = ElectronCandidates[0].electronID("cutBasedElectronID-Fall17-94X-V1-loose");
            lepton2_electron_cutBasedID_medium   = ElectronCandidates[0].electronID("cutBasedElectronID-Fall17-94X-V1-medium");
            lepton2_electron_cutBasedID_tight    = ElectronCandidates[0].electronID("cutBasedElectronID-Fall17-94X-V1-tight");
            lepton2_electron_mvaNoIsoID_loose    = ElectronCandidates[0].electronID("mvaEleID-Fall17-noIso-V1-wpLoose");
            lepton2_electron_mvaNoIsoID_wp80     = ElectronCandidates[0].electronID("mvaEleID-Fall17-noIso-V1-wp80");
            lepton2_electron_mvaNoIsoID_wp90     = ElectronCandidates[0].electronID("mvaEleID-Fall17-noIso-V1-wp90");
        }
        //
        LeptonFound = true;
        if (ElectronCandidates.size() > 1) {
            double lepton1dR1 = sqrt(deltaR2(ElectronCandidates[0].eta(), ElectronCandidates[0].phi(), lepton1_eta, lepton1_phi));
            double lepton2dR1 = sqrt(deltaR2(ElectronCandidates[0].eta(), ElectronCandidates[0].phi(), lepton2_eta, lepton2_phi));
            double lepton1dR2 = sqrt(deltaR2(ElectronCandidates[1].eta(), ElectronCandidates[1].phi(), lepton1_eta, lepton1_phi));
            double lepton2dR2 = sqrt(deltaR2(ElectronCandidates[1].eta(), ElectronCandidates[1].phi(), lepton2_eta, lepton2_phi));
            if (lepton1dR1 < lepton2dR1 && lepton1dR1 < 0.1) {
                if (monitoringLeptons) std::cout << "lepton 1 is electron 1 (among 2)" << std::endl;
                //
                if (lepton2dR2 < 0.1) {
                    if (monitoringLeptons) std::cout << "lepton 2 is electron 2 (among 2)" << std::endl;
                    lepton2_electron_cutBasedID_loose    = ElectronCandidates[1].electronID("cutBasedElectronID-Fall17-94X-V1-loose");
                    lepton2_electron_cutBasedID_medium   = ElectronCandidates[1].electronID("cutBasedElectronID-Fall17-94X-V1-medium");
                    lepton2_electron_cutBasedID_tight    = ElectronCandidates[1].electronID("cutBasedElectronID-Fall17-94X-V1-tight");
                    lepton2_electron_mvaNoIsoID_loose    = ElectronCandidates[1].electronID("mvaEleID-Fall17-noIso-V1-wpLoose");
                    lepton2_electron_mvaNoIsoID_wp80     = ElectronCandidates[1].electronID("mvaEleID-Fall17-noIso-V1-wp80");
                    lepton2_electron_mvaNoIsoID_wp90     = ElectronCandidates[1].electronID("mvaEleID-Fall17-noIso-V1-wp90");
                }
            } else if (lepton2dR1 < 0.1) {
                if (monitoringLeptons) std::cout << "lepton 2 is electron 1 (among 2)" << std::endl;
                //
                if (lepton1dR2 < 0.1) {
                    if (monitoringLeptons) std::cout << "lepton 1 is electron 2 (among 2)" << std::endl;
                    lepton1_electron_cutBasedID_loose    = ElectronCandidates[1].electronID("cutBasedElectronID-Fall17-94X-V1-loose");
                    lepton1_electron_cutBasedID_medium   = ElectronCandidates[1].electronID("cutBasedElectronID-Fall17-94X-V1-medium");
                    lepton1_electron_cutBasedID_tight    = ElectronCandidates[1].electronID("cutBasedElectronID-Fall17-94X-V1-tight");
                    lepton1_electron_mvaNoIsoID_loose    = ElectronCandidates[1].electronID("mvaEleID-Fall17-noIso-V1-wpLoose");
                    lepton1_electron_mvaNoIsoID_wp80     = ElectronCandidates[1].electronID("mvaEleID-Fall17-noIso-V1-wp80");
                    lepton1_electron_mvaNoIsoID_wp90     = ElectronCandidates[1].electronID("mvaEleID-Fall17-noIso-V1-wp90");
                }
            }
        }
    }
    if (MuonCandidates.size() > 0) {
        SortMuons(MuonCandidates);
        if (monitoringLeptons) {
            printMuonCandidate(MuonCandidates[0]);
            std::cout << "delta(Tau) = " << sqrt(deltaR2(tau_eta, tau_phi, MuonCandidates[0].eta(), MuonCandidates[0].phi())) << std::endl;
        }
        double lepton1dR = sqrt(deltaR2(MuonCandidates[0].eta(), MuonCandidates[0].phi(), lepton1_eta, lepton1_phi));
        double lepton2dR = sqrt(deltaR2(MuonCandidates[0].eta(), MuonCandidates[0].phi(), lepton2_eta, lepton2_phi));
        if ((lepton1dR < lepton2dR && lepton1dR < 0.1) || lepton2_pt < 0) {
            if (monitoringLeptons) std::cout << "lepton 1 is muon 1" << std::endl;
            lepton1_muon_CutBasedIdLoose = MuonCandidates[0].passed(reco::Muon::CutBasedIdLoose);
            lepton1_muon_CutBasedIdMedium = MuonCandidates[0].passed(reco::Muon::CutBasedIdMedium);
            lepton1_muon_CutBasedIdTight = MuonCandidates[0].passed(reco::Muon::CutBasedIdTight);
            lepton1_muon_CutBasedIdGlobalHighPt = MuonCandidates[0].passed(reco::Muon::CutBasedIdGlobalHighPt);
            lepton1_muon_CutBasedIdTrkHighPt = MuonCandidates[0].passed(reco::Muon::CutBasedIdTrkHighPt);
            lepton1_muon_PFIsoLoose = MuonCandidates[0].passed(reco::Muon::PFIsoLoose);
            lepton1_muon_PFIsoMedium = MuonCandidates[0].passed(reco::Muon::PFIsoMedium);
            lepton1_muon_PFIsoTight = MuonCandidates[0].passed(reco::Muon::PFIsoTight);
            lepton1_muon_TkIsoLoose = MuonCandidates[0].passed(reco::Muon::TkIsoLoose);
            lepton1_muon_TkIsoTight = MuonCandidates[0].passed(reco::Muon::TkIsoTight);
            lepton1_muon_MvaLoose = MuonCandidates[0].passed(reco::Muon::MvaLoose);
            lepton1_muon_MvaMedium = MuonCandidates[0].passed(reco::Muon::MvaMedium);
            lepton1_muon_MvaTight = MuonCandidates[0].passed(reco::Muon::MvaTight);
        } else if (lepton2dR < 0.1) {
            if (monitoringLeptons) std::cout << "lepton 2 is muon 1" << std::endl;
            lepton2_muon_CutBasedIdLoose = MuonCandidates[0].passed(reco::Muon::CutBasedIdLoose);
            lepton2_muon_CutBasedIdMedium = MuonCandidates[0].passed(reco::Muon::CutBasedIdMedium);
            lepton2_muon_CutBasedIdTight = MuonCandidates[0].passed(reco::Muon::CutBasedIdTight);
            lepton2_muon_CutBasedIdGlobalHighPt = MuonCandidates[0].passed(reco::Muon::CutBasedIdGlobalHighPt);
            lepton2_muon_CutBasedIdTrkHighPt = MuonCandidates[0].passed(reco::Muon::CutBasedIdTrkHighPt);
            lepton2_muon_PFIsoLoose = MuonCandidates[0].passed(reco::Muon::PFIsoLoose);
            lepton2_muon_PFIsoMedium = MuonCandidates[0].passed(reco::Muon::PFIsoMedium);
            lepton2_muon_PFIsoTight = MuonCandidates[0].passed(reco::Muon::PFIsoTight);
            lepton2_muon_TkIsoLoose = MuonCandidates[0].passed(reco::Muon::TkIsoLoose);
            lepton2_muon_TkIsoTight = MuonCandidates[0].passed(reco::Muon::TkIsoTight);
            lepton2_muon_MvaLoose = MuonCandidates[0].passed(reco::Muon::MvaLoose);
            lepton2_muon_MvaMedium = MuonCandidates[0].passed(reco::Muon::MvaMedium);
            lepton2_muon_MvaTight = MuonCandidates[0].passed(reco::Muon::MvaTight);
        }
        if (MuonCandidates.size() > 1) {
            double lepton1dR1 = sqrt(deltaR2(MuonCandidates[0].eta(), MuonCandidates[0].phi(), lepton1_eta, lepton1_phi));
            double lepton2dR1 = sqrt(deltaR2(MuonCandidates[0].eta(), MuonCandidates[0].phi(), lepton2_eta, lepton2_phi));
            double lepton1dR2 = sqrt(deltaR2(MuonCandidates[1].eta(), MuonCandidates[1].phi(), lepton1_eta, lepton1_phi));
            double lepton2dR2 = sqrt(deltaR2(MuonCandidates[1].eta(), MuonCandidates[1].phi(), lepton2_eta, lepton2_phi));
            if (lepton1dR1 < lepton2dR1 && lepton1dR1 < 0.1) {
                if (monitoringLeptons) std::cout << "lepton 1 is muon 1 (among 2)" << std::endl;
                //
                if (lepton2dR2 < 0.1) {
                    if (monitoringLeptons) std::cout << "lepton 2 is muon 2 (among 2)" << std::endl;
                    lepton2_muon_CutBasedIdLoose = MuonCandidates[1].passed(reco::Muon::CutBasedIdLoose);
                    lepton2_muon_CutBasedIdMedium = MuonCandidates[1].passed(reco::Muon::CutBasedIdMedium);
                    lepton2_muon_CutBasedIdTight = MuonCandidates[1].passed(reco::Muon::CutBasedIdTight);
                    lepton2_muon_CutBasedIdGlobalHighPt = MuonCandidates[1].passed(reco::Muon::CutBasedIdGlobalHighPt);
                    lepton2_muon_CutBasedIdTrkHighPt = MuonCandidates[1].passed(reco::Muon::CutBasedIdTrkHighPt);
                    lepton2_muon_PFIsoLoose = MuonCandidates[1].passed(reco::Muon::PFIsoLoose);
                    lepton2_muon_PFIsoMedium = MuonCandidates[1].passed(reco::Muon::PFIsoMedium);
                    lepton2_muon_PFIsoTight = MuonCandidates[1].passed(reco::Muon::PFIsoTight);
                    lepton2_muon_TkIsoLoose = MuonCandidates[1].passed(reco::Muon::TkIsoLoose);
                    lepton2_muon_TkIsoTight = MuonCandidates[1].passed(reco::Muon::TkIsoTight);
                    lepton2_muon_MvaLoose = MuonCandidates[1].passed(reco::Muon::MvaLoose);
                    lepton2_muon_MvaMedium = MuonCandidates[1].passed(reco::Muon::MvaMedium);
                    lepton2_muon_MvaTight = MuonCandidates[1].passed(reco::Muon::MvaTight); 
                }
            } else if (lepton2dR1 < 0.1) {
                if (monitoringLeptons) std::cout << "lepton 2 is muon 1 (among 2)" << std::endl;
                //
                if (lepton1dR2 < 0.1) {
                    if (monitoringLeptons) std::cout << "lepton 1 is muon 2 (among 2)" << std::endl;
                    lepton1_muon_CutBasedIdLoose = MuonCandidates[1].passed(reco::Muon::CutBasedIdLoose);
                    lepton1_muon_CutBasedIdMedium = MuonCandidates[1].passed(reco::Muon::CutBasedIdMedium);
                    lepton1_muon_CutBasedIdTight = MuonCandidates[1].passed(reco::Muon::CutBasedIdTight);
                    lepton1_muon_CutBasedIdGlobalHighPt = MuonCandidates[1].passed(reco::Muon::CutBasedIdGlobalHighPt);
                    lepton1_muon_CutBasedIdTrkHighPt = MuonCandidates[1].passed(reco::Muon::CutBasedIdTrkHighPt);
                    lepton1_muon_PFIsoLoose = MuonCandidates[1].passed(reco::Muon::PFIsoLoose);
                    lepton1_muon_PFIsoMedium = MuonCandidates[1].passed(reco::Muon::PFIsoMedium);
                    lepton1_muon_PFIsoTight = MuonCandidates[1].passed(reco::Muon::PFIsoTight);
                    lepton1_muon_TkIsoLoose = MuonCandidates[1].passed(reco::Muon::TkIsoLoose);
                    lepton1_muon_TkIsoTight = MuonCandidates[1].passed(reco::Muon::TkIsoTight);
                    lepton1_muon_MvaLoose = MuonCandidates[1].passed(reco::Muon::MvaLoose);
                    lepton1_muon_MvaMedium = MuonCandidates[1].passed(reco::Muon::MvaMedium);
                    lepton1_muon_MvaTight = MuonCandidates[1].passed(reco::Muon::MvaTight);
                }
            }
        }
        //
        LeptonFound = true;
    }
    if (TauCandidates.size() > 0) {
        SortTaus(TauCandidates);
        if (monitoringLeptons) {
            TauMonitor (TauCandidates[0], DeepTau, pv_position);
            std::cout << "delta(Tau)   = " << sqrt(deltaR2(tau_eta, tau_phi, TauCandidates[0].eta(), TauCandidates[0].phi())) << std::endl;
        }
        double lepton1dR = sqrt(deltaR2(MuonCandidates[0].eta(), MuonCandidates[0].phi(), lepton1_eta, lepton1_phi));
        double lepton2dR = sqrt(deltaR2(MuonCandidates[0].eta(), MuonCandidates[0].phi(), lepton2_eta, lepton2_phi));
        if ((lepton1dR < lepton2dR && lepton1dR < 0.1) || lepton2_pt < 0) {
            if (monitoringLeptons) std::cout << "lepton 1 is tau 1" << std::endl;
            lepton1_tau_dm                      = TauCandidates[0].decayMode();
            lepton1_tau_m                       = TauCandidates[0].mass();
            lepton1_tau_absIso                  = TauCandidates[0].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
            lepton1_tau_looseCombinedIso        = TauCandidates[0].tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
            lepton1_tau_tightCombinedIso        = TauCandidates[0].tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
            lepton1_tau_looseMvaIso             = TauCandidates[0].tauID("byLooseIsolationMVArun2v1DBoldDMwLT"); // there are more possible discriminators 
            lepton1_tau_mediumMvaIso            = TauCandidates[0].tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
            lepton1_tau_tightMvaIso             = TauCandidates[0].tauID("byTightIsolationMVArun2v1DBoldDMwLT");
            lepton1_tau_VtightMvaIso            = TauCandidates[0].tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
            lepton1_tau_VVtightMvaIso           = TauCandidates[0].tauID("byVVTightIsolationMVArun2v1DBoldDMwLT");
            // New Raw discriminators
            lepton1_tau_againstElectronRaw            = TauCandidates[0].tauID("againstElectronMVA6Raw");
            lepton1_tau_IsoMVArun2v1DBdR03oldDMwLTraw = TauCandidates[0].tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw");
            //
            lepton1_tau_tightMuonRejection      = TauCandidates[0].tauID("againstMuonTight3");
            lepton1_tau_looseElectronRejection  = TauCandidates[0].tauID("againstElectronLooseMVA6");
            lepton1_tau_mediumElectronRejection = TauCandidates[0].tauID("againstElectronMediumMVA6");
            lepton1_tau_tightElectronRejection  = TauCandidates[0].tauID("againstElectronTightMVA6");
            lepton1_tau_VtightElectronRejection = TauCandidates[0].tauID("againstElectronVTightMVA6");
        } else if (lepton2dR < 0.1) {
            if (monitoringLeptons) std::cout << "lepton 2 is tau 1" << std::endl;
            lepton2_tau_dm                      = TauCandidates[0].decayMode();
            lepton2_tau_m                       = TauCandidates[0].mass();
            lepton2_tau_absIso                  = TauCandidates[0].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
            lepton2_tau_looseCombinedIso        = TauCandidates[0].tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
            lepton2_tau_tightCombinedIso        = TauCandidates[0].tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
            lepton2_tau_looseMvaIso             = TauCandidates[0].tauID("byLooseIsolationMVArun2v1DBoldDMwLT"); // there are more possible discriminators 
            lepton2_tau_mediumMvaIso            = TauCandidates[0].tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
            lepton2_tau_tightMvaIso             = TauCandidates[0].tauID("byTightIsolationMVArun2v1DBoldDMwLT");
            lepton2_tau_VtightMvaIso            = TauCandidates[0].tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
            lepton2_tau_VVtightMvaIso           = TauCandidates[0].tauID("byVVTightIsolationMVArun2v1DBoldDMwLT");
            // New Raw discriminators
            lepton2_tau_againstElectronRaw            = TauCandidates[0].tauID("againstElectronMVA6Raw");
            lepton2_tau_IsoMVArun2v1DBdR03oldDMwLTraw = TauCandidates[0].tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw");
            //
            lepton2_tau_tightMuonRejection      = TauCandidates[0].tauID("againstMuonTight3");
            lepton2_tau_looseElectronRejection  = TauCandidates[0].tauID("againstElectronLooseMVA6");
            lepton2_tau_mediumElectronRejection = TauCandidates[0].tauID("againstElectronMediumMVA6");
            lepton2_tau_tightElectronRejection  = TauCandidates[0].tauID("againstElectronTightMVA6");
            lepton2_tau_VtightElectronRejection = TauCandidates[0].tauID("againstElectronVTightMVA6");
        }
        LeptonFound = true;
        if (TauCandidates.size() > 1) {
            double lepton1dR1 = sqrt(deltaR2(TauCandidates[0].eta(), TauCandidates[0].phi(), lepton1_eta, lepton1_phi));
            double lepton2dR1 = sqrt(deltaR2(TauCandidates[0].eta(), TauCandidates[0].phi(), lepton2_eta, lepton2_phi));
            double lepton1dR2 = sqrt(deltaR2(TauCandidates[1].eta(), TauCandidates[1].phi(), lepton1_eta, lepton1_phi));
            double lepton2dR2 = sqrt(deltaR2(TauCandidates[1].eta(), TauCandidates[1].phi(), lepton2_eta, lepton2_phi));
            if (lepton1dR1 < lepton2dR1 && lepton1dR1 < 0.1) {
                if (monitoringLeptons) std::cout << "lepton 1 is tau 1 (among 2)" << std::endl;
                //
                if (lepton2dR2 < 0.1) {
                    if (monitoringLeptons) std::cout << "lepton 2 is tau 2 (among 2)" << std::endl;
                    lepton2_tau_dm                      = TauCandidates[1].decayMode();
                    lepton2_tau_m                       = TauCandidates[1].mass();
                    lepton2_tau_absIso                  = TauCandidates[1].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
                    lepton2_tau_looseCombinedIso        = TauCandidates[1].tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
                    lepton2_tau_tightCombinedIso        = TauCandidates[1].tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
                    lepton2_tau_looseMvaIso             = TauCandidates[1].tauID("byLooseIsolationMVArun2v1DBoldDMwLT"); // there are more possible discriminators 
                    lepton2_tau_mediumMvaIso            = TauCandidates[1].tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
                    lepton2_tau_tightMvaIso             = TauCandidates[1].tauID("byTightIsolationMVArun2v1DBoldDMwLT");
                    lepton2_tau_VtightMvaIso            = TauCandidates[1].tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
                    lepton2_tau_VVtightMvaIso           = TauCandidates[1].tauID("byVVTightIsolationMVArun2v1DBoldDMwLT");
                    // New Raw discriminators
                    lepton2_tau_againstElectronRaw            = TauCandidates[1].tauID("againstElectronMVA6Raw");
                    lepton2_tau_IsoMVArun2v1DBdR03oldDMwLTraw = TauCandidates[1].tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw");
                    //
                    lepton2_tau_tightMuonRejection      = TauCandidates[1].tauID("againstMuonTight3");
                    lepton2_tau_looseElectronRejection  = TauCandidates[1].tauID("againstElectronLooseMVA6");
                    lepton2_tau_mediumElectronRejection = TauCandidates[1].tauID("againstElectronMediumMVA6");
                    lepton2_tau_tightElectronRejection  = TauCandidates[1].tauID("againstElectronTightMVA6");
                    lepton2_tau_VtightElectronRejection = TauCandidates[1].tauID("againstElectronVTightMVA6");
                }
            } else if (lepton2dR1 < 0.1) {
                if (monitoringLeptons) std::cout << "lepton 2 is tau 1 (among 2)" << std::endl;
                //
                if (lepton1dR2 < 0.1) {
                    if (monitoringLeptons) std::cout << "lepton 1 is tau 2 (among 2)" << std::endl;
                    lepton1_tau_dm                      = TauCandidates[1].decayMode();
                    lepton1_tau_m                       = TauCandidates[1].mass();
                    lepton1_tau_absIso                  = TauCandidates[1].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
                    lepton1_tau_looseCombinedIso        = TauCandidates[1].tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
                    lepton1_tau_tightCombinedIso        = TauCandidates[1].tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
                    lepton1_tau_looseMvaIso             = TauCandidates[1].tauID("byLooseIsolationMVArun2v1DBoldDMwLT"); // there are more possible discriminators 
                    lepton1_tau_mediumMvaIso            = TauCandidates[1].tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
                    lepton1_tau_tightMvaIso             = TauCandidates[1].tauID("byTightIsolationMVArun2v1DBoldDMwLT");
                    lepton1_tau_VtightMvaIso            = TauCandidates[1].tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
                    lepton1_tau_VVtightMvaIso           = TauCandidates[1].tauID("byVVTightIsolationMVArun2v1DBoldDMwLT");
                    // New Raw discriminators
                    lepton1_tau_againstElectronRaw            = TauCandidates[1].tauID("againstElectronMVA6Raw");
                    lepton1_tau_IsoMVArun2v1DBdR03oldDMwLTraw = TauCandidates[1].tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw");
                    //
                    lepton1_tau_tightMuonRejection      = TauCandidates[1].tauID("againstMuonTight3");
                    lepton1_tau_looseElectronRejection  = TauCandidates[1].tauID("againstElectronLooseMVA6");
                    lepton1_tau_mediumElectronRejection = TauCandidates[1].tauID("againstElectronMediumMVA6");
                    lepton1_tau_tightElectronRejection  = TauCandidates[1].tauID("againstElectronTightMVA6");
                    lepton1_tau_VtightElectronRejection = TauCandidates[1].tauID("againstElectronVTightMVA6");
                }
            }
        }
    }
//    }

    if (LeptonRequired) {
        return (LeptonFound);
    } else {
        return true;        
    } 
};

void TTbarTauLepton::CountTracks(const edm::Event& event) {
    nTrks = 0;
    edm::Handle<pat::IsolatedTrackCollection> tracks;
    event.getByToken(TrackToken_, tracks);
    if (!tracks.isValid()) return;
    nTrks = tracks->size();
};

void TTbarTauLepton::AddWT(const edm::Event& iEvent) {

    edm::Handle<double> WTHandle;
    iEvent.getByToken(TauSpinnerWTToken_, WTHandle);
    edm::Handle<double> WTFlipHandle;
    iEvent.getByToken(TauSpinnerWTFlipToken_, WTFlipHandle);
    edm::Handle<double> WThminusHandle;
    iEvent.getByToken(TauSpinnerWThminusToken_, WThminusHandle);
    edm::Handle<double> WThplusHandle;
    iEvent.getByToken(TauSpinnerWThplusToken_, WThplusHandle);
    edm::Handle<bool>   WTisValidHandle;
    iEvent.getByToken(TauSpinnerWTisValidToken_, WTisValidHandle);
    edm::Handle<int>   TauMotherHandle;
    iEvent.getByToken(TauSpinnerMotherToken_, TauMotherHandle);

    if (!WTisValidHandle.isValid()) {
        if (monitoring) std::cout << "WTisValidHandle is not valid" << std::endl;
        return;
    }
    WTisValid = *WTisValidHandle;
    if (WTisValid) {
        WT        = *WTHandle;
        WTFlip    = *WTFlipHandle;
        WThminus  = *WThminusHandle;
        WThplus   = *WThplusHandle;
        TauSpinnerMother = *TauMotherHandle;
    } else {
        if (monitoring) std::cout << "WT Collections are not valid" << std::endl;
        WT        = null;
        WTFlip    = null;
        WThminus  = null;
        WThplus   = null;
        TauSpinnerMother = null;
    }
};

void TTbarTauLepton::AddWTData(const edm::Event& iEvent) {

    int One = 1;
    WTisValid = One;
    WT        = One;
    WTFlip    = One;
    WThminus  = One;
    WThplus   = One;
    TauSpinnerMother = null;
};

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TTbarTauLepton::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.add("TTbarTauLepton", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TTbarTauLepton);