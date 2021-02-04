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
	//void AddGenMET         (const edm::Event&);
	bool AddLepton         (const edm::Event&);
	bool AddMET            (const edm::Event&);
	bool AddVertex         (const edm::Event&);
	void CountTracks       (const edm::Event&);
	bool JetPtSum          (const edm::Event&);
	void AddWT             (const edm::Event&);
	void AddWTData         (const edm::Event&);
	void GenTauDM          (const edm::Event&);
	bool TriggerOK         (const edm::Event&);
	void AddEmptyTau       (const edm::Event&);

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
	std::vector<std::string>  trigNamesTau;
	std::vector<std::string>  trigNamesJetHT;
	std::vector<std::string>  trigNamesMET;
	std::vector<std::string>  trigNamesBTagCSV;
	std::vector<std::string>  trigNamesBTagMu;
	std::vector<std::string>  trigNamesSingleMuon;
	std::vector<std::string>  trigNamesSingleElectron;
	std::vector<std::string>  trigNamesTarget;
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
	int tau_VVlooseMvaIso;
	int tau_VlooseMvaIso;
	int tau_looseMvaIso;
	int tau_mediumMvaIso;
	int tau_tightMvaIso;
	int tau_VtightMvaIso;
	int tau_VVtightMvaIso;
	int tau_looseMuonRejection;
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
	double genb1Fromt;
	double genb2Fromt;
	int genb1Mother;
	int genb2Mother;

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
	double LeadingJet_pt;
	double LeadingJet_eta;
	double LeadingJet_phi;
	double LeadingJet_m;
	double LeadingJet_btag;
	double SubLeadingJet_pt;
	double SubLeadingJet_eta;
	double SubLeadingJet_phi;
	double SubLeadingJet_m;
	double SubLeadingJet_btag;

	// BJets
	double BJet1_bprob;
	double BJet1_pt;
	double BJet1_eta;
	double BJet1_phi;
	double BJet1_E;
	double BJet2_bprob;
	double BJet2_pt;
	double BJet2_eta;
	double BJet2_phi;
	double BJet2_E;
	double TauJet_pt;
	double TauJet_eta;
	double TauJet_phi;
	double TauJet_E;

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
	double lepton1_tauAbsIso;
	int    nLeptonCandidates;

	// Electron candidate
	double electron_pt;
	double electron_eta;
	double electron_phi;
	double electron_E;
	int    electron_charge;
	float  electron_trackIso;
	float  electron_sumPuppiIso;
	float  electron_sumPuppiNoLeptonIso;
	float  electron_caloIso;
	int    electron_cutBasedID_loose;
	int    electron_cutBasedID_medium;
	int    electron_cutBasedID_tight;
	//int    electron_mvaIsoID_loose;
	//int    electron_mvaIsoID_wp80;
	//int    electron_mvaIsoID_wp90;
	int    electron_mvaNoIsoID_loose;
	int    electron_mvaNoIsoID_wp80;
	int    electron_mvaNoIsoID_wp90;
	//int    electron_cutBasedID;
	//int    electron_mvaIsoID;
	//int    electron_mvaNoIsoID;

	// Muon candidate
	double muon_pt;
	double muon_eta;
	double muon_phi;
	double muon_E;
	int    muon_charge;
	float  muon_trackIso;
	float  muon_sumPuppiIso;
	float  muon_sumPuppiNoLeptonIso;
	float  muon_caloIso;
	int    muon_CutBasedIdLoose;
	int    muon_CutBasedIdMedium;
	int    muon_CutBasedIdTight;
	int    muon_CutBasedIdGlobalHighPt;
	int    muon_CutBasedIdTrkHighPt;
	int    muon_PFIsoLoose;
	int    muon_PFIsoMedium;
	int    muon_PFIsoTight;
	int    muon_TkIsoLoose;
	int    muon_TkIsoTight;
	int    muon_MvaLoose;
	int    muon_MvaMedium;
	int    muon_MvaTight;

	// tau 2
	double tau2_pt;
	double tau2_eta;
	double tau2_phi;
	double tau2_dm;
	double tau2_q;
	double tau2_m;
	double tau2_dz;
	
	// MVA raw
	double tau2_absIso;
	double tau2_againstElectronRaw;
	double tau2_IsoMVArun2v1DBdR03oldDMwLTraw;
	double tau2_IsoMVArun2v1DBnewDMwLTraw;
	double tau2_IsoMVArun2v1DBoldDMwLTraw;
	double tau2_IsoMVArun2v1PWdR03oldDMwLTraw;
	double tau2_IsoMVArun2v1PWnewDMwLTraw;
	double tau2_IsoMVArun2v1PWoldDMwLTraw;

	// Deep raw
	double tau2_Deep2017v2ElectronRejection;
	double tau2_Deep2017v2MuonRejection;
	double tau2_Deep2017v2JetRejection;

	// Deep WP
	int tau2_VVLooseDeepTau2017v2VSjet;
	int tau2_VLooseDeepTau2017v2VSjet;
	int tau2_LooseDeepTau2017v2VSjet;
	int tau2_MediumDeepTau2017v2VSjet;
	int tau2_TightDeepTau2017v2VSjet;
	int tau2_VTightDeepTau2017v2VSjet;
	int tau2_VVTightDeepTau2017v2VSjet;

	int tau2_LooseDeepTau2017v2VSmu;
	int tau2_MediumDeepTau2017v2VSmu;
	int tau2_TightDeepTau2017v2VSmu;

	int tau2_VVLooseDeepTau2017v2VSe;
	int tau2_VLooseDeepTau2017v2VSe;
	int tau2_LooseDeepTau2017v2VSe;
	int tau2_MediumDeepTau2017v2VSe;
	int tau2_TightDeepTau2017v2VSe;
	int tau2_VTightDeepTau2017v2VSe;
	int tau2_VVTightDeepTau2017v2VSe;

	// MVA WP
	int tau2_looseCombinedIso;
	int tau2_mediumCombinedIso;
	int tau2_tightCombinedIso;
	int tau2_VVlooseMvaIso;
	int tau2_VlooseMvaIso;
	int tau2_looseMvaIso;
	int tau2_mediumMvaIso;
	int tau2_tightMvaIso;
	int tau2_VtightMvaIso;
	int tau2_VVtightMvaIso;
	int tau2_looseMuonRejection;
	int tau2_tightMuonRejection;
	int tau2_looseElectronRejection;
	int tau2_mediumElectronRejection;
	int tau2_tightElectronRejection;
	int tau2_VtightElectronRejection;

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

	//Gen MET
	double TrueMetPt;
	double TrueMetEta;
	double TrueMetPhi;
	double TrueMetEnergy;

	// Trigger matching
	int nTauTriggers;
	int nJetHTTriggers;
	int nMETTriggers;
	int nBTagCSVTriggers;
	int nBTagMuTriggers;
	int nSingleMuonTriggers;
	int nSingleElectronTriggers;
	int nTargetTriggers;
	
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
	bool fullMC;
	bool useHLT;
	bool useTargetHLT;
	bool TauSpinnerOn;
	double tauPtMin;
	double piPtMin;
	double tauEtaMax;
	double tauDzMax;
	double METcut;
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
	int NrequiredBJets;
	int NrequiredJets;
	int NrequiredLeptons;
	int requiredLeptonPDGID;

	int iT;
	int nTrigger = 0;
	int nTau1 = 0;
	int nTauJet = 0;
	int nJet = 0;
	int nLepton = 0;
	int nEvent = 0;
	int nPassed = 0;
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

	isMC						= iConfig.getParameter<bool>("isMC");
	fullMC						= iConfig.getParameter<bool>("fullMC");
	useHLT                      = iConfig.getParameter<bool>("useHLT");
	useTargetHLT                = iConfig.getParameter<bool>("useTargetHLT");
	monitoring					= iConfig.getParameter<bool>("monitoring");
	monitoringHLT				= iConfig.getParameter<bool>("monitoringHLT");
	monitoringTau				= iConfig.getParameter<bool>("monitoringTau");
	monitoringGen				= iConfig.getParameter<bool>("monitoringGen");
	monitoringJets				= iConfig.getParameter<bool>("monitoringJets");
	monitoringBJets				= iConfig.getParameter<bool>("monitoringBJets");
	monitoringLeptons			= iConfig.getParameter<bool>("monitoringLeptons");
	monitoringMET               = iConfig.getParameter<bool>("monitoringMET");
	TauSpinnerOn                = iConfig.getParameter<bool>("TauSpinnerOn");
	tauPtMin 					= iConfig.getParameter<double>("tauPtMin");
	piPtMin 					= iConfig.getParameter<double>("piPtMin");
	tauEtaMax 					= iConfig.getParameter<double>("tauEtaMax");
	tauDzMax 					= iConfig.getParameter<double>("tauDzMax");
	METcut                      = iConfig.getParameter<double>("METcut");
	looseTauID					= iConfig.getParameter<bool>("looseTauID");
	BJetPtMin                   = iConfig.getParameter<double>("BJetPtMin"); // 30
	JetEtaMax                   = iConfig.getParameter<double>("JetEtaMax");
	MuElePtMin                  = iConfig.getParameter<double>("MuElePtMin"); // 20
	EtaMax                      = iConfig.getParameter<double>("EtaMax"); // 2.4
	DeepTau  					= iConfig.getParameter<bool>("DeepTau");
	null 						= iConfig.getParameter<double>("null");
	NrequiredBJets              = iConfig.getParameter<int>("NrequiredBJets"); // Loose b-tagged jets (1 or 2 for ttbar)
	NrequiredJets               = iConfig.getParameter<int>("NrequiredJets"); // at least 2 i think
	NrequiredLeptons            = iConfig.getParameter<int>("NrequiredLeptons");
	requiredLeptonPDGID         = iConfig.getParameter<int>("requiredLeptonPDGID"); // 11, 13 or 15

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
	theTriggerResultsLabel		      = edm::InputTag("TriggerResults","","HLT");
	trigNamesTau                      = iConfig.getParameter<std::vector<std::string>>("TauTriggers");
	trigNamesJetHT                    = iConfig.getParameter<std::vector<std::string>>("JetHTTriggers");
	trigNamesMET                      = iConfig.getParameter<std::vector<std::string>>("METTriggers");
	trigNamesBTagCSV                  = iConfig.getParameter<std::vector<std::string>>("BTagCSVTriggers");
	trigNamesBTagMu                   = iConfig.getParameter<std::vector<std::string>>("BTagMuTriggers");
	trigNamesSingleMuon               = iConfig.getParameter<std::vector<std::string>>("SingleMuonTriggers");
	trigNamesSingleElectron           = iConfig.getParameter<std::vector<std::string>>("SingleElectronTriggers");
	trigNamesTarget                   = iConfig.getParameter<std::vector<std::string>>("TriggersTarget");
	std::string PackedCandidateCollection = iConfig.getParameter<std::string>("PackedCandidateCollection");
	
	TauCollectionToken_ 		= consumes<pat::TauCollection>(edm::InputTag(tauCollection));
	MuonCollectionToken_ 		= consumes<pat::MuonCollection>(edm::InputTag(muonCollection));
	ElectronCollectionToken_	= consumes<pat::ElectronCollection>(edm::InputTag(electronCollection));
	PuppiJetCollectionToken_    = consumes<pat::JetCollection>(edm::InputTag(PuppijetCollection));
	MetCollectionToken_ 		= consumes<pat::METCollection>(edm::InputTag(metCollection));
	PuppiMetCollectionToken_    = consumes<pat::METCollection>(edm::InputTag(PuppimetCollection));
	PVToken_ 			        = consumes<reco::VertexCollection>(edm::InputTag(vertexCollection));
	SVToken_                    = consumes<reco::VertexCompositePtrCandidateCollection>(edm::InputTag(SVCollection));
	GenParticleToken_           = consumes<reco::GenParticleCollection>(edm::InputTag(genParticleCollection));
	TrackToken_			        = consumes<pat::IsolatedTrackCollection>(edm::InputTag(trackCollection));
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

	// Operation with triggers (filtering, scaling)
	if (!TriggerOK(event)) {
		nTrigger++;
		if (monitoring) std::cout << "Trigger" << std::endl;
		if (!fullMC) return;
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
		if (!fullMC) {
			if (monitoring) std::cout << "Tau" << std::endl;
			return;
		} else if (isMC) AddEmptyTau(event);
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

	std::cout << "Events  = " << nEvent << std::endl;
	std::cout << "Trigger = " << nTrigger << std::endl;
	std::cout << "Tau     = " << nTau1 << std::endl;
	std::cout << "Lepton  = " << nLepton << std::endl;
	std::cout << "Jet     = " << nJet << std::endl;
	std::cout << "Passed  = " << nPassed << std::endl;

	TLorentzVector pi0, pi1;
	pi0.SetPtEtaPhiM(piZero_pt, piZero_eta, piZero_phi, piZero_m);
	pi1.SetPtEtaPhiM(piChar_pt, piChar_eta, piChar_phi, piChar_m);

	pipiMass = (pi0 + pi1).M();
	//upsilon  = (pi1.E() - pi0.E()) / tau_pt;
	upsilon  = 2 * piChar_pt / tau_pt - 1;
	dPhi	        = dphi(tau_phi, met_phi);
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
	tree->Branch("tau_VVlooseMvaIso",&tau_VVlooseMvaIso,"tau_VVlooseMvaIso/I");
	tree->Branch("tau_VlooseMvaIso",&tau_VlooseMvaIso,"tau_VlooseMvaIso/I");
	tree->Branch("tau_looseMvaIso",&tau_looseMvaIso,"tau_looseMvaIso/I");
	tree->Branch("tau_mediumMvaIso",&tau_mediumMvaIso,"tau_mediumMvaIso/I");
	tree->Branch("tau_tightMvaIso",&tau_tightMvaIso,"tau_tightMvaIso/I");
	tree->Branch("tau_VtightMvaIso",&tau_VtightMvaIso,"tau_VtightMvaIso/I");
	tree->Branch("tau_VVtightMvaIso",&tau_VVtightMvaIso,"tau_VVtightMvaIso/I");
	tree->Branch("tau_looseMuonRejection",&tau_looseMuonRejection,"tau_looseMuonRejection/I");
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
	tree->Branch("genLeptonFromWFromt", &genLeptonFromWFromt, "genLeptonFromWFromt/D");
	tree->Branch("genb1Fromt", &genb1Fromt, "genb1Fromt/D");
	tree->Branch("genb2Fromt", &genb2Fromt, "genb2Fromt/D");
	tree->Branch("genb1Mother", &genb1Mother, "genb1Mother/I");
	tree->Branch("genb2Mother", &genb2Mother, "genb2Mother/I");
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
	// PuppiJets
	tree->Branch("PuppijetPtSum30", &PuppijetPtSum30, "PuppijetPtSum30/D");
	tree->Branch("PuppijetPtSum20", &PuppijetPtSum20, "PuppijetPtSum20/D");
	tree->Branch("PuppijetPtSum30PV", &PuppijetPtSum30PV, "PuppijetPtSum30PV/D");
	tree->Branch("PuppijetPtSum20PV", &PuppijetPtSum20PV, "PuppijetPtSum20PV/D");
	tree->Branch("nPuppiJets20", &nPuppiJets20, "nPuppiJets20/I");
	tree->Branch("nPuppiJets20PV", &nPuppiJets20PV, "nPuppiJets20PV/I");
	tree->Branch("nPuppiJets30", &nPuppiJets30, "nPuppiJets30/I");
	tree->Branch("nPuppiJets30PV", &nPuppiJets30PV, "nPuppiJets30PV/I");
	tree->Branch("nLooseBtagedPuppiJets", &nLooseBtagedPuppiJets, "nLooseBtagedPuppiJets/I");
	tree->Branch("nMediumBtagedPuppiJets", &nMediumBtagedPuppiJets, "nMediumBtagedPuppiJets/I");
	tree->Branch("nTightBtagedPuppiJets", &nTightBtagedPuppiJets, "nTightBtagedPuppiJets/I");
	tree->Branch("nLooseBtagedPuppiJetsPV", &nLooseBtagedPuppiJetsPV, "nLooseBtagedPuppiJetsPV/I");
	tree->Branch("nMediumBtagedPuppiJetsPV", &nMediumBtagedPuppiJetsPV, "nMediumBtagedPuppiJetsPV/I");
	tree->Branch("nTightBtagedPuppiJetsPV", &nTightBtagedPuppiJetsPV, "nTightBtagedPuppiJetsPV/I");

	// Leading, subleading, bJet parameters
	tree->Branch("LeadingJet_pt", &LeadingJet_pt, "LeadingJet_pt/D");
	tree->Branch("LeadingJet_eta", &LeadingJet_eta, "LeadingJet_eta/D");
	tree->Branch("LeadingJet_phi", &LeadingJet_phi, "LeadingJet_phi/D");
	tree->Branch("LeadingJet_m", &LeadingJet_m, "LeadingJet_m/D");
	tree->Branch("LeadingJet_btag", &LeadingJet_btag, "LeadingJet_btag/D");
	tree->Branch("SubLeadingJet_pt", &SubLeadingJet_pt, "SubLeadingJet_pt/D");
	tree->Branch("SubSubLeadingJet_eta", &SubLeadingJet_eta, "SubLeadingJet_eta/D");
	tree->Branch("SubLeadingJet_phi", &SubLeadingJet_phi, "SubLeadingJet_phi/D");
	tree->Branch("SubLeadingJet_m", &SubLeadingJet_m, "SubLeadingJet_m/D");
	tree->Branch("SubLeadingJet_btag", &SubLeadingJet_btag, "SubLeadingJet_btag/D");

	//BJets and leptons
	if (NrequiredLeptons >= 0) {
		tree->Branch("lepton1_pt", &lepton1_pt, "lepton1_pt/D");
		tree->Branch("lepton1_E", &lepton1_E, "lepton1_E/D");
		tree->Branch("lepton1_eta", &lepton1_eta, "lepton1_eta/D");
		tree->Branch("lepton1_phi", &lepton1_phi, "lepton1_phi/D");
		tree->Branch("lepton1_dz", &lepton1_dz, "lepton1_dz/D");
		tree->Branch("lepton1_flavor", &lepton1_flavor, "lepton1_flavor/I");
		tree->Branch("lepton1_charge", &lepton1_charge, "lepton1_charge/I");
		tree->Branch("lepton1_trackIso", &lepton1_trackIso, "lepton1_trackIso/F");
		tree->Branch("lepton1_sumPuppiIso", &lepton1_sumPuppiIso, "lepton1_sumPuppiIso/F");
		tree->Branch("lepton1_tauAbsIso", &lepton1_tauAbsIso, "lepton1_tauAbsIso/D");
		tree->Branch("lepton1_sumPuppiNoLeptonIso", &lepton1_sumPuppiNoLeptonIso, "lepton1_sumPuppiNoLeptonIso/F");
	}
	if (NrequiredBJets >= 0) {
		tree->Branch("BJet1_pt", &BJet1_pt, "BJet1_pt/D");
		tree->Branch("BJet1_E", &BJet1_E, "BJet1_E/D");
		tree->Branch("BJet1_eta", &BJet1_eta, "BJet1_eta/D");
		tree->Branch("BJet1_phi", &BJet1_phi, "BJet1_phi/D");
		tree->Branch("BJet1_bprob", &BJet1_bprob, "BJet1_bprob/D");
		tree->Branch("BJet2_pt", &BJet2_pt, "BJet2_pt/D");
		tree->Branch("BJet2_E", &BJet2_E, "BJet2_E/D");
		tree->Branch("BJet2_eta", &BJet2_eta, "BJet2_eta/D");
		tree->Branch("BJet2_phi", &BJet2_phi, "BJet2_phi/D");
		tree->Branch("BJet2_bprob", &BJet2_bprob, "BJet2_bprob/D");
	}
	if (requiredLeptonPDGID == 11) {
		tree->Branch("electron_pt", &electron_pt, "electron_pt/D");
		tree->Branch("electron_E", &electron_E, "electron_E/D");
		tree->Branch("electron_eta", &electron_eta, "electron_eta/D");
		tree->Branch("electron_phi", &electron_phi, "electron_phi/D");
		tree->Branch("electron_charge", &electron_charge, "electron_charge/I");
		tree->Branch("electron_trackIso", &electron_trackIso, "electron_trackIso/F");
		tree->Branch("electron_sumPuppiIso", &electron_sumPuppiIso, "electron_sumPuppiIso/F");
		tree->Branch("electron_sumPuppiNoLeptonIso", &electron_sumPuppiNoLeptonIso, "electron_sumPuppiNoLeptonIso/F");
		tree->Branch("electron_caloIso", &electron_caloIso, "electron_caloIso/F");
		tree->Branch("electron_cutBasedID_loose", &electron_cutBasedID_loose, "electron_cutBasedID_loose/I");
		tree->Branch("electron_cutBasedID_medium", &electron_cutBasedID_medium, "electron_cutBasedID_medium/I");
		tree->Branch("electron_cutBasedID_tight", &electron_cutBasedID_tight, "electron_cutBasedID_tight/I");
		//tree->Branch("electron_mvaIsoID_loose", &electron_mvaIsoID_loose, "electron_mvaIsoID_loose/I");
		//tree->Branch("electron_mvaIsoID_wp80", &electron_mvaIsoID_wp80, "electron_mvaIsoID_wp80/I");
		//tree->Branch("electron_mvaIsoID_wp90", &electron_mvaIsoID_wp90, "electron_mvaIsoID_wp90/I");
		tree->Branch("electron_mvaNoIsoID_loose", &electron_mvaNoIsoID_loose, "electron_mvaNoIsoID_loose/I");
		tree->Branch("electron_mvaNoIsoID_wp80", &electron_mvaNoIsoID_wp80, "electron_mvaNoIsoID_wp80/I");
		tree->Branch("electron_mvaNoIsoID_wp90", &electron_mvaNoIsoID_wp90, "electron_mvaNoIsoID_wp90/I");
		//tree->Branch("electron_cutBasedID", &electron_cutBasedID, "electron_cutBasedID/I");
		//tree->Branch("electron_mvaIsoID", &electron_mvaIsoID, "electron_mvaIsoID/I");
		//tree->Branch("electron_mvaNoIsoID", &electron_mvaNoIsoID, "electron_mvaNoIsoID/I");
	} else if (requiredLeptonPDGID == 13) {
		tree->Branch("muon_pt", &muon_pt, "muon_pt/D");
		tree->Branch("muon_E", &muon_E, "muon_E/D");
		tree->Branch("muon_eta", &muon_eta, "muon_eta/D");
		tree->Branch("muon_phi", &muon_phi, "muon_phi/D");
		tree->Branch("muon_charge", &muon_charge, "muon_charge/I");
		tree->Branch("muon_trackIso", &muon_trackIso, "muon_trackIso/F");
		tree->Branch("muon_sumPuppiIso", &muon_sumPuppiIso, "muon_sumPuppiIso/F");
		tree->Branch("muon_sumPuppiNoLeptonIso", &muon_sumPuppiNoLeptonIso, "muon_sumPuppiNoLeptonIso/F");
		tree->Branch("muon_caloIso", &muon_caloIso, "muon_caloIso/F");
		tree->Branch("muon_CutBasedIdLoose", &muon_CutBasedIdLoose, "muon_CutBasedIdLoose/I");
		tree->Branch("muon_CutBasedIdMedium", &muon_CutBasedIdMedium, "muon_CutBasedIdMedium/I");
		tree->Branch("muon_CutBasedIdTight", &muon_CutBasedIdTight, "muon_CutBasedIdTight/I");
		tree->Branch("muon_CutBasedIdGlobalHighPt", &muon_CutBasedIdGlobalHighPt, "muon_CutBasedIdGlobalHighPt/I");
		tree->Branch("muon_CutBasedIdTrkHighPt", &muon_CutBasedIdTrkHighPt, "muon_CutBasedIdTrkHighPt/I");
		tree->Branch("muon_PFIsoLoose", &muon_PFIsoLoose, "muon_PFIsoLoose/I");
		tree->Branch("muon_PFIsoMedium", &muon_PFIsoMedium, "muon_PFIsoMedium/I");
		tree->Branch("muon_PFIsoTight", &muon_PFIsoTight, "muon_PFIsoTight/I");
		tree->Branch("muon_TkIsoLoose", &muon_TkIsoLoose, "muon_TkIsoLoose/I");
		tree->Branch("muon_MvaLoose", &muon_MvaLoose, "muon_MvaLoose/I");
		tree->Branch("muon_MvaMedium", &muon_MvaMedium, "muon_MvaMedium/I");
		tree->Branch("muon_MvaTight", &muon_MvaTight, "muon_MvaTight/I");
	} else if (requiredLeptonPDGID == 15) {
		// Tau kinematic parameters
		tree->Branch("tau2_pt",&tau2_pt,"tau2_pt/D");
		tree->Branch("tau2_eta",&tau2_eta,"tau2_eta/D");
		tree->Branch("tau2_phi",&tau2_phi,"tau2_phi/D");
		tree->Branch("tau2_q",&tau2_q,"tau2_q/D");
		tree->Branch("tau2_m",&tau2_m,"tau2_m/D");
		tree->Branch("tau2_dm",&tau2_dm,"tau2_dm/D");
		tree->Branch("tau2_dz",&tau2_dz,"tau2_dz/D");
		// Raw discriminators
		tree->Branch("tau2_absIso",&tau2_absIso,"tau2_absIso/D");
		tree->Branch("tau2_againstElectronRaw",&tau2_againstElectronRaw,"tau2_againstElectronRaw/D");
		tree->Branch("tau2_IsoMVArun2v1DBdR03oldDMwLTraw",&tau2_IsoMVArun2v1DBdR03oldDMwLTraw,"tau2_IsoMVArun2v1DBdR03oldDMwLTraw/D");
		tree->Branch("tau2_IsoMVArun2v1DBnewDMwLTraw",&tau2_IsoMVArun2v1DBnewDMwLTraw,"tau2_IsoMVArun2v1DBnewDMwLTraw/D");
		tree->Branch("tau2_IsoMVArun2v1DBoldDMwLTraw",&tau2_IsoMVArun2v1DBoldDMwLTraw,"tau2_IsoMVArun2v1DBoldDMwLTraw/D");
		tree->Branch("tau2_IsoMVArun2v1PWdR03oldDMwLTraw",&tau2_IsoMVArun2v1PWdR03oldDMwLTraw,"tau2_IsoMVArun2v1PWdR03oldDMwLTraw/D");
		tree->Branch("tau2_IsoMVArun2v1PWnewDMwLTraw",&tau2_IsoMVArun2v1PWnewDMwLTraw,"tau2_IsoMVArun2v1PWnewDMwLTraw/D");
		tree->Branch("tau2_IsoMVArun2v1PWoldDMwLTraw",&tau2_IsoMVArun2v1PWoldDMwLTraw,"tau2_IsoMVArun2v1PWoldDMwLTraw/D");

		// Discriminators
		tree->Branch("tau2_looseCombinedIso",&tau2_looseCombinedIso,"tau2_looseCombinedIso/I");
		tree->Branch("tau2_mediumCombinedIso",&tau2_mediumCombinedIso,"tau2_mediumCombinedIso/I");
		tree->Branch("tau2_tightCombinedIso",&tau2_tightCombinedIso,"tau2_tightCombinedIso/I");
		tree->Branch("tau2_VVlooseMvaIso",&tau2_VVlooseMvaIso,"tau2_VVlooseMvaIso/I");
		tree->Branch("tau2_VlooseMvaIso",&tau2_VlooseMvaIso,"tau2_VlooseMvaIso/I");
		tree->Branch("tau2_looseMvaIso",&tau2_looseMvaIso,"tau2_looseMvaIso/I");
		tree->Branch("tau2_mediumMvaIso",&tau2_mediumMvaIso,"tau2_mediumMvaIso/I");
		tree->Branch("tau2_tightMvaIso",&tau2_tightMvaIso,"tau2_tightMvaIso/I");
		tree->Branch("tau2_VtightMvaIso",&tau2_VtightMvaIso,"tau2_VtightMvaIso/I");
		tree->Branch("tau2_VVtightMvaIso",&tau2_VVtightMvaIso,"tau2_VVtightMvaIso/I");
		tree->Branch("tau2_looseMuonRejection",&tau2_looseMuonRejection,"tau2_looseMuonRejection/I");
		tree->Branch("tau2_tightMuonRejection",&tau2_tightMuonRejection,"tau2_tightMuonRejection/I");
		tree->Branch("tau2_looseElectronRejection",&tau2_looseElectronRejection,"tau2_looseElectronRejection/I");
		tree->Branch("tau2_mediumElectronRejection",&tau2_mediumElectronRejection,"tau2_mediumElectronRejection/I");
		tree->Branch("tau2_tightElectronRejection",&tau2_tightElectronRejection,"tau2_tightElectronRejection/I");
		tree->Branch("tau2_VtightElectronRejection",&tau2_VtightElectronRejection,"tau2_VtightElectronRejection/I");

		// Deep 2017v2
		if (DeepTau) {
			tree->Branch("tau2_Deep2017v2ElectronRejection",&tau2_Deep2017v2ElectronRejection,"tau2_Deep2017v2ElectronRejection/D");
			tree->Branch("tau2_Deep2017v2MuonRejection",&tau2_Deep2017v2MuonRejection,"tau2_Deep2017v2MuonRejection/D");
			tree->Branch("tau2_Deep2017v2JetRejection",&tau2_Deep2017v2JetRejection,"tau2_Deep2017v2JetRejection/D");

			tree->Branch("tau2_VVLooseDeepTau2017v2VSjet",&tau2_VVLooseDeepTau2017v2VSjet,"tau2_VVLooseDeepTau2017v2VSjet/I");
			tree->Branch("tau2_VLooseDeepTau2017v2VSjet",&tau2_VLooseDeepTau2017v2VSjet,"tau2_VLooseDeepTau2017v2VSjet/I");
			tree->Branch("tau2_LooseDeepTau2017v2VSjet",&tau2_LooseDeepTau2017v2VSjet,"tau2_LooseDeepTau2017v2VSjet/I");
			tree->Branch("tau2_MediumDeepTau2017v2VSjet",&tau2_MediumDeepTau2017v2VSjet,"tau2_MediumDeepTau2017v2VSjet/I");
			tree->Branch("tau2_TightDeepTau2017v2VSjet",&tau2_TightDeepTau2017v2VSjet,"tau2_TightDeepTau2017v2VSjet/I");
			tree->Branch("tau2_VTightDeepTau2017v2VSjet",&tau2_VTightDeepTau2017v2VSjet,"tau2_VTightDeepTau2017v2VSjet/I");
			tree->Branch("tau2_VVTightDeepTau2017v2VSjet",&tau2_VVTightDeepTau2017v2VSjet,"tau2_VVTightDeepTau2017v2VSjet/I");
			tree->Branch("tau2_LooseDeepTau2017v2VSmu",&tau2_LooseDeepTau2017v2VSmu,"tau2_LooseDeepTau2017v2VSmu/I");
			tree->Branch("tau2_MediumDeepTau2017v2VSmu",&tau2_MediumDeepTau2017v2VSmu,"tau2_MediumDeepTau2017v2VSmu/I");
			tree->Branch("tau2_TightDeepTau2017v2VSmu",&tau2_TightDeepTau2017v2VSmu,"tau2_TightDeepTau2017v2VSmu/I");

			tree->Branch("tau2_VVLooseDeepTau2017v2VSe",&tau2_VVLooseDeepTau2017v2VSe,"tau2_VVLooseDeepTau2017v2VSe/I");
			tree->Branch("tau2_VLooseDeepTau2017v2VSe",&tau2_VLooseDeepTau2017v2VSe,"tau2_VLooseDeepTau2017v2VSe/I");
			tree->Branch("tau2_LooseDeepTau2017v2VSe",&tau2_LooseDeepTau2017v2VSe,"tau2_LooseDeepTau2017v2VSe/I");
			tree->Branch("tau2_MediumDeepTau2017v2VSe",&tau2_MediumDeepTau2017v2VSe,"tau2_MediumDeepTau2017v2VSe/I");
			tree->Branch("tau2_TightDeepTau2017v2VSe",&tau2_TightDeepTau2017v2VSe,"tau2_TightDeepTau2017v2VSe/I");
			tree->Branch("tau2_VTightDeepTau2017v2VSe",&tau2_VTightDeepTau2017v2VSe,"tau2_VTightDeepTau2017v2VSe/I");
			tree->Branch("tau2_VVTightDeepTau2017v2VSe",&tau2_VVTightDeepTau2017v2VSe,"tau2_VVTightDeepTau2017v2VSe/I");
		}
	}
	// Tau Jet
	tree->Branch("TauJet_pt", &TauJet_pt, "TauJet_pt/D");
	tree->Branch("TauJet_eta", &TauJet_eta, "TauJet_eta/D");
	tree->Branch("TauJet_phi", &TauJet_phi, "TauJet_phi/D");
	tree->Branch("TauJet_E", &TauJet_E, "TauJet_E/D");

	// MET
	tree->Branch("met", &met, "met/D");
	tree->Branch("met_phi", &met_phi, "met_phi/D");
	tree->Branch("met_eta", &met_eta, "met_eta/D");
	tree->Branch("met_significance", &met_significance, "met_significance/D");
	tree->Branch("met_mEtSig", &met_mEtSig, "met_mEtSig/D");
	tree->Branch("met_energy", &met_energy, "met_energy/D");
	tree->Branch("tauMET_mass", &tauMET_mass, "tauMET_mass/D");
	tree->Branch("m_t", &m_t, "m_t/D");
	tree->Branch("dPhi", &dPhi, "dPhi/D");

	// Puppi MET
	tree->Branch("Puppimet", &Puppimet, "Puppimet/D");
	tree->Branch("Puppimet_phi", &Puppimet_phi, "Puppimet_phi/D");
	tree->Branch("Puppimet_eta", &Puppimet_eta, "Puppimet_eta/D");
	tree->Branch("Puppimet_significance", &Puppimet_significance, "Puppimet_significance/D");
	tree->Branch("Puppimet_metSig", &Puppimet_metSig, "Puppimet_metSig/D");
	tree->Branch("Puppimet_energy", &Puppimet_energy, "Puppimet_energy/D");
	tree->Branch("tauPuppimet_mass", &tauPuppimet_mass, "tauPuppimet_mass/D");
	tree->Branch("dPhiPuppimetTau", &dPhiPuppimetTau, "dPhiPuppimetTau/D");
	//tree->Branch("m_t", &m_t, "m_t/D");

	// True MET from MC
	tree->Branch("TrueMetPt", &TrueMetPt, "TrueMetPt/D");
	tree->Branch("TrueMetEta", &TrueMetEta, "TrueMetEta/D");
	tree->Branch("TrueMetPhi", &TrueMetPhi, "TrueMetPhi/D");
	tree->Branch("TrueMetEnergy", &TrueMetEnergy, "TrueMetEnergy/D");
	
	// Vertices, tracks, tau candidates
	tree->Branch("nVtx",&nVtx,"nVtx/I");
	tree->Branch("nTrks",&nTrks,"nTrks/I");
	tree->Branch("nTau",&nTau,"nTau/I");
	tree->Branch("nTauC",&nTauC,"nTauC/I");

	tree->Branch("nGamma",&nGamma,"nGamma/I");
	tree->Branch("nPhotons",&nPhotons,"nPhotons/I");

	// TauSpinner
	tree->Branch("WT", &WT, "WT/D");
	tree->Branch("WTFlip", &WTFlip, "WTFlip/D");
	tree->Branch("WThminus", &WThminus, "WThminus/D");
	tree->Branch("WThplus", &WThplus, "WThplus/D");
	tree->Branch("TauSpinnerMother", &TauSpinnerMother, "TauSpinnerMother/I");
	tree->Branch("nTauTriggers", &nTauTriggers, "nTauTriggers/I");
	tree->Branch("nJetHTTriggers", &nJetHTTriggers, "nJetHTTriggers/I");
	tree->Branch("nMETTriggers", &nMETTriggers, "nMETTriggers/I");
	tree->Branch("nBTagCSVTriggers", &nBTagCSVTriggers, "nBTagCSVTriggers/I");
	tree->Branch("nBTagMuTriggers", &nBTagMuTriggers, "nBTagMuTriggers/I");
	tree->Branch("nSingleMuonTriggers", &nSingleMuonTriggers, "nSingleMuonTriggers/I");
	tree->Branch("nSingleElectronTriggers", &nSingleElectronTriggers, "nSingleElectronTriggers/I");
	tree->Branch("nTargetTriggers", &nTargetTriggers, "nTargetTriggers/I");
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

	/////////////////////////////////// New method //////////////////////////////////////////
	/*
	const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);
    std::cout << "\n === TRIGGER PATHS === " << std::endl;
    for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) {
        std::cout << "Trigger " << names.triggerName(i) << 
            ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
            ": " << (triggerResults->accept(i) ? "PASS" : "fail (or not run)") 
            << std::endl;
    }
    std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
        obj.unpackPathNames(names);
        std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
        // Print trigger object collection and type
        std::cout << "\t   Collection: " << obj.collection() << std::endl;
        std::cout << "\t   Type IDs:   ";
        for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
        std::cout << std::endl;
        // Print associated trigger filters
        std::cout << "\t   Filters:    ";
        for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
        std::cout << std::endl;
        std::vector<std::string> pathNamesAll  = obj.pathNames(false);
        std::vector<std::string> pathNamesLast = obj.pathNames(true);
        // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
        // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
        // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
        std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
        for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
            bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
            bool isL3   = obj.hasPathName( pathNamesAll[h], false, true ); 
            bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
            bool isNone = obj.hasPathName( pathNamesAll[h], false, false ); 
            std::cout << "   " << pathNamesAll[h];
            if (isBoth) std::cout << "(L,3)";
            if (isL3 && !isBoth) std::cout << "(*,3)";
            if (isLF && !isBoth) std::cout << "(L,*)";
            if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    */

    /////////////////////////////////////////////////////////////////

	std::vector<std::string> trigNameVec;
	std::vector<bool> trigPassVec;
	std::vector<int> trigPsVec;
	std::vector<int> trigL1minPsVec;
	std::vector<int> trigL1maxPsVec;
	std::vector<int> trigPrescaleVec;
	int finalPrescale = 0;

	bool triggerOK = false;
	nTauTriggers = 0;
	nJetHTTriggers = 0;
	nMETTriggers = 0;
	nBTagCSVTriggers = 0;
	nBTagMuTriggers = 0;
	nSingleMuonTriggers = 0;
	nSingleElectronTriggers = 0;
	nTargetTriggers = 0;
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
			if (hlt > 0 && monitoringHLT) {
				std::cout << "Name    = " << trigName << std::endl;
		        std::cout << "ps      = " << ps << std::endl;
		        std::cout << "psL1min = " << psL1min << std::endl;
		        std::cout << "psL1max = " << psL1max << std::endl;
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
				for ( unsigned int i=0; i<trigNamesTau.size(); ++i ) {
					if ( triggerNames_[iHLT].find(trigNamesTau[i].c_str())!= std::string::npos ) {
						nTauTriggers++;
						if (monitoringHLT) std::cout << "Tau Trigger" << std::endl;
					}
				}
				for ( unsigned int i=0; i<trigNamesJetHT.size(); ++i ) {
					if ( triggerNames_[iHLT].find(trigNamesJetHT[i].c_str())!= std::string::npos ) {
						nJetHTTriggers++;
						if (monitoringHLT) std::cout << "JetHT Trigger" << std::endl;
					}
				}
				for ( unsigned int i=0; i<trigNamesMET.size(); ++i ) {
					if ( triggerNames_[iHLT].find(trigNamesMET[i].c_str())!= std::string::npos ) {
						nMETTriggers++;
						if (monitoringHLT) std::cout << "MET Trigger" << std::endl;
					}
				}
				for ( unsigned int i=0; i<trigNamesBTagCSV.size(); ++i ) {
					if ( triggerNames_[iHLT].find(trigNamesBTagCSV[i].c_str())!= std::string::npos ) {
						nBTagCSVTriggers++;
						if (monitoringHLT) std::cout << "BTagCSV Trigger" << std::endl;
					}
				}
				for ( unsigned int i=0; i<trigNamesBTagMu.size(); ++i ) {
					if ( triggerNames_[iHLT].find(trigNamesBTagMu[i].c_str())!= std::string::npos ) {
						nBTagMuTriggers++;
						if (monitoringHLT) std::cout << "BTagMu Trigger" << std::endl;
					}
				}
				for ( unsigned int i=0; i<trigNamesSingleMuon.size(); ++i ) {
					if ( triggerNames_[iHLT].find(trigNamesSingleMuon[i].c_str())!= std::string::npos ) {
						nSingleMuonTriggers++;
						if (monitoringHLT) std::cout << "# " << triggerNames_[iHLT] << std::endl;
					}
				}
				for ( unsigned int i=0; i<trigNamesSingleElectron.size(); ++i ) {
					if ( triggerNames_[iHLT].find(trigNamesSingleElectron[i].c_str())!= std::string::npos ) {
						nSingleElectronTriggers++;
						if (monitoringHLT) std::cout << "SingleElectron Trigger" << std::endl;
					}
				}
				for ( unsigned int i=0; i<trigNamesTarget.size(); ++i ) {
					if ( triggerNames_[iHLT].find(trigNamesTarget[i].c_str())!= std::string::npos ) {
						nTargetTriggers++;
						if (monitoringHLT) std::cout << "# " << triggerNames_[iHLT] << std::endl;
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
    	if (monitoringHLT) {
			std::cout << "###################" << std::endl;
			std::cout << "current Index    = " << l << std::endl;
			std::cout << "current Name     = " << trigNameVec[l] << std::endl;
			std::cout << "current Prescale = " << trigPrescaleVec[l] << std::endl;
    	}
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

	if (nTauTriggers + nJetHTTriggers + nMETTriggers + nBTagCSVTriggers + nSingleMuonTriggers + nSingleElectronTriggers > 0 && !useTargetHLT) {
		triggerOK = true;
	} else if (nTargetTriggers > 0 && useTargetHLT) {
		triggerOK = true;
	}
	if (useHLT) return triggerOK;
	else return true;
}

//-----------------------------------------------------------------------------------------

bool TTbarTauLepton::AddTau(const edm::Event& event) {

	if (monitoringTau) std::cout << std::endl << "!!! Tau !!!" << std::endl;

	tau_found = 0;
	nTauC     = 0;
	math::XYZTLorentzVector DiPhoton_p4;
	int nPhotons_temp = 0;

	edm::Handle<pat::TauCollection> taus;
	event.getByToken(TauCollectionToken_, taus);
	if (!taus.isValid() || taus->size() == 0) return false;

	edm::Handle<reco::VertexCompositePtrCandidateCollection> SecondaryVertices;
	event.getByToken(SVToken_, SecondaryVertices);

	// search for the tau candidate with the minimum isolation and the
	// maximum transverse momentum
	size_t index  = taus->size();;
	//int index_2 = null;
	//std::map <int, std::pair <double, double>> TauIdMap;

	// Tau index and pass of cuts
	std::map <int, bool> TauIndexMap;
	for (size_t i = 0; i < taus->size(); ++i) {
#define cut(condition) if (!(condition)) continue;
		TauIndexMap.insert(std::pair<int, bool>(i, false));
		auto& tau = (*taus)[i];
		if (monitoringTau) {
			std::cout << "######### Tau candidate number " << i << " ################" << std::endl;
			std::cout << "pt = " << tau.pt() << ", eta = " << tau.eta() << ", "
			<< "IsoIdMVA = " << tau.tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw") << ", "
			<< "Jet = " << tau.tauID("byVVLooseIsolationMVArun2v1DBoldDMwLT") << ", "
			<< "Ele = " << tau.tauID("againstElectronVLooseMVA6") << ", "
			<< "Mu = " << tau.tauID("againstMuonLoose3") << ", "
			<< "dZ = " << (pv_position - tau.vertex()).R() << std::endl;
		}
		cut(tau.pt() > tauPtMin); // tau transverse momentum
		cut(TMath::Abs(tau.eta()) < tauEtaMax); // tau pseudorapidity
		if (looseTauID && !DeepTau) {
			cut(tau.tauID("byVVLooseIsolationMVArun2v1DBoldDMwLT")); // at least VVloose Iso MVA
			cut(tau.tauID("againstElectronVLooseMVA6")); // at least loose Iso MVA
			cut(tau.tauID("againstMuonLoose3")); // at least loose Iso MVA
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
		if (monitoringTau) std::cout << "PV passed" << std::endl;
		if (monitoringTau) std::cout << "candidate passed cuts" << std::endl;
		TauIndexMap.erase(i);
		TauIndexMap.insert(std::pair<int, bool>(i, true));
		if (monitoringTau) std::cout << "Index map pair inseted" << std::endl;

		if (index < taus->size()) {
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


	if (monitoringTau) {
		std::cout << "Tau index (size) = " << index << "(" << taus->size() << ")" << std::endl;
		//for (auto iter = TauIndexMap.begin(); iter != TauIndexMap.end(); ++iter) {
		//	std::cout << "tau index = " << iter->first << ", value = " << iter->second << std::endl;
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
	tau_VVlooseMvaIso           = tau.tauID("byVVLooseIsolationMVArun2v1DBoldDMwLT");
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
	tau_looseMuonRejection      = tau.tauID("againstMuonLoose3");
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


void TTbarTauLepton::AddEmptyTau(const edm::Event& event) {
	tau_pt   = null;
	tau_eta  = null;
	tau_phi  = null;
	tau_dm   = null;
	tau_dz   = null;
	tau_q    = null;
	tau_m    = null;
	
	tau_absIso                  = null;
	tau_looseCombinedIso        = null;
	tau_mediumCombinedIso       = null;
	tau_tightCombinedIso        = null;
	tau_VVlooseMvaIso           = null;
	tau_VlooseMvaIso            = null;
	tau_looseMvaIso             = null;
	tau_mediumMvaIso            = null;
	tau_tightMvaIso             = null;
	tau_VtightMvaIso            = null;
	tau_VVtightMvaIso           = null;
	// New Raw discriminators
	tau_againstElectronRaw            = null;
	tau_IsoMVArun2v1DBdR03oldDMwLTraw = null;
	tau_IsoMVArun2v1DBnewDMwLTraw     = null;
	tau_IsoMVArun2v1DBoldDMwLTraw     = null;
	tau_IsoMVArun2v1PWdR03oldDMwLTraw = null;
	tau_IsoMVArun2v1PWnewDMwLTraw     = null;
	tau_IsoMVArun2v1PWoldDMwLTraw     = null;
	//
	tau_looseMuonRejection      = null;
	tau_tightMuonRejection      = null;
	tau_looseElectronRejection  = null;
	tau_mediumElectronRejection = null;
	tau_tightElectronRejection  = null;
	tau_VtightElectronRejection = null;
	decayModeFindingNewDMs      = null;
	decayModeFinding            = null;

	tau_Deep2017v2ElectronRejection = null;
	tau_Deep2017v2MuonRejection     = null;
	tau_Deep2017v2JetRejection      = null;
	
	tau_VVLooseDeepTau2017v2VSjet   = null;
	tau_VLooseDeepTau2017v2VSjet    = null;
	tau_LooseDeepTau2017v2VSjet     = null;
    tau_MediumDeepTau2017v2VSjet    = null;
    tau_TightDeepTau2017v2VSjet     = null;
    tau_VTightDeepTau2017v2VSjet    = null;
    tau_VVTightDeepTau2017v2VSjet   = null;

    tau_LooseDeepTau2017v2VSmu      = null;
    tau_MediumDeepTau2017v2VSmu     = null;
    tau_TightDeepTau2017v2VSmu      = null;

    tau_VVLooseDeepTau2017v2VSe     = null;
    tau_VLooseDeepTau2017v2VSe      = null;
    tau_LooseDeepTau2017v2VSe       = null;
    tau_MediumDeepTau2017v2VSe      = null;
    tau_TightDeepTau2017v2VSe       = null;
    tau_VTightDeepTau2017v2VSe      = null;
    tau_VVTightDeepTau2017v2VSe     = null;

	piChar_pt = null;
	piChar_eta = null;
	piChar_phi = null;
	piChar_m = null;
	piChar_q = null;

    piZero_pt  = null;
	piZero_eta = null;
	piZero_phi = null;
	piZero_m   = null;

	DiPhoton_pt  = null;
	DiPhoton_eta = null;
	DiPhoton_phi = null;
	DiPhoton_m   = null;

	TauPhoton1_pt      = null;
	TauPhoton1_eta     = null;
	TauPhoton1_phi     = null;
	TauPhoton1_energy  = null;
	TauPhoton2_pt      = null;
	TauPhoton2_eta     = null;
	TauPhoton2_phi     = null;
	TauPhoton2_energy  = null;

	tau_hasSV     = null;
	tau_SVdR      = null;
	SV_Chi2       = null;
	SV_Chi2NDF    = null;
	tauPVtoSVdPhi = null;
	tau_dR2       = null;
}

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
	// test
	//double dRmin = null;
	double dRmin       = 100;
	double leptondRmin = 100;
	double b1dRmin     = 100;
	double b2dRmin     = 100;

	math::XYZTLorentzVector TauVisP4;
	math::XYZTLorentzVector SumNuP4;

	std::vector <reco::GenParticle> GenTauCandidates;

	// Look for tau among generated particles
	for (auto& particle: *genParticles) {
#define cut(condition) if (!(condition)) continue;
		// look for the tau -> pi+ pi0 neutrino decay most oriented towards
		// the reconstructed tau (if present)
		cut(abs(particle.pdgId()) == pdg_tau);
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

		if (!fullMC) {
			cut(particle.pt() > tauPtMin);
			cut(TMath::Abs(particle.eta()) < tauEtaMax);
			cut((pv_position - particle.vertex()).R() < tauDzMax);
		}

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
		} else if (fullMC) {
			if (monitoringGen) {
				std::cout << "GenTau candidate" << std::endl << "daughters: ";
				for (unsigned i = 0; i < particle.numberOfDaughters(); ++i) {
					const reco::Candidate* daughter = particle.daughter(i);
					std::cout << daughter->pdgId() << ", ";
				}
				std::cout << std::endl << "pt = " << particle.pt() << std::endl;
			}
			GenTauCandidates.push_back(particle);
		}
#undef cut
		//if (monitoringGen) GenTauMonitor(particle);
	};

	if (fullMC) {
		SortGenTaus(GenTauCandidates);
		for (unsigned n = 0; n < GenTauCandidates.size(); n++) {
			double lepton1dR = sqrt(deltaR2(GenTauCandidates[n].eta(), GenTauCandidates[n].phi(), lepton1_eta, lepton1_phi));
			if (lepton1dR < 0.1) {
				continue;
			} else {
				tau = &GenTauCandidates[n];
				break;
			}
		}
		// Selected tau candidate
		if (monitoringGen && tau) {
			std::cout << "GenTau found" << std::endl << "daughters: ";
			for (unsigned i = 0; i < tau->numberOfDaughters(); ++i) {
				const reco::Candidate* daughter = tau->daughter(i);
				std::cout << daughter->pdgId() << ", ";
			}
			std::cout << std::endl << "pt = " << tau->pt() << std::endl;
		}
	}

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
		if (abs(particle.pdgId()) == pdg_mu || abs(particle.pdgId()) == pdg_electron || abs(particle.pdgId()) == pdg_tau) {
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
			double bquark1dR_ = sqrt(deltaR2(particle.eta(), particle.phi(), BJet1_eta, BJet1_phi));
			if ( !bquark1 || bquark1dR_ < b1dRmin) {
				bquark1 = &particle;
				b1dRmin = bquark1dR_;
				bquark1dR = b1dRmin;
				if (monitoringGen) {
					std::cout << "---------- b-quark 1" << std::endl;
					std::cout << "dR(calc) = " << bquark1dR << std::endl;
					GenRecoMonitor *b1Gen = new GenRecoMonitor(particle, pdg_bquark, BJet1_pt, BJet1_eta, BJet1_phi);
					b1Gen->PrintComp(true, true);
					delete b1Gen;
				}
			}
		} else continue;
	}

	for (auto& particle: *genParticles) {
		if (abs(particle.pdgId()) == pdg_bquark) {
			double bquark2dR_ = sqrt(deltaR2(particle.eta(), particle.phi(), BJet2_eta, BJet2_phi));
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
		} else continue;
	}

	if (bquark1) {
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

	//return;
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

	if (met < METcut) return false;
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
	int nPuppiJets30PVpassed = 0;
	nPuppiJets20PV    = 0;
	nLooseBtagedPuppiJets    = 0;
	nMediumBtagedPuppiJets   = 0;
	nTightBtagedPuppiJets    = 0;
	nLooseBtagedPuppiJetsPV  = 0;
	nMediumBtagedPuppiJetsPV = 0;
	nTightBtagedPuppiJetsPV  = 0;
	unsigned NrequiredBJets_unsigned = 0;
	if (NrequiredBJets > 0) NrequiredBJets_unsigned = NrequiredBJets;

	LeadingJet_pt   = null;
	LeadingJet_eta  = null;
	LeadingJet_phi  = null;
	LeadingJet_m    = null;
	LeadingJet_btag = null;
	SubLeadingJet_pt   = null;
	SubLeadingJet_eta  = null;
	SubLeadingJet_phi  = null;
	SubLeadingJet_m    = null;
	SubLeadingJet_btag = null;

	ttbarEvent  = null;
	BJet1_pt    = null;
	BJet1_eta   = null;
	BJet1_phi   = null;
	BJet1_bprob = null;
	BJet1_E     = null;
	BJet2_pt    = null;
	BJet2_eta   = null;
	BJet2_phi   = null;
	BJet2_bprob = null;
	BJet2_E     = null;

	TauJet_pt = null;
	TauJet_eta = null;
	TauJet_phi = null;
	TauJet_E = null;

	int nLooseBJets = 0;

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

	// Selection of two bJets for ttbar analysis
	// Also study of multijets and leading/subleading jets in event

	// Vector of bJet candidates
	std::vector <BJetCandidate> looseBJets;

	edm::Handle<pat::JetCollection> Puppijets;
	event.getByToken(PuppiJetCollectionToken_, Puppijets);
	if (!Puppijets.isValid()) return false;

	int j = 0;
	// Lopp over Jets
	for (auto& jet: *Puppijets) {
		if (TMath::Abs(jet.eta()) > JetEtaMax) continue;
		///if (TMath::Abs(jet.eta()) > 3.) continue;
		if (jet.pt() > 20) {///
		///if (jet.pt() > 15) {
			j++;
			if (monitoringJets) {
				std::cout << std::endl;
				std::cout << "PuppiJet[" << j << "] Pt = " << jet.pt() << std::endl;
			}
			// Search for tau jet
			double deltaRTauJet = sqrt(deltaR2(tau_eta, tau_phi, jet.eta(), jet.phi()));
			if (deltaRTauJet < 0.1) {
				if (monitoringJets) std::cout << "This is Tau jet, dR = " << deltaRTauJet << std::endl;
				TauJet_pt = jet.pt();
				TauJet_eta = jet.eta();
				TauJet_phi = jet.phi();
				TauJet_E = jet.energy();
				//continue; 
			}
			// Exclude jets close to lepton 1
			double deltaRLepton1 = sqrt(deltaR2(jet.eta(), jet.phi(), lepton1_eta, lepton1_phi));
			if (deltaRLepton1 < 0.1) {
				if (monitoringJets) std::cout << "Jet with leton inside cond dR < 0.3" << std::endl;
				///continue; ///
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

			// Fill vector of jets with at least loose btag
			if ((jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_loose) &&
				(jet.pt() > BJetPtMin) && (TMath::Abs(jet.eta()) < JetEtaMax)) {
				///if ((jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_loose)) {
				// BJetFound = true;
				///if (TMath::Abs(jet.eta()) > EtaMax) continue;
				///if (jet.pt() < BJetPtMin) continue;
				double deltaRTauJet = sqrt(deltaR2(tau_eta, tau_phi, jet.eta(), jet.phi()));
				if (deltaRTauJet < 0.1) {
					nTauJet++;
					///continue;
				}
				nLooseBJets++;
				BJetCandidate BJet(jet, JetFromPV);
				looseBJets.push_back(BJet);
			}

    		if (j == 1) {
    		    LeadingJet_pt   = jet.pt();
    		    LeadingJet_eta  = jet.eta();
    		    LeadingJet_phi  = jet.phi();
    		    LeadingJet_m    = jet.mass();
    		    LeadingJet_btag = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb");
    		} else if (j == 2) {
    		    SubLeadingJet_pt   = jet.pt();
    		    SubLeadingJet_eta  = jet.eta();
    		    SubLeadingJet_phi  = jet.phi();
    		    SubLeadingJet_m    = jet.mass();
    		    SubLeadingJet_btag = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb");
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
					if (jet.neutralHadronEnergy() / jet.energy() > 0.9) continue;
					if (jet.neutralEmEnergyFraction() > 0.9) continue;
					//if (jet.getPFConstituents().size() <= 1) continue;
					if (jet.muonEnergyFraction() > 0.8) continue;
					if (jet.chargedHadronEnergyFraction() < 0.001) continue;
					if (jet.chargedMultiplicity() <= 0) continue;
					if (jet.chargedEmEnergyFraction() > 0.8) continue;
					nPuppiJets30PVpassed++;
				}
			}
		}
	}
	if (monitoringBJets) std::cout << "number of 30 GeV Jets from PV (passed selections) = " << nPuppiJets30PV
	<< "(" << nPuppiJets30PV << ")" << std::endl;

	if (monitoringBJets) std::cout << std::endl << "!!! bTagged Jets !!!" << std::endl;

	SortJets(looseBJets);

	if (monitoringBJets) std::cout << "List of loose bJets (" << looseBJets.size() << ", " << nLooseBJets << "):" << std::endl;
	for(unsigned ijet = 0; ijet < looseBJets.size(); ijet++) {
		if (monitoringBJets) {
			std::cout << "jet " << ijet << std::endl;
			looseBJets[ijet].Monitoring();
		}
	}

	if (looseBJets.size() < NrequiredBJets_unsigned) {
		std::cout << "Loose b-jets required " << NrequiredBJets_unsigned << ", found " << looseBJets.size() << std::endl;
		return false;
	}

	if (looseBJets.size() > 1) {
		ttbarEvent  = 1;
		BJet1_pt    = looseBJets[0].Pt;
		BJet1_eta   = looseBJets[0].Eta;
		BJet1_phi   = looseBJets[0].Phi;
		BJet1_bprob = looseBJets[0].Bprob;
		BJet1_E     = looseBJets[0].FourMomentum.E();
		BJet2_pt    = looseBJets[1].Pt;
		BJet2_eta   = looseBJets[1].Eta;
		BJet2_phi   = looseBJets[1].Phi;
		BJet2_bprob = looseBJets[1].Bprob;
		BJet2_E     = looseBJets[1].FourMomentum.E();
    } else if (looseBJets.size() > 0) {
    	BJet1_pt    = looseBJets[0].Pt;
		BJet1_eta   = looseBJets[0].Eta;
		BJet1_phi   = looseBJets[0].Phi;
		BJet1_bprob = looseBJets[0].Bprob;
		BJet1_E     = looseBJets[0].FourMomentum.E();
    }

	looseBJets.clear();

	return nPuppiJets20 > 0;
};

bool TTbarTauLepton::AddLepton (const edm::Event& event) {

	if (monitoringLeptons) std::cout << std::endl << "!!! Leptons !!!" << std::endl;

	nLeptonCandidates = 0;
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
	lepton1_tauAbsIso = null;

	electron_pt = null;
	electron_eta = null;
	electron_phi = null;
	electron_E = null;
	electron_charge = null;
	electron_trackIso = null;
	electron_caloIso = null;
	electron_sumPuppiIso = null;
	electron_sumPuppiNoLeptonIso = null;
	electron_cutBasedID_loose    = null;
	electron_cutBasedID_medium   = null;
	electron_cutBasedID_tight    = null;
	//electron_mvaIsoID_loose      = null;
	//electron_mvaIsoID_wp80       = null;
	//electron_mvaIsoID_wp90       = null;
	electron_mvaNoIsoID_loose    = null;
	electron_mvaNoIsoID_wp80     = null;
	electron_mvaNoIsoID_wp90     = null;
	//electron_cutBasedID          = null;
	//electron_mvaIsoID            = null;
	//electron_mvaNoIsoID          = null;

	muon_pt = null;
	muon_eta = null;
	muon_phi = null;
	muon_E = null;
	muon_charge = null;
	muon_trackIso = null;
	muon_sumPuppiIso = null;
	muon_sumPuppiNoLeptonIso = null;
	muon_caloIso = null;
	muon_CutBasedIdLoose = null;
	muon_CutBasedIdMedium = null;
	muon_CutBasedIdTight = null;
	muon_CutBasedIdGlobalHighPt = null;
	muon_CutBasedIdTrkHighPt = null;
	muon_PFIsoLoose = null;
	muon_PFIsoMedium = null;
	muon_PFIsoTight = null;
	muon_TkIsoLoose = null;
	muon_TkIsoTight = null;
	muon_MvaLoose = null;
	muon_MvaMedium = null;
	muon_MvaTight = null;

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
	int nLeptons = 0;

	// Select electrons and muons which satisfy following criteria
	// Not in bJets, pt above threshold, eta below threshold
	std::vector<pat::Electron> ElectronCandidates;
	if (electrons.isValid() && !(*electrons).empty()) {
		//
		if (monitoringLeptons) std::cout << "Isolated electrons:" << std::endl;
		for (auto& electron: *electrons) {
			nLeptons++;
			if (monitoringLeptons) std::cout << "Electron " << nLeptons << std::endl;
#define cut(condition) if (!(condition)) continue;
			cut(electron.pt() > MuElePtMin);
			// Discard events with electron with pt > 15 GeV
			if (NrequiredLeptons < 0) return false;
			cut(abs(electron.eta()) < EtaMax);
			reco::TrackRef trackRef;
			trackRef = electron.track();
			if (trackRef.isNonnull()) {
				if (abs(electron.eta()) < 1.3) {
					cut(std::abs(trackRef->dxy(pv_position)) <= 0.05);
					cut(std::abs(trackRef->dz(pv_position)) <= 0.1);
				} else {
					cut(std::abs(trackRef->dxy(pv_position)) <= 0.1);
					cut(std::abs(trackRef->dz(pv_position)) <= 0.2);
				}
			}
			cut(electron.electronID("cutBasedElectronID-Fall17-94X-V1-loose")); // loose for a start (MVA id is also availible)
			//cut(electron.electronID("cutBasedElectronID-Fall17-94X-V1-tight"));
			//
			if (monitoringLeptons) std::cout << "Passed cuts" << std::endl;
			//double deltaRJet1 = sqrt(deltaR2(BJet1_eta, BJet1_phi, electron.eta(), electron.phi()));
			//double deltaRJet2 = sqrt(deltaR2(BJet2_eta, BJet2_phi, electron.eta(), electron.phi()));
			double deltaRTau  = sqrt(deltaR2(tau_eta, tau_phi, electron.eta(), electron.phi()));

			if (deltaRTau < 0.3) {
				if (monitoringLeptons) std::cout << "In Tau cone: dR = " << deltaRTau << std::endl;
				return false;
			}
			//if (deltaRJet1 > 0.4 && deltaRJet2 > 0.4) {
			if (deltaRTau > 0.3) {
				bool SameLepton = false;
				if(!LepCandidates.empty()) {
					for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
						if (abs(LepCandidates[ilep].Phi - electron.phi()) < 0.3) {
							if (monitoringLeptons) std::cout << "Same lepton" << std::endl;
							SameLepton = true;
						}
					}
				}
				cut(!SameLepton);
				nLeptonCandidates++;
				LeptonCandidate Lepton(electron, pv_position);
				LepCandidates.push_back(Lepton);
				ElectronCandidates.push_back(electron);
				if (monitoringLeptons) std::cout << "Added to vector" << std::endl;
	//delete Lepton;
			}
#undef cut
		}
	}

	// Muons, minor cuts commented
	// declare the vector of muons
	std::vector<pat::Muon> MuonCandidates;
	if (muons.isValid() && !(*muons).empty()) {
		if (monitoringLeptons) std::cout << "Isolated muons:" << std::endl;
		for (auto& muon: *muons) {
			nLeptons++;
			if (monitoringLeptons) std::cout << "Muon " << nLeptons << std::endl;
#define cut(condition) if (!(condition)) continue;
			cut(muon.pt() > MuElePtMin);
			// Discard events with muon with pt > 15 GeV
			if (NrequiredLeptons < 0) return false;
			cut(abs(muon.eta()) < EtaMax);
			cut(muon.isPFMuon());
			cut(muon.isGlobalMuon());
			cut(muon.numberOfMatchedStations() > 1);
			// check muon dXY and dZ
			reco::TrackRef trackRef;
			trackRef = muon.innerTrack();
			if (trackRef.isNonnull()) {
				cut(std::abs(trackRef->dxy(pv_position)) <= 0.05);
				cut(std::abs(trackRef->dz(pv_position)) <= 0.1);
				//cut(trackRef->found() > 5);  // more than 5 hits in inner tracker
			}
			//cut(muon.passed(reco::Muon::PFIsoLoose)); // I_rel < 0.25 (particle flow)
			if (monitoringLeptons) std::cout << "Passed cuts" << std::endl;
			double deltaRTau  = sqrt(deltaR2(tau_eta, tau_phi, muon.eta(), muon.phi()));
			if (deltaRTau < 0.3) {
				if (monitoringLeptons) std::cout << "In Tau cone: dR = " << deltaRTau << std::endl;
				return false;
			}
			if (deltaRTau > 0.3) {
				bool SameLepton = false;
				if(!LepCandidates.empty()) {
					for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
						if (abs(LepCandidates[ilep].Phi - muon.phi()) < 0.3) {
							if (monitoringLeptons) std::cout << "Same lepton" << std::endl;
							SameLepton = true;
						}
					}
				}
				cut(!SameLepton);
				LeptonCandidate Lepton(muon, pv_position);
				nLeptonCandidates++;
				LepCandidates.push_back(Lepton);
				MuonCandidates.push_back(muon);
				if (monitoringLeptons) std::cout << "Added to vector" << std::endl;
	//delete Lepton;
			}
#undef cut
		}
	}

	// tau as second lepton candidate
	std::vector<pat::Tau> TauCandidates;
	if (taus.isValid() && !(*taus).empty()) {
		//
		// declare the vector of muons
		if (monitoringLeptons) std::cout << "Isolated taus:" << std::endl;
		for (auto& tau: *taus) {
			nLeptons++;
			if (monitoringLeptons) std::cout << "Tau " << nLeptons << std::endl;
#define cut(condition) if (!(condition)) continue;
			cut(tau.pt() > MuElePtMin);
			cut(tau.tauID("byVVLooseIsolationMVArun2v1DBoldDMwLT"));
			cut(abs(tau.eta()) < EtaMax);
			cut((pv_position - tau.vertex()).R() < tauDzMax);
			if (monitoringLeptons) std::cout << "Passed cuts" << std::endl;
			//double deltaRJet1 = sqrt(deltaR2(BJet1_eta, BJet1_phi, tau.eta(), tau.phi()));
			//double deltaRJet2 = sqrt(deltaR2(BJet2_eta, BJet2_phi, tau.eta(), tau.phi()));
			double deltaRTau  = sqrt(deltaR2(tau_eta, tau_phi, tau.eta(), tau.phi()));
			cut(deltaRTau > 0.3);
			if (deltaRTau > 0.3) {
				bool SameLepton = false;
				if(!LepCandidates.empty()) {
					for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
						if (abs(LepCandidates[ilep].Phi - tau.phi()) < 0.3) {
							if (monitoringLeptons) std::cout << "Same lepton" << std::endl;
							SameLepton = true;
						}
					}
				}
				cut(!SameLepton);
				nLeptonCandidates++;
				LeptonCandidate Lepton(tau, pv_position);
				LepCandidates.push_back(Lepton);
				TauCandidates.push_back(tau);
				if (monitoringLeptons) std::cout << "Added to vector" << std::endl;
	//delete Lepton;
			}
#undef cut
		}
	}

	if (nLeptonCandidates < NrequiredLeptons) {
		if (monitoringLeptons) std::cout << "Number of lepton candiates = " << LepCandidates.size() << ", " << nLeptonCandidates << std::endl;
		return false;
	}

	SortLeptons(LepCandidates);

	if (monitoringLeptons) std::cout << "List of selected leptons:" << std::endl;
	for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
		if (monitoringLeptons) {
			std::cout << "lepton " << ilep << std::endl;
			LepCandidates[ilep].Monitoring();
			//std::cout << "delta(BJet1) = " << sqrt(deltaR2(BJet1_eta, BJet1_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
			//std::cout << "delta(BJet2) = " << sqrt(deltaR2(BJet2_eta, BJet2_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
			std::cout << "delta(Tau)   = " << sqrt(deltaR2(tau_eta, tau_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
		}
	}

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
		lepton1_tauAbsIso = LepCandidates[0].tauAbsIso;
	}

	LepCandidates.clear();

	if (requiredLeptonPDGID == 11 && ElectronCandidates.size() > 0) {
		SortElectrons(ElectronCandidates);
		if (monitoringLeptons) printElectronCandidate(ElectronCandidates[0]);
		electron_pt = ElectronCandidates[0].pt();
		electron_eta = ElectronCandidates[0].eta();
		electron_phi = ElectronCandidates[0].phi();
		electron_E = ElectronCandidates[0].energy();
		electron_charge = ElectronCandidates[0].charge();
		electron_trackIso = ElectronCandidates[0].trackIso();
		electron_caloIso = ElectronCandidates[0].caloIso();
		electron_sumPuppiIso = ElectronCandidates[0].puppiChargedHadronIso() + ElectronCandidates[0].puppiNeutralHadronIso() + ElectronCandidates[0].puppiPhotonIso();
		electron_sumPuppiNoLeptonIso = ElectronCandidates[0].puppiNoLeptonsChargedHadronIso() + ElectronCandidates[0].puppiNoLeptonsNeutralHadronIso() + ElectronCandidates[0].puppiNoLeptonsPhotonIso();
		electron_cutBasedID_loose    = ElectronCandidates[0].electronID("cutBasedElectronID-Fall17-94X-V1-loose");
		electron_cutBasedID_medium   = ElectronCandidates[0].electronID("cutBasedElectronID-Fall17-94X-V1-medium");
		electron_cutBasedID_tight    = ElectronCandidates[0].electronID("cutBasedElectronID-Fall17-94X-V1-tight");
		//electron_mvaIsoID_loose      = ElectronCandidates[0].electronID("mvaEleID-Fall17-iso-V1-wpLoose");
		//electron_mvaIsoID_wp80       = ElectronCandidates[0].electronID("mvaEleID-Fall17-iso-V1-wp80");
		//electron_mvaIsoID_wp90       = ElectronCandidates[0].electronID("mvaEleID-Fall17-iso-V1-wp90");
		electron_mvaNoIsoID_loose    = ElectronCandidates[0].electronID("mvaEleID-Fall17-noIso-V1-wpLoose");
		electron_mvaNoIsoID_wp80     = ElectronCandidates[0].electronID("mvaEleID-Fall17-noIso-V1-wp80");
		electron_mvaNoIsoID_wp90     = ElectronCandidates[0].electronID("mvaEleID-Fall17-noIso-V1-wp90");
		//electron_cutBasedID          = electron_cutBasedID_loose + electron_cutBasedID_medium + electron_cutBasedID_tight;
		//electron_mvaIsoID            = electron_mvaIsoID_loose + electron_mvaIsoID_wp90 + electron_mvaIsoID_wp90;
		//electron_mvaNoIsoID          = electron_mvaNoIsoID_loose + electron_mvaNoIsoID_wp80 + electron_mvaNoIsoID_wp90;
		//
		LeptonFound = true;
	} else if (requiredLeptonPDGID == 13 && MuonCandidates.size() > 0) {
		SortMuons(MuonCandidates);
		if (monitoringLeptons) printMuonCandidate(MuonCandidates[0]);
		muon_pt     = MuonCandidates[0].pt();
		muon_eta    = MuonCandidates[0].eta();
		muon_phi    = MuonCandidates[0].phi();
		muon_E      = MuonCandidates[0].energy();
		muon_charge = MuonCandidates[0].charge();
		muon_trackIso            = MuonCandidates[0].trackIso();
		muon_sumPuppiIso         = MuonCandidates[0].puppiChargedHadronIso() + MuonCandidates[0].puppiNeutralHadronIso() + MuonCandidates[0].puppiPhotonIso();
		muon_sumPuppiNoLeptonIso = MuonCandidates[0].puppiNoLeptonsChargedHadronIso() + MuonCandidates[0].puppiNoLeptonsNeutralHadronIso() + MuonCandidates[0].puppiNoLeptonsPhotonIso();
		muon_caloIso             = MuonCandidates[0].caloIso();
		muon_CutBasedIdLoose = MuonCandidates[0].passed(reco::Muon::CutBasedIdLoose);
		muon_CutBasedIdMedium = MuonCandidates[0].passed(reco::Muon::CutBasedIdMedium);
		muon_CutBasedIdTight = MuonCandidates[0].passed(reco::Muon::CutBasedIdTight);
		muon_CutBasedIdGlobalHighPt = MuonCandidates[0].passed(reco::Muon::CutBasedIdGlobalHighPt);
		muon_CutBasedIdTrkHighPt = MuonCandidates[0].passed(reco::Muon::CutBasedIdTrkHighPt);
		muon_PFIsoLoose = MuonCandidates[0].passed(reco::Muon::PFIsoLoose);
		muon_PFIsoMedium = MuonCandidates[0].passed(reco::Muon::PFIsoMedium);
		muon_PFIsoTight = MuonCandidates[0].passed(reco::Muon::PFIsoTight);
		muon_TkIsoLoose = MuonCandidates[0].passed(reco::Muon::TkIsoLoose);
		muon_TkIsoTight = MuonCandidates[0].passed(reco::Muon::TkIsoTight);
		muon_MvaLoose = MuonCandidates[0].passed(reco::Muon::MvaLoose);
		muon_MvaMedium = MuonCandidates[0].passed(reco::Muon::MvaMedium);
		muon_MvaTight = MuonCandidates[0].passed(reco::Muon::MvaTight);
		//
		LeptonFound = true;
	} else if (requiredLeptonPDGID == 15 && TauCandidates.size() > 0) {
		SortTaus(TauCandidates);
		if (monitoringLeptons) TauMonitor (TauCandidates[0], DeepTau, pv_position);
		tau2_pt                     = TauCandidates[0].pt();
		tau2_eta                    = TauCandidates[0].eta();
		tau2_phi                    = TauCandidates[0].phi();
		tau2_dm                     = TauCandidates[0].decayMode();
		tau2_q                      = TauCandidates[0].charge();
		tau2_m                      = TauCandidates[0].mass();
		
		tau2_absIso                  = TauCandidates[0].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
		tau2_looseCombinedIso        = TauCandidates[0].tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
		tau2_mediumCombinedIso       = TauCandidates[0].tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
		tau2_tightCombinedIso        = TauCandidates[0].tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
		tau2_VVlooseMvaIso           = TauCandidates[0].tauID("byVVLooseIsolationMVArun2v1DBoldDMwLT");
		tau2_VlooseMvaIso            = TauCandidates[0].tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
		tau2_looseMvaIso             = TauCandidates[0].tauID("byLooseIsolationMVArun2v1DBoldDMwLT"); // there are more possible discriminators 
		tau2_mediumMvaIso            = TauCandidates[0].tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
		tau2_tightMvaIso             = TauCandidates[0].tauID("byTightIsolationMVArun2v1DBoldDMwLT");
		tau2_VtightMvaIso            = TauCandidates[0].tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
		tau2_VVtightMvaIso           = TauCandidates[0].tauID("byVVTightIsolationMVArun2v1DBoldDMwLT");
		// New Raw discriminators
		tau2_againstElectronRaw            = TauCandidates[0].tauID("againstElectronMVA6Raw");
		tau2_IsoMVArun2v1DBdR03oldDMwLTraw = TauCandidates[0].tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw");
		tau2_IsoMVArun2v1DBnewDMwLTraw     = TauCandidates[0].tauID("byIsolationMVArun2v1DBnewDMwLTraw");
		tau2_IsoMVArun2v1DBoldDMwLTraw     = TauCandidates[0].tauID("byIsolationMVArun2v1DBoldDMwLTraw");
		tau2_IsoMVArun2v1PWdR03oldDMwLTraw = TauCandidates[0].tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw");
		tau2_IsoMVArun2v1PWnewDMwLTraw     = TauCandidates[0].tauID("byIsolationMVArun2v1PWnewDMwLTraw");
		tau2_IsoMVArun2v1PWoldDMwLTraw     = TauCandidates[0].tauID("byIsolationMVArun2v1PWoldDMwLTraw");
		//
		tau2_looseMuonRejection      = TauCandidates[0].tauID("againstMuonLoose3");
		tau2_tightMuonRejection      = TauCandidates[0].tauID("againstMuonTight3");
		tau2_looseElectronRejection  = TauCandidates[0].tauID("againstElectronLooseMVA6");
		tau2_mediumElectronRejection = TauCandidates[0].tauID("againstElectronMediumMVA6");
		tau2_tightElectronRejection  = TauCandidates[0].tauID("againstElectronTightMVA6");
		tau2_VtightElectronRejection = TauCandidates[0].tauID("againstElectronVTightMVA6");
		LeptonFound = true;
	}

	//return (nLeptonCandidates >= NrequiredLeptons);
	return LeptonFound;
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