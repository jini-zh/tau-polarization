// -*- C++ -*-
//
// Package:    Tau/TreeMaker
// Class:      TreeMaker
//
/**\class TreeMaker TreeMaker.cc Tau/TreeMaker/plugins/TreeMaker.cc

 Description: [one line class summary]

 Implementation:
		 [Notes on implementation]
*/
//
// Original Author:  Dmitry Kondratyev
//         Created:  Wed, 19 Sep 2016
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"

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
#include "DataFormats/JetReco/interface/JetTracksAssociation.h"

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
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedRefCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedRefCandidateFwd.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"

#include "Tau/TauAnalyzer/plugins/MySimpleParticle.h"

static inline double sqr(double x) {
	return x * x;
}

static double dphi(double phi1, double phi2) {
	double result = TMath::Abs(phi1 - phi2);
	if (result < TMath::Pi()) return result;
	return 2 * TMath::Pi() - result;
}

//
// class declaration
//

class TreeMaker : public edm::EDAnalyzer {
public:
	explicit TreeMaker(const edm::ParameterSet&);
	~TreeMaker();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;
	
	bool TriggerOK    (const edm::Event&);
	bool AddTau       (const edm::Event&);
	bool FindGenTau   (const edm::Event&);
	bool CheckMuon    (const edm::Event&);
	bool CheckElectron(const edm::Event&);
	bool AddMET       (const edm::Event&);
	bool AddVertex    (const edm::Event&);
	void CountTracks  (const edm::Event&);
	bool JetPtSum     (const edm::Event&);
	void AddWT        (const edm::Event&);
	void AddWTData    (const edm::Event&);
	bool DecaychannelMatch (std::vector<MySimpleParticle> &particles, int p1, int p2=0, int p3=0, int p4=0, int p5=0, int p6=0);
	void GenTauDM          (const edm::Event&);

	virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	
	// ----------member data ---------------------------


  	HLTConfigProvider             hltConfig_;
  	std::vector<std::string>      trigNames, HLTNames;

	edm::InputTag								triggerEvent_, theTriggerResultsLabel;
	edm::EDGetTokenT<trigger::TriggerEvent>		tok_trigEvt;
	edm::EDGetTokenT<edm::TriggerResults>		tok_trigRes;

	edm::EDGetTokenT<reco::PFTauCollection> TauCollectionToken_;
	edm::EDGetTokenT<reco::MuonCollection> MuonCollectionToken_;
	edm::EDGetTokenT<reco::GsfElectronCollection> ElectronCollectionToken_;
	edm::EDGetTokenT<reco::PFJetCollection> JetCollectionToken_;
	edm::EDGetTokenT<reco::PFJetCollection> CHSJetCollectionToken_;
	edm::EDGetTokenT<reco::PFMETCollection> MetCollectionToken_;
	edm::EDGetTokenT<reco::VertexCollection> PVToken_;
	
	edm::EDGetTokenT<reco::TrackCollection> TrackToken_;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_absIso;
	edm::EDGetTokenT<reco::GenParticleCollection> GenParticleToken_;
	edm::EDGetTokenT<reco::JetTracksAssociationCollection> Token_JetTrack;
	edm::EDGetTokenT<reco::RecoChargedRefCandidateCollection> TrackRefsJets_;
	edm::EDGetTokenT<reco::JetTagCollection> tok_pfDeepCSVJetTags_b;
  	edm::EDGetTokenT<reco::JetTagCollection> tok_pfDeepCSVJetTags_bb;

	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_looseCombinedIso;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_mediumCombinedIso;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_tightCombinedIso;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_looseMvaIso;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_mediumMvaIso;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_tightMvaIso;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_looseMuonRejection;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_tightMuonRejection;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_looseElectronRejection;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_tightElectronRejection;

	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_pfTausDiscriminationByDecayModeFinding;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_hpsPFTauDiscriminationByDecayModeFinding;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_hpsPFTauDiscriminationByDecayModeFindingNewDMs;
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_hpsPFTauDiscriminationByDecayModeFindingOldDMs;

	edm::EDGetTokenT<double> TauSpinnerWTToken_;
	edm::EDGetTokenT<double> TauSpinnerWTFlipToken_;
	edm::EDGetTokenT<double> TauSpinnerWThminusToken_;
	edm::EDGetTokenT<double> TauSpinnerWThplusToken_;
	edm::EDGetTokenT<bool>   TauSpinnerWTisValidToken_;
	edm::EDGetTokenT<int>    TauSpinnerMotherToken_;
	edm::EDGetTokenT<bool>   Photon1Token_;
	edm::EDGetTokenT<bool>   Photon2Token_;
	
	
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
	
	double tau_absIso;
	
	double tau_looseCombinedIso;
	double tau_mediumCombinedIso;
	double tau_tightCombinedIso;
	double tau_looseMvaIso;
	double tau_mediumMvaIso;
	double tau_tightMvaIso;
	double tau_looseMuonRejection;
	double tau_tightMuonRejection;
	double tau_looseElectronRejection;
	double tau_tightElectronRejection;
	double tau_DecayModeFindingNewDMs;
	

	double tau_found;
	double gentau_found;
	double dR;
	double genTauFromW;
	int genTauMother;
	double W_pt;

	double piChar_pt;
	double piChar_eta;
	double piChar_phi;
	double piChar_q;
	double piChar_m;

	double piZero_pt;
	double piZero_eta;
	double piZero_phi;
	double piZero_m;
	
	double pipiMass;
	double upsilon;

	double jetPtSum15, jetPtSum20, jetPtSum15PV, jetPtSum20PV;
	int nJets20, nJets20PV;
	double CHSjetPtSum15, CHSjetPtSum20, CHSjetPtSum15PV, CHSjetPtSum20PV;
	int nCHSJets20, nCHSJets20PV;
	int nLooseBtagedCHSJets, nMediumBtagedCHSJets, nTightBtagedCHSJets;
	int nLooseBtagedCHSJetsPV, nMediumBtagedCHSJetsPV, nTightBtagedCHSJetsPV;
	int nLooseBtagedJets, nMediumBtagedJets, nTightBtagedJets;
	int nLooseBtagedJetsPV, nMediumBtagedJetsPV, nTightBtagedJetsPV;
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
	double BJet_pt;
	double BJet_eta;
	double BJet_phi;
	double BJet_m;
	double BJet_btag;

	int nPi0;
	double met;
	double met_phi;
	double met_eta;
	double met_significance;
	double met_mEtSig;
	
	double m_t;
	
	double dPhi;
	
	math::XYZPoint pv_position;

	int GenHadronDecayChannel;

	double WT;
	double WTFlip;
	double WThminus;
	double WThplus;
	bool WTisValid;
	int  TauSpinnerMother;
	bool PhotonFlag1, PhotonFlag2;

	// Generated tau parameters
	double gentau_pt, gentau_px, gentau_py, gentau_pz, gentau_eta, gentau_phi;
	double genPiChar_pt, genPiChar_px, genPiChar_py, genPiChar_pz, genPiChar_eta, genPiChar_phi;
	double genPi0_pt, genPi0_px, genPi0_py, genPi0_pz, genPi0_eta, genPi0_phi;
	double nuW_pt, nuW_px, nuW_py, nuW_pz, nuW_eta, nuW_phi;
	double nutau_pt, nutau_px, nutau_py, nutau_pz, nutau_eta, nutau_phi;
	double gentau_dm, gentau_nPi0;
	double nunu_pt;

	//////////////////////////////////////////////////////
	bool isMC;
	double tauPtMin;
	double piPtMin;
	double tauEtaMax;
	double tauDzMax;
	double null;
	bool monitoring;
	bool monitoringGen;
	double METcut;

	int iT;
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
TreeMaker::TreeMaker(const edm::ParameterSet& iConfig) {
	//now do what ever initialization is needed
	iT =0;

	isMC						= iConfig.getParameter<bool>("isMC");
	monitoring					= iConfig.getParameter<bool>("monitoring");
	monitoringGen					= iConfig.getParameter<bool>("monitoringGen");
	tauPtMin 					= iConfig.getParameter<double>("tauPtMin");
	piPtMin 					= iConfig.getParameter<double>("piPtMin");
	tauEtaMax 					= iConfig.getParameter<double>("tauEtaMax");
	tauDzMax 					= iConfig.getParameter<double>("tauDzMax");
	null 						= iConfig.getParameter<double>("null");
	METcut                                          = iConfig.getParameter<double>("METcut");
	
	trigNames 					= iConfig.getParameter<std::vector<std::string>>("Triggers");
	triggerEvent_				= edm::InputTag("hltTriggerSummaryAOD","","HLT");
	theTriggerResultsLabel		= edm::InputTag("TriggerResults","","HLT");

	std::string tauCollection         = iConfig.getParameter<std::string>("tauCollection");
	std::string muonCollection        = iConfig.getParameter<std::string>("muonCollection");
	std::string electronCollection    = iConfig.getParameter<std::string>("electronCollection");
	std::string jetCollection         = iConfig.getParameter<std::string>("jetCollection");
	std::string CHSjetCollection         = iConfig.getParameter<std::string>("CHSjetCollection");
	std::string metCollection         = iConfig.getParameter<std::string>("metCollection");
	std::string vertexCollection      = iConfig.getParameter<std::string>("vertexCollection");
	std::string genParticleCollection = iConfig.getParameter<std::string>("genParticleCollection");
	std::string trackCollection       = iConfig.getParameter<std::string>("trackCollection");
	
	std::string absIsoDiscriminator					 = iConfig.getParameter<std::string>("absIsoDiscriminator");
	std::string looseCombinedIsoDiscriminator		 = iConfig.getParameter<std::string>("looseCombinedIsoDiscriminator");
	std::string mediumCombinedIsoDiscriminator		 = iConfig.getParameter<std::string>("mediumCombinedIsoDiscriminator");
	std::string tightCombinedIsoDiscriminator		 = iConfig.getParameter<std::string>("tightCombinedIsoDiscriminator");
	std::string looseMvaIsoDiscriminator			 = iConfig.getParameter<std::string>("looseMvaIsoDiscriminator");
	std::string mediumMvaIsoDiscriminator			 = iConfig.getParameter<std::string>("mediumMvaIsoDiscriminator");
	std::string tightMvaIsoDiscriminator			 = iConfig.getParameter<std::string>("tightMvaIsoDiscriminator");
	std::string looseMuonRejectionDiscriminator		 = iConfig.getParameter<std::string>("looseMuonRejectionDiscriminator");
	std::string tightMuonRejectionDiscriminator		 = iConfig.getParameter<std::string>("tightMuonRejectionDiscriminator");
	std::string looseElectronRejectionDiscriminator	 = iConfig.getParameter<std::string>("looseElectronRejectionDiscriminator");
	std::string tightElectronRejectionDiscriminator	 = iConfig.getParameter<std::string>("tightElectronRejectionDiscriminator");
	std::string pfTauDMFindingDiscriminator       	 = iConfig.getParameter<std::string>("pfTauDMFindingDiscriminator");
	std::string hpspfTauDMFindingDiscriminator	     = iConfig.getParameter<std::string>("hpspfTauDMFindingDiscriminator");
	std::string hpspfTauDMFindingDiscriminatorNewDMs = iConfig.getParameter<std::string>("hpspfTauDMFindingDiscriminatorNewDMs");
	std::string hpspfTauDMFindingDiscriminatorOldDMs = iConfig.getParameter<std::string>("hpspfTauDMFindingDiscriminatorOldDMs");
	
	
	tok_trigEvt					= consumes<trigger::TriggerEvent>(triggerEvent_);
	tok_trigRes					= consumes<edm::TriggerResults>(theTriggerResultsLabel);
	TauCollectionToken_ 		= consumes<reco::PFTauCollection>(edm::InputTag(tauCollection));
	MuonCollectionToken_ 		= consumes<reco::MuonCollection>(edm::InputTag(muonCollection));
	ElectronCollectionToken_	= consumes<reco::GsfElectronCollection>(edm::InputTag(electronCollection));
	JetCollectionToken_ 		= consumes<reco::PFJetCollection>(edm::InputTag(jetCollection));
	CHSJetCollectionToken_          = consumes<reco::PFJetCollection>(edm::InputTag(CHSjetCollection));
	MetCollectionToken_ 		= consumes<reco::PFMETCollection>(edm::InputTag(metCollection));
	PVToken_ 		    = consumes<reco::VertexCollection>(edm::InputTag(vertexCollection));
	GenParticleToken_ 	    = consumes<reco::GenParticleCollection>(edm::InputTag(genParticleCollection));
	TrackToken_		    = consumes<reco::TrackCollection>(edm::InputTag(trackCollection));
	Token_JetTrack              = consumes<reco::JetTracksAssociationCollection>(iConfig.getParameter<edm::InputTag>("ak4JetTracksAssociatorAtVertexPF"));
	TrackRefsJets_              = consumes<reco::RecoChargedRefCandidateCollection>(iConfig.getParameter<edm::InputTag>("trackRefsForJets"));
	// pfDeepCSVJetTags:probb + pfDeepCSVJetTags:probbb working points
	// loose	0.1522
 	// medium	0.4941	 	 
 	// tight	0.8001
	tok_pfDeepCSVJetTags_b      = consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("pfDeepCSVJetTags_b"));
    tok_pfDeepCSVJetTags_bb     = consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("pfDeepCSVJetTags_bb"));

	Token_absIso 				    = consumes<reco::PFTauDiscriminator>(edm::InputTag(absIsoDiscriminator));
	Token_looseCombinedIso 		  	= consumes<reco::PFTauDiscriminator>(edm::InputTag(looseCombinedIsoDiscriminator));
	Token_mediumCombinedIso 	  	= consumes<reco::PFTauDiscriminator>(edm::InputTag(mediumCombinedIsoDiscriminator));
	Token_tightCombinedIso 		 	= consumes<reco::PFTauDiscriminator>(edm::InputTag(tightCombinedIsoDiscriminator));
	Token_looseMvaIso 			    = consumes<reco::PFTauDiscriminator>(edm::InputTag(looseMvaIsoDiscriminator));
	Token_mediumMvaIso 			    = consumes<reco::PFTauDiscriminator>(edm::InputTag(mediumMvaIsoDiscriminator));
	Token_tightMvaIso 			    = consumes<reco::PFTauDiscriminator>(edm::InputTag(tightMvaIsoDiscriminator));
	Token_looseMuonRejection  		= consumes<reco::PFTauDiscriminator>(edm::InputTag(looseMuonRejectionDiscriminator));
	Token_tightMuonRejection 	  	= consumes<reco::PFTauDiscriminator>(edm::InputTag(tightMuonRejectionDiscriminator));
	Token_looseElectronRejection	= consumes<reco::PFTauDiscriminator>(edm::InputTag(looseElectronRejectionDiscriminator));
	Token_tightElectronRejection	= consumes<reco::PFTauDiscriminator>(edm::InputTag(tightElectronRejectionDiscriminator));
	Token_pfTausDiscriminationByDecayModeFinding         = consumes<reco::PFTauDiscriminator>(edm::InputTag(pfTauDMFindingDiscriminator));
	Token_hpsPFTauDiscriminationByDecayModeFinding       = consumes<reco::PFTauDiscriminator>(edm::InputTag(hpspfTauDMFindingDiscriminator));
	Token_hpsPFTauDiscriminationByDecayModeFindingNewDMs = consumes<reco::PFTauDiscriminator>(edm::InputTag(hpspfTauDMFindingDiscriminatorNewDMs));
	Token_hpsPFTauDiscriminationByDecayModeFindingOldDMs = consumes<reco::PFTauDiscriminator>(edm::InputTag(hpspfTauDMFindingDiscriminatorOldDMs));

	if (isMC) {
		TauSpinnerWTToken_        = consumes<double>(iConfig.getParameter<edm::InputTag>("WTCollection"));
		TauSpinnerWTFlipToken_    = consumes<double>(iConfig.getParameter<edm::InputTag>("WTFlipCollection"));
		TauSpinnerWThminusToken_  = consumes<double>(iConfig.getParameter<edm::InputTag>("WThminusCollection"));
		TauSpinnerWThplusToken_   = consumes<double>(iConfig.getParameter<edm::InputTag>("WThplusCollection"));
		TauSpinnerWTisValidToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("WTisValidCollection"));
		TauSpinnerMotherToken_    = consumes<int>(iConfig.getParameter<edm::InputTag>("MotherCollection"));
		Photon1Token_             = consumes<bool>(iConfig.getParameter<edm::InputTag>("Photon1Flag"));
		Photon2Token_             = consumes<bool>(iConfig.getParameter<edm::InputTag>("Photon2Flag"));
	}

}


TreeMaker::~TreeMaker() {
	 // do anything here that needs to be done at desctruction time
	 // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//GenTauDM

void TreeMaker::analyze(const edm::Event& event, const edm::EventSetup&) {
	t_Run   = event.id().run();
	t_Event = event.id().event();

	if (!TriggerOK(event)) {
		if (monitoring) std::cout << "Trigger" << std::endl;
		return;
	}
	if (!AddVertex(event)) {
		if (monitoring) std::cout << "Vertex" << std::endl;
		return;
	}
	if (CheckMuon(event)) {
		if (monitoring) std::cout << "Muon" << std::endl;
		return;
	}
	if (CheckElectron(event)) {
		if (monitoring) std::cout << "Electron" << std::endl;
		return;
	}
	if (!AddMET(event)) {
		if (monitoring) std::cout << "MET" << std::endl;
		return;
	}
	if (!AddTau(event)) {
		if (monitoring) std::cout << "Tau" << std::endl;
		return;
	}
	FindGenTau(event);
	GenTauDM(event);
	CountTracks(event);
	if (!JetPtSum(event)) {
		if (monitoring) std::cout << "Jet" << std::endl;
		return;
	}
	if (isMC) {
		AddWT(event);
	} else {
		AddWTData(event);
	}

	/*
	if (monitoring) {
		std::cout << "Gentau Decay channel = " << GenHadronDecayChannel << std::endl;
		std::cout << "WT           = " << WT << std::endl;
		std::cout << "gentau_found = " << gentau_found << std::endl;
		std::cout << "tau_found = " << tau_found << std::endl;
	}
	*/

	TLorentzVector pi0, pi1;
	pi0.SetPtEtaPhiM(piZero_pt, piZero_eta, piZero_phi, piZero_m);
	pi1.SetPtEtaPhiM(piChar_pt, piChar_eta, piChar_phi, piChar_m);

	pipiMass = (pi0 + pi1).M();
//	upsilon  = (pi1.E() - pi0.E()) / tau_pt;
	upsilon  = 2 * piChar_pt / tau_pt - 1;
	dPhi	 = dphi(tau_phi, met_phi);
	m_t      = sqrt(2 * tau_pt * met * (1 - cos(dPhi)));

	if (monitoring) {
		std::cout << "!!! Event was writen to the tree" << std::endl;
		std::cout << "tau_pt       = " << tau_pt << std::endl;
		std::cout << "tau_dm       = " << tau_dm << std::endl;
		std::cout << "gentau_found = " << gentau_found << std::endl;
		std::cout << "gentau_pt    = " << gentau_pt << std::endl;
		std::cout << "gentau_dm    = " << gentau_dm << std::endl;
	}

	tree->Fill();
};


// ------------ method called once each job just before starting event loop  ------------
void TreeMaker::beginJob() {
	// declaring the tree and its branches.
	edm::Service<TFileService> FS;
	tree = FS->make<TTree>("tree", "tree", 1);
	
	tree->Branch("t_Run",  &t_Run,  "t_Run/I");
	tree->Branch("t_Event",&t_Event,"t_Event/I");
	
	tree->Branch("tau_pt",&tau_pt,"tau_pt/D");
	tree->Branch("tau_eta",&tau_eta,"tau_eta/D");
	tree->Branch("tau_phi",&tau_phi,"tau_phi/D");
	tree->Branch("tau_q",&tau_q,"tau_q/D");
	tree->Branch("tau_m",&tau_m,"tau_m/D");
	tree->Branch("tau_dm",&tau_dm,"tau_dm/D");
	tree->Branch("tau_dz",&tau_dz,"tau_dz/D");
	tree->Branch("tau_absIso",&tau_absIso,"tau_absIso/D");
	tree->Branch("tau_looseCombinedIso",&tau_looseCombinedIso,"tau_looseCombinedIso/D");
	tree->Branch("tau_mediumCombinedIso",&tau_mediumCombinedIso,"tau_mediumCombinedIso/D");
	tree->Branch("tau_tightCombinedIso",&tau_tightCombinedIso,"tau_tightCombinedIso/D");
	tree->Branch("tau_looseMvaIso",&tau_looseMvaIso,"tau_looseMvaIso/D");
	tree->Branch("tau_mediumMvaIso",&tau_mediumMvaIso,"tau_mediumMvaIso/D");
	tree->Branch("tau_tightMvaIso",&tau_tightMvaIso,"tau_tightMvaIso/D");
	tree->Branch("tau_looseMuonRejection",&tau_looseMuonRejection,"tau_looseMuonRejection/D");
	tree->Branch("tau_tightMuonRejection",&tau_tightMuonRejection,"tau_tightMuonRejection/D");
	tree->Branch("tau_looseElectronRejection",&tau_looseElectronRejection,"tau_looseElectronRejection/D");
	tree->Branch("tau_tightElectronRejection",&tau_tightElectronRejection,"tau_tightElectronRejection/D");
	tree->Branch("tau_DecayModeFindingNewDMs",&tau_DecayModeFindingNewDMs,"tau_DecayModeFindingNewDMs/D");
	
	tree->Branch("piChar_pt", &piChar_pt, "piChar_pt/D");
	tree->Branch("piChar_eta", &piChar_eta, "piChar_eta/D");
	tree->Branch("piChar_phi", &piChar_phi, "piChar_phi/D");
	tree->Branch("piChar_q", &piChar_q, "piChar_q/D");
	tree->Branch("piChar_m", &piChar_m, "piChar_m/D");

	tree->Branch("piZero_pt", &piZero_pt, "piZero_pt/D");
	tree->Branch("piZero_eta", &piZero_eta, "piZero_eta/D");
	tree->Branch("piZero_phi", &piZero_phi, "piZero_phi/D");
	tree->Branch("piZero_m", &piZero_m, "piZero_m/D");

	tree->Branch("pipiMass", &pipiMass, "pipiMass/D");
	tree->Branch("upsilon", &upsilon, "upsilon/D");

	tree->Branch("tau_found",&tau_found, "tau_found/D");
	tree->Branch("gentau_found",&gentau_found, "gentau_found/D");
	tree->Branch("dR", &dR, "dR/D");
	tree->Branch("genTauFromW", &genTauFromW, "genTauFromW/D");
	tree->Branch("W_pt", &W_pt, "W_pt/D");

	tree->Branch("genTauMother", &genTauMother, "genTauMother/I");

	// Jets
	tree->Branch("jetPtSum15", &jetPtSum15, "jetPtSum15/D");
	tree->Branch("jetPtSum20", &jetPtSum20, "jetPtSum20/D");
	tree->Branch("jetPtSum15PV", &jetPtSum15PV, "jetPtSum15PV/D");
	tree->Branch("jetPtSum20PV", &jetPtSum20PV, "jetPtSum20PV/D");
	tree->Branch("nJets20", &nJets20, "nJets20/I");
	tree->Branch("nJets20PV", &nJets20PV, "nJets20PV/I");
	tree->Branch("nLooseBtagedJets", &nLooseBtagedJets, "nLooseBtagedJets/I");
	tree->Branch("nMediumBtagedJets", &nMediumBtagedJets, "nMediumBtagedJets/I");
	tree->Branch("nTightBtagedJets", &nTightBtagedJets, "nTightBtagedJets/I");
	tree->Branch("nLooseBtagedJetsPV", &nLooseBtagedJetsPV, "nLooseBtagedJetsPV/I");
	tree->Branch("nMediumBtagedJetsPV", &nMediumBtagedJetsPV, "nMediumBtagedJetsPV/I");
	tree->Branch("nTightBtagedJetsPV", &nTightBtagedJetsPV, "nTightBtagedJetsPV/I");

	// CHSJets
	tree->Branch("CHSjetPtSum15", &CHSjetPtSum15, "CHSjetPtSum15/D");
	tree->Branch("CHSjetPtSum20", &CHSjetPtSum20, "CHSjetPtSum20/D");
	tree->Branch("CHSjetPtSum15PV", &CHSjetPtSum15PV, "CHSjetPtSum15PV/D");
	tree->Branch("CHSjetPtSum20PV", &CHSjetPtSum20PV, "CHSjetPtSum20PV/D");
	tree->Branch("nCHSJets20", &nCHSJets20, "nCHSJets20/I");
	tree->Branch("nCHSJets20PV", &nCHSJets20PV, "nCHSJets20PV/I");
	tree->Branch("nLooseBtagedCHSJets", &nLooseBtagedCHSJets, "nLooseBtagedCHSJets/I");
	tree->Branch("nMediumBtagedCHSJets", &nMediumBtagedCHSJets, "nMediumBtagedCHSJets/I");
	tree->Branch("nTightBtagedCHSJets", &nTightBtagedCHSJets, "nTightBtagedCHSJets/I");
	tree->Branch("nLooseBtagedCHSJetsPV", &nLooseBtagedCHSJetsPV, "nLooseBtagedCHSJetsPV/I");
	tree->Branch("nMediumBtagedCHSJetsPV", &nMediumBtagedCHSJetsPV, "nMediumBtagedCHSJetsPV/I");
	tree->Branch("nTightBtagedCHSJetsPV", &nTightBtagedCHSJetsPV, "nTightBtagedCHSJetsPV/I");

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
	tree->Branch("BJet_pt", &BJet_pt, "BJet_pt/D");
	tree->Branch("BJet_eta", &BJet_eta, "BJet_eta/D");
	tree->Branch("BJet_phi", &BJet_phi, "BJet_phi/D");
	tree->Branch("BJet_m", &BJet_m, "BJet_m/D");
	tree->Branch("BJet_btag", &BJet_btag, "BJet_btag/D");
	
	tree->Branch("met", &met, "met/D");
	tree->Branch("met_phi", &met_phi, "met_phi/D");
	tree->Branch("met_eta", &met_eta, "met_eta/D");
	tree->Branch("met_significance", &met_significance, "met_significance/D");
	tree->Branch("met_mEtSig", &met_mEtSig, "met_mEtSig/D");
	tree->Branch("m_t", &m_t, "m_t/D");
	tree->Branch("dPhi", &dPhi, "dPhi/D");
	
	tree->Branch("nVtx",&nVtx,"nVtx/I");
	tree->Branch("nTrks",&nTrks,"nTrks/I");
	tree->Branch("nTau",&nTau,"nTau/I");
	tree->Branch("nTauC",&nTauC,"nTauC/I");

	tree->Branch("nPi0",&nPi0,"nPi0/I");

	tree->Branch("WT", &WT, "WT/D");
	tree->Branch("WTFlip", &WTFlip, "WTFlip/D");
	tree->Branch("WThminus", &WThminus, "WThminus/D");
	tree->Branch("WThplus", &WThplus, "WThplus/D");
	tree->Branch("TauSpinnerMother", &TauSpinnerMother, "TauSpinnerMother/I");
	tree->Branch("PhotonFlag1", &PhotonFlag1, "PhotonFlag1/B");
	tree->Branch("PhotonFlag2", &PhotonFlag2, "PhotonFlag2/B");
	tree->Branch("GenHadronDecayChannel", &GenHadronDecayChannel, "GenHadronDecayChannel/I");

	tree->Branch("gentau_pt",&gentau_pt,"gentau_pt/D");
	tree->Branch("genPiChar_pt",&genPiChar_pt,"genPiChar_pt/D");
	tree->Branch("genPi0_pt",&genPi0_pt,"genPi0_pt/D");
	tree->Branch("nutau_pt",&nutau_pt,"nutau_pt/D");
	tree->Branch("nuW_pt",&nuW_pt,"nuW_pt/D");
	tree->Branch("nunu_pt",&nunu_pt,"nunu_pt/D");
	tree->Branch("gentau_dm",&gentau_dm,"gentau_dm/D");
	tree->Branch("gentau_nPi0",&gentau_nPi0,"gentau_nPi0/D");
	//tree->Branch("WTisValid", &WTisValid, "WTisValid/D")
	// add more branches
	
	allTauPt = FS->make<TH1F>("allTauPt","allTauPt",300,0,300);
}

// ------------ method called once each job just after ending the event loop  ------------
void TreeMaker::endJob() {

}

// ------------ method called when starting to processes a run  ------------

void TreeMaker::beginRun(const edm::Run& Run, const edm::EventSetup& Setup) {
}

// ------------ method called when ending the processing of a run  ------------
/*
void
TreeMaker::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TreeMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

//---------------------------------TRIGGER-------------------------------------------------
bool TreeMaker::TriggerOK(const edm::Event& iEvent){
	bool triggerOK = false;
	if (isMC) {
		triggerOK = true; // ignore HLT for single pion MC
	} else {
		trigger::TriggerEvent triggerEvent;
		edm::Handle<trigger::TriggerEvent> triggerEventHandle;
		iEvent.getByToken(tok_trigEvt, triggerEventHandle);
		if (triggerEventHandle.isValid()) {
			triggerEvent = *(triggerEventHandle.product());
			const trigger::TriggerObjectCollection& TOC(triggerEvent.getObjects());
    /////////////////////////////TriggerResults////////////////////////////////////
			edm::Handle<edm::TriggerResults> triggerResults;
			iEvent.getByToken(tok_trigRes, triggerResults);
			if (triggerResults.isValid()) {
				const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
				const std::vector<std::string> & triggerNames_ = triggerNames.triggerNames();
				for ( unsigned int iHLT=0; iHLT<triggerResults->size(); iHLT++ ) {
					int hlt    = triggerResults->accept(iHLT);
					if ( hlt > 0 ) {
						for ( unsigned int i=0; i<trigNames.size(); ++i ) {
							if ( triggerNames_[iHLT].find(trigNames[i].c_str())!= std::string::npos ) {
								triggerOK = true;
							}
						}
					}
				}
			}
		}
	}
	return triggerOK;
}


//-----------------------------------------------------------------------------------------

bool TreeMaker::AddTau(const edm::Event& event) {
	tau_found = 0;
	nTauC     = 0;

	edm::Handle<reco::PFTauCollection> taus;
	event.getByToken(TauCollectionToken_, taus);
	if (!taus.isValid() || taus->size() == 0) return false;

	edm::Handle<reco::PFTauDiscriminator> absIso;
	edm::Handle<reco::PFTauDiscriminator> looseCombinedIso;
	edm::Handle<reco::PFTauDiscriminator> mediumCombinedIso;
	edm::Handle<reco::PFTauDiscriminator> tightCombinedIso;
	edm::Handle<reco::PFTauDiscriminator> looseMvaIso;
	edm::Handle<reco::PFTauDiscriminator> mediumMvaIso;
	edm::Handle<reco::PFTauDiscriminator> tightMvaIso;
	edm::Handle<reco::PFTauDiscriminator> looseMuonRejection;
	edm::Handle<reco::PFTauDiscriminator> tightMuonRejection;
	edm::Handle<reco::PFTauDiscriminator> looseElectronRejection;
	edm::Handle<reco::PFTauDiscriminator> tightElectronRejection;
	edm::Handle<reco::PFTauDiscriminator> pfTausDiscriminationByDecayModeFinding;
	edm::Handle<reco::PFTauDiscriminator> hpsPFTauDiscriminationByDecayModeFinding;
	edm::Handle<reco::PFTauDiscriminator> hpsPFTauDiscriminationByDecayModeFindingNewDMs;
	edm::Handle<reco::PFTauDiscriminator> hpsPFTauDiscriminationByDecayModeFindingOldDMs;

	event.getByToken(Token_absIso, absIso);
	event.getByToken(Token_looseCombinedIso, looseCombinedIso);
	event.getByToken(Token_mediumCombinedIso, mediumCombinedIso);
	event.getByToken(Token_tightCombinedIso, tightCombinedIso);
	event.getByToken(Token_looseMvaIso, looseMvaIso);
	event.getByToken(Token_mediumMvaIso, mediumMvaIso);
	event.getByToken(Token_tightMvaIso, tightMvaIso);
	event.getByToken(Token_looseMuonRejection, looseMuonRejection);
	event.getByToken(Token_tightMuonRejection, tightMuonRejection);
	event.getByToken(Token_looseElectronRejection, looseElectronRejection);
	event.getByToken(Token_tightElectronRejection, tightElectronRejection);
	event.getByToken(Token_pfTausDiscriminationByDecayModeFinding, pfTausDiscriminationByDecayModeFinding);
	event.getByToken(Token_hpsPFTauDiscriminationByDecayModeFinding, hpsPFTauDiscriminationByDecayModeFinding);
	event.getByToken(Token_hpsPFTauDiscriminationByDecayModeFindingNewDMs, hpsPFTauDiscriminationByDecayModeFindingNewDMs);
	event.getByToken(Token_hpsPFTauDiscriminationByDecayModeFindingOldDMs, hpsPFTauDiscriminationByDecayModeFindingOldDMs);

	// search for the tau candidate with the minimum isolation and the
	// maximum transverse momentum
	size_t index  = taus->size();;
	for (size_t i = 0; i < taus->size(); ++i) {
#define cut(condition) if (!(condition)) continue;
		auto& tau = (*taus)[i];
		cut(tau.pt() > tauPtMin); // tau transverse momentum
		cut(TMath::Abs(tau.eta()) < tauEtaMax); // tau pseudorapidity

		auto& chargedHadrons = tau.signalPFChargedHadrCands();
		cut(chargedHadrons.size() == 1); // single charged hadron (pi+)
		auto& pi0s = tau.signalPiZeroCandidates();
		cut(pi0s.size() > 0); // at least one pi0

		auto& pi0 = pi0s.front();
		auto& pi1 = *chargedHadrons.front();

		cut(pi1.pt() > piPtMin); // pi+ transverse momentum

		++nTauC;
		allTauPt->Fill(tau.pt());

		cut((pv_position - tau.vertex()).R() < tauDzMax); // tau vertex displacement

		if (index < taus->size()) {
			double iso = absIso->value(i);
			double minIso = absIso->value(index);
			cut(
			      iso < minIso
			      || iso == minIso && tau.pt() > (*taus)[index].pt()
			)
		}

		index = i;
#undef cut
	};

	if (index == taus->size()) return false;

	auto& tau = (*taus)[index];
	auto& pi0 = tau.signalPiZeroCandidates().front();
	auto& pi1 = *tau.signalPFChargedHadrCands().front();
	tau_pt                     = tau.pt();
	tau_eta                    = tau.eta();
	tau_phi                    = tau.phi();
	tau_dm                     = tau.decayMode();
	tau_dz                     = (pv_position - tau.vertex()).R();
	tau_q                      = tau.charge();
	tau_m                      = tau.mass();
	tau_absIso                 = absIso->value(index);
	tau_looseCombinedIso       = looseCombinedIso->value(index);
	tau_mediumCombinedIso      = mediumCombinedIso->value(index);
	tau_tightCombinedIso       = tightCombinedIso->value(index);
	tau_looseMvaIso            = looseMvaIso->value(index);
	tau_mediumMvaIso           = mediumMvaIso->value(index);
	tau_tightMvaIso            = tightMvaIso->value(index);
	tau_looseMuonRejection     = looseMuonRejection->value(index);
	tau_tightMuonRejection     = tightMuonRejection->value(index);
	tau_looseElectronRejection = looseElectronRejection->value(index);
	tau_tightElectronRejection = tightElectronRejection->value(index);
	tau_DecayModeFindingNewDMs = hpsPFTauDiscriminationByDecayModeFindingNewDMs->value(index);

	/*
	if (tau_dm > 1) {
	  std::cout << "Tau pt = " << tau_pt << std::endl;
	  std::cout << "tau_dm = " << tau_dm << std::endl; 
	  std::cout << "pfTausDiscriminationByDecayModeFinding         = " << pfTausDiscriminationByDecayModeFinding->value(index) << std::endl;
	  std::cout << "hpsPFTauDiscriminationByDecayModeFinding       = " << hpsPFTauDiscriminationByDecayModeFinding->value(index) << std::endl;
	  std::cout << "hpsPFTauDiscriminationByDecayModeFindingNewDMs = " << hpsPFTauDiscriminationByDecayModeFindingNewDMs->value(index) << std::endl;
	  std::cout << "hpsPFTauDiscriminationByDecayModeFindingOldDMs = " << hpsPFTauDiscriminationByDecayModeFindingOldDMs->value(index) << std::endl;
	}
	*/

	piZero_pt  = pi0.pt();
	piZero_eta = pi0.eta();
	piZero_phi = pi0.phi();
	piZero_m   = pi0.mass();

	piChar_pt  = pi1.pt();
	piChar_eta = pi1.eta();
	piChar_phi = pi1.phi();
	piChar_q   = pi1.charge();
	piChar_m   = pi1.mass();

	nPi0      = tau.signalPiZeroCandidates().size();
	nTau      = taus->size();
	tau_found = 1;

	return true;
};

bool TreeMaker::FindGenTau(const edm::Event& event) {
	gentau_found  = 0;
	genTauFromW   = null;
	genTauMother  = null;
	dR            = null;
	gentau_pt     = null;
	genPiChar_pt  = null;
	genPi0_pt     = null;
	nutau_pt      = null;
	nuW_pt        = null;
	nunu_pt       = null;
	double nuW_m  = 0.;

	edm::Handle<reco::GenParticleCollection> genParticles;
	event.getByToken(GenParticleToken_, genParticles);
	if (!genParticles.isValid()) return false;

	if (monitoringGen) std::cout << "GenParticles is valid" << std::endl;

	const int pdg_tau    = 15;
	const int pdg_pi0    = 111;
	const int pdg_pi1    = 211;
	const int pdg_W      = 24;
	const int pdg_nu_tau = 16;

	const reco::Candidate* tau    = nullptr;
	const reco::Candidate* pi0    = nullptr;
	const reco::Candidate* pi1    = nullptr;
	const reco::Candidate* nu_W   = nullptr;
	const reco::Candidate* nu_tau = nullptr;
	double dRmin = null;
	for (auto& particle: *genParticles) {
		// look for the tau -> pi+ pi0 neutrino decay most oriented towards
		// the reconstructed tau (if present)
#define cut(condition) if (!(condition)) continue;
		cut(abs(particle.pdgId()) == pdg_tau);
		if (monitoringGen) {
			std::cout << "Tau found in GP collection (pt = " << particle.pt() << ", eta = " << particle.eta() << ")" << std::endl;
			std::cout << "Tau has " << particle.numberOfDaughters() << " daughters:" << std::endl;
			for (unsigned j = 0; j < particle.numberOfDaughters(); ++j) {
				const reco::Candidate* daughter = particle.daughter(j);
				int id = abs(daughter->pdgId());
				std::cout << "pdgId(" << j << ") = " << id << std::endl;
			}
		}
		cut(particle.numberOfDaughters() == 3);
		for (int i = 0; i < 3; ++i) {
			const reco::Candidate* daughter = particle.daughter(i);
			int id = abs(daughter->pdgId());
			if (id == pdg_pi0)
				pi0 = daughter;
			else if (id == pdg_pi1)
				pi1 = daughter;
			else if (id == pdg_nu_tau)
				nu_tau = daughter;
		};
		cut(pi0 && pi1);
		if (monitoringGen) std::cout << "PiChar and Pi0 among daughters" << std::endl;

		cut(particle.pt() > tauPtMin);
		if (monitoringGen) std::cout << "pass pt cut" << std::endl;
		cut(TMath::Abs(particle.eta()) < tauEtaMax);
		if (monitoringGen) std::cout << "pass eta cut" << std::endl;
		cut((pv_position - particle.vertex()).R() < tauDzMax);
		if (monitoringGen) std::cout << "pass Dz cut" << std::endl;
		cut(pi1->pt() > piPtMin);
		if (monitoringGen) std::cout << "pass PiChar_pt cut" << std::endl;

		double dR_ = null;
		if (tau_found) {
			dR_ = TMath::Sqrt(
					  sqr(dphi(particle.phi(), tau_phi))
					+ sqr(particle.eta() - tau_eta)
			);
			if (!tau || dR_ < dRmin) {
				tau = &particle;
				dRmin = dR_;
			};
		};
		if (monitoringGen) std::cout << "distance from reco tau = " << dRmin << std::endl;
#undef cut
	};
	if (!tau) {
		if (monitoringGen) std::cout << "gen has not been written to Tree" << std::endl;
		return false;
	}

	gentau_pt  = tau->pt();
	gentau_px  = tau->px();
	gentau_py  = tau->py();
	gentau_pz  = tau->pz();
	gentau_eta = tau->eta();
	gentau_phi = tau->phi();

	dR           = dRmin;
	gentau_found = 1;
	genTauFromW  = 0;
	for (auto p = tau->mother(); p; p = p->mother()) {
		genTauMother = p->pdgId();
		if (abs(p->pdgId()) == pdg_W) {
			genTauFromW = 1;
			W_pt        = p->pt();
			if (monitoringGen) std::cout << "W was found, daughters:" << std::endl;
			for (unsigned l = 0; l < p->numberOfDaughters(); l++) {
				const reco::Candidate* Wdaughter= p->daughter(l);
				if (monitoringGen) std::cout << "W_daughter(" << l << ") = " << Wdaughter->pdgId() << std::endl;
				if (abs(Wdaughter->pdgId()) == pdg_nu_tau) {
					nu_W = Wdaughter;
				} else continue;
			}
			break;
		};

	};

	// Parameters of tau daughters from generator
	gentau_pt  = tau->pt();
	gentau_px  = tau->px();
	gentau_py  = tau->py();
	gentau_pz  = tau->pz();
	gentau_eta = tau->eta();
	gentau_phi = tau->phi();

	genPiChar_pt  = pi1->pt();
	genPiChar_px  = pi1->px();
	genPiChar_py  = pi1->py();
	genPiChar_pz  = pi1->pz();
	genPiChar_eta = pi1->eta();
	genPiChar_phi = pi1->phi();
	double genPiChar_m = pi1->mass();

	genPi0_pt  = pi0->pt();
	genPi0_px  = pi0->px();
	genPi0_py  = pi0->py();
	genPi0_pz  = pi0->pz();
	genPi0_eta = pi0->eta();
	genPi0_phi = pi0->phi();
	double genPi0_m = pi0->mass();

	nutau_pt  = nu_tau->pt();
	nutau_px  = nu_tau->px();
	nutau_py  = nu_tau->py();
	nutau_pz  = nu_tau->pz();
	nutau_eta = nu_tau->eta();
	nutau_phi = nu_tau->phi();
	double nutau_m = nu_tau->mass();
	if (genTauFromW > 0) {
		nuW_pt  = nu_W->pt();
		nuW_px  = nu_W->px();
		nuW_py  = nu_W->py();
		nuW_pz  = nu_W->pz();
		nuW_eta = nu_W->eta();
		nuW_phi = nu_W->phi();
		double nuW_m = nu_W->mass();
	}

	if (genTauFromW > 0) {
		TLorentzVector genpi0, genpi1, gennutau, gennuW;
		genpi1.SetPtEtaPhiM(genPiChar_pt, genPiChar_eta, genPiChar_phi, genPiChar_m);
		genpi0.SetPtEtaPhiM(genPi0_pt, genPi0_eta, genPi0_phi, genPi0_m);
		gennutau.SetPtEtaPhiM(nutau_pt, nutau_eta, nutau_phi, nutau_m);
		gennuW.SetPtEtaPhiM(nuW_pt, nuW_eta, nuW_phi, nuW_m);

		nunu_pt = (gennutau + gennuW).Pt();
	}
	
	if (monitoringGen) {
		std::cout << "Final result Gen:" << std::endl;
		std::cout << "nutau_pt          = " << nutau_pt << std::endl;
		std::cout << "nuW_pt            = " << nuW_pt << std::endl;
		std::cout << "nunu_pt           = " << nunu_pt << std::endl;
	}
	
	return true;
};

bool TreeMaker::CheckMuon(const edm::Event& event) {
	edm::Handle<reco::MuonCollection> muons;
	event.getByToken(MuonCollectionToken_, muons);
	if (!muons.isValid()) return false;
	for (auto& muon: *muons) if (muon.pt() > 15) return true;
	return false;
}

bool TreeMaker::CheckElectron(const edm::Event& event) {
	edm::Handle<reco::GsfElectronCollection> electrons;
	event.getByToken(ElectronCollectionToken_, electrons);
	if (!electrons.isValid()) return false;
	for (auto& electron: *electrons) if (electron.pt() > 15) return true;
	return false;
}

bool TreeMaker::AddMET(const edm::Event& event) {
	edm::Handle<reco::PFMETCollection> mets;
	event.getByToken(MetCollectionToken_, mets);
	if (!mets.isValid() || !mets->size()) return false;

	auto& MET = mets->front();
	met = MET.pt();
	met_phi = MET.phi();
	met_eta = MET.eta();
	met_significance = MET.significance();
	met_mEtSig       = MET.mEtSig();
	if (met < METcut) return false;
	return true;
}

bool TreeMaker::AddVertex(const edm::Event& event) {
	edm::Handle<reco::VertexCollection> vertices;
	event.getByToken(PVToken_, vertices);
	if (!vertices.isValid()) return false;

	nVtx = vertices->size();
	if (nVtx == 0) return false;
	pv_position = vertices->front().position();
	return true;
};

// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookAssociationVector
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMRef

// Products to test:
// trackRefsForJets
// ak4TrackJets

bool TreeMaker::JetPtSum(const edm::Event& event) {
	// Jets
	jetPtSum15 = 0;
	jetPtSum20 = 0;
	nJets20    = 0;
	jetPtSum15PV = 0;
	jetPtSum20PV = 0;
	nJets20PV    = 0;
	nLooseBtagedJets    = 0;
	nMediumBtagedJets   = 0;
	nTightBtagedJets    = 0;
	nLooseBtagedJetsPV  = 0;
	nMediumBtagedJetsPV = 0;
	nTightBtagedJetsPV  = 0;
	// CHS Jets
	CHSjetPtSum15 = 0;
	CHSjetPtSum20 = 0;
	nCHSJets20    = 0;
	CHSjetPtSum15PV = 0;
	CHSjetPtSum20PV = 0;
	nCHSJets20PV    = 0;
	nLooseBtagedCHSJets    = 0;
	nMediumBtagedCHSJets   = 0;
	nTightBtagedCHSJets    = 0;
	nLooseBtagedCHSJetsPV  = 0;
	nMediumBtagedCHSJetsPV = 0;
	nTightBtagedCHSJetsPV  = 0;

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

	BJet_pt   = null;
	BJet_eta  = null;
	BJet_phi  = null;
	BJet_m    = null;
	BJet_btag = null;

	edm::Handle<reco::PFJetCollection> jets;
	event.getByToken(JetCollectionToken_, jets);
	if (!jets.isValid()) return false;
	edm::Handle<reco::PFJetCollection> CHSjets;
	event.getByToken(CHSJetCollectionToken_, CHSjets);
	if (!CHSjets.isValid()) return false;
	edm::Handle<reco::TrackCollection> tracks;
	event.getByToken(TrackToken_, tracks);
	if (!tracks.isValid()) return false;
	edm::Handle<reco::VertexCollection> vertices;
	event.getByToken(PVToken_, vertices);
	if (!vertices.isValid() || vertices->size() == 0) return false;

	edm::Handle<reco::JetTagCollection> pfDeepCSVJetTags_b;
    	event.getByToken(tok_pfDeepCSVJetTags_b, pfDeepCSVJetTags_b);
  	edm::Handle<reco::JetTagCollection> pfDeepCSVJetTags_bb;
  	event.getByToken(tok_pfDeepCSVJetTags_bb, pfDeepCSVJetTags_bb);

  	double WPBTag_loose   = 0.1522;
 	double WPBTag_medium  = 0.4941; 	 
 	double WPBTag_tight   = 0.8001;

  	const reco::JetTagCollection & bTag  = *(pfDeepCSVJetTags_b.product());
  	const reco::JetTagCollection & bbTag = *(pfDeepCSVJetTags_bb.product());

	// Primary vertex or 0 element of vertices collection
	if (monitoring) {
		std::cout << "pv position = (" << pv_position.X() << ", " << pv_position.Y() << ", " << pv_position.Z() << ")" << std::endl;
		std::cout << "Jets:"<< std::endl;
	}

	int i = 0;
	for (auto& jet: *jets) {
		if (TMath::Abs(jet.eta()) > 3) continue;
		if (jet.pt() > 15) {
			i++;
			if (monitoring) std::cout << std::endl << "Jet[" << i << "] Pt = " << jet.pt() << std::endl;
			reco::TrackRefVector JetTracksRef = jet.getTrackRefs();
			std::map <int, int> VertexTracksMap;
			std::map <int, std::pair<int, double>> VertexTracksPtMap;
			double JetDeepCSV_prob_b = 0;
			double JetDeepCSV_prob_bb = 0;
			int NtracksFromPV = 0;
			double SumtracksPtFromPV = 0;
			bool JetFromPV = true;
			int ntrk1 = 0;
			std::vector<int> TrackCounter;
			std::vector<double> TracksPtSum;
			//int TrackCounter[vertices->size()];
			//double TracksPtSum[vertices->size()];
			for (unsigned count = 0; count < vertices->size(); count++) {
				//TrackCounter[count] = 0;
				//TracksPtSum[count] = 0;
				TrackCounter.push_back(0);
				TracksPtSum.push_back(0);
			}
			// Loop over tracks from Jet
			for (reco::TrackRefVector::const_iterator i_trk = JetTracksRef.begin(); i_trk != JetTracksRef.end(); i_trk++, ntrk1++) {
				// Loop over vrtices
				for (unsigned nvtx = 0; nvtx < vertices->size(); nvtx++) {
					// number of tracks from this vertex in jet
					VertexTracksMap.insert(std::pair<int, int> (nvtx, TrackCounter.at(nvtx)));
					//std::pair <int, double> Pair1 (TrackCounter[nvtx], TracksPtSum[nvtx]);
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
			for (auto Map_iter1 = VertexTracksPtMap.begin(); Map_iter1 != VertexTracksPtMap.end(); ++Map_iter1) {
				if ((*Map_iter1).second.first < 1) continue;
				if ((*Map_iter1).first == 0) SumtracksPtFromPV = (*Map_iter1).second.second;
				if ((*Map_iter1).first > 0 && (*Map_iter1).second.second > SumtracksPtFromPV) JetFromPV = false;
				if (monitoring) std::cout << "Vertex number   = " << (*Map_iter1).first << "  entries = " << (*Map_iter1).second.first << "  SumPt = " << (*Map_iter1).second.second << std::endl;
			}

			// Loop over bTags
    		for (unsigned i = 0; i != bTag.size(); i++) {
    			if (abs(bTag[i].first->eta() - jet.eta()) < 0.025 && abs(bTag[i].first->phi() - jet.phi()) < 0.025 ){
          			JetDeepCSV_prob_b = bTag[i].second;
          			if (monitoring) {
            			std::cout << "##### Btag associated with jet number " << std::endl;
            			std::cout << "BTagCollection   eta = " << bTag[i].first->eta() << "   phi = " << bTag[i].first->phi() << std::endl;
            			std::cout << "PFCollection     eta = " << jet.eta() << "   phi = " << jet.phi() << std::endl;
            			std::cout << "Pt (BTagCollection)  = " << bTag[i].first->pt() << std::endl;
            			std::cout << "bTag                 = " << bTag[i].second << std::endl;
          			}
        		} else continue;
    		}
    		for (unsigned i = 0; i != bbTag.size(); i++) {
        		if (abs(bbTag[i].first->eta() - jet.eta()) < 0.025 && abs(bbTag[i].first->phi() - jet.phi()) < 0.025 ){
          			JetDeepCSV_prob_bb = bbTag[i].second;
          			if (monitoring) {
            			std::cout << "bbTag                = " << bbTag[i].second << std::endl;
          			}
        		} else continue;
    		}

    		jetPtSum15 += jet.pt();
			if (JetFromPV) {
				jetPtSum15PV += jet.pt();
				if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_loose) {
    				nLooseBtagedJetsPV++;
    				if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_medium) {
    					nMediumBtagedJetsPV++;
    					if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_tight) {
    						nTightBtagedJetsPV++;
    					}
    				}
    			}
			}
			if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_loose) {
    			nLooseBtagedJets++;
    			if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_medium) {
    				nMediumBtagedJets++;
    				if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_tight) {
    					nTightBtagedJets++;
    				}
    			}
    		}

			if (jet.pt() > 20) {
				jetPtSum20 += jet.pt();
				nJets20++;
				if (JetFromPV) {
					nJets20PV++;
					jetPtSum20PV += jet.pt();
				}
			}
		}
	}

	if (monitoring) std::cout << "CHSJets:"<< std::endl;

	int j = 0;
	// Lopp over Jets
	for (auto& jet: *CHSjets) {
		if (TMath::Abs(jet.eta()) > 3) continue;
		if (jet.pt() > 15) {
			j++;
			if (monitoring) {
				std::cout << std::endl;
				std::cout << "CHSJet[" << j << "] Pt = " << jet.pt() << std::endl;
			}
			reco::TrackRefVector JetTracksRef = jet.getTrackRefs();
			std::map <int, int> VertexTracksMap;
			std::map <int, std::pair<int, double>> VertexTracksPtMap;
			double JetDeepCSV_prob_b = 0;
			double JetDeepCSV_prob_bb = 0;
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
				// Loop over vrtices
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
			if (monitoring) std::cout << std::endl;
			for (auto Map_iter = VertexTracksMap.begin(); Map_iter != VertexTracksMap.end(); ++Map_iter) {
				if ((*Map_iter).first == 0) NtracksFromPV = (*Map_iter).second;
				if ((*Map_iter).first > 0 && (*Map_iter).second > NtracksFromPV) JetFromPV = false;
				if ((*Map_iter).second < 1) continue;
			}
			for (auto Map_iter1 = VertexTracksPtMap.begin(); Map_iter1 != VertexTracksPtMap.end(); ++Map_iter1) {
				if ((*Map_iter1).second.first < 1) continue;
				if (monitoring) std::cout << "Vertex number   = " << (*Map_iter1).first << "  entries = " << (*Map_iter1).second.first << "  SumPt = " << (*Map_iter1).second.second << std::endl;
			}

			// Loop over bTags
    		for (unsigned i = 0; i != bTag.size(); i++) {
    			if (abs(bTag[i].first->eta() - jet.eta()) < 0.001 && abs(bTag[i].first->phi() - jet.phi()) < 0.001 ){
          			JetDeepCSV_prob_b = bTag[i].second;
          			if (monitoring) {
            			std::cout << "##### Btag associated with jet number " << std::endl;
            			std::cout << "BTagCollection   eta = " << bTag[i].first->eta() << "   phi = " << bTag[i].first->phi() << std::endl;
            			std::cout << "PFCollection     eta = " << jet.eta() << "   phi = " << jet.phi() << std::endl;
            			std::cout << "Pt (BTagCollection)  = " << bTag[i].first->pt() << std::endl;
            			std::cout << "bTag                 = " << bTag[i].second << std::endl;
          			}
        		} else continue;
    		}
    		for (unsigned i = 0; i != bbTag.size(); i++) {
        		if (abs(bbTag[i].first->eta() - jet.eta()) < 0.001 && abs(bbTag[i].first->phi() - jet.phi()) < 0.001 ){
          			JetDeepCSV_prob_bb = bbTag[i].second;
          			if (monitoring) {
            			std::cout << "bbTag                = " << bbTag[i].second << std::endl;
          			}
        		} else continue;
    		}

    		if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > BJet_btag) {
    			BJet_pt   = jet.pt();
	            BJet_eta  = jet.eta();
				BJet_phi  = jet.phi();
				BJet_m    = jet.mass();
				BJet_btag = JetDeepCSV_prob_b + JetDeepCSV_prob_bb;
    		}

    		if (j == 1) {
    		    LeadingJet_pt   = jet.pt();
    		    LeadingJet_eta  = jet.eta();
    		    LeadingJet_phi  = jet.phi();
    		    LeadingJet_m    = jet.mass();
    		    LeadingJet_btag = JetDeepCSV_prob_b + JetDeepCSV_prob_bb;
    		} else if (j == 2) {
    		    SubLeadingJet_pt   = jet.pt();
    		    SubLeadingJet_eta  = jet.eta();
    		    SubLeadingJet_phi  = jet.phi();
    		    SubLeadingJet_m    = jet.mass();
    		    SubLeadingJet_btag = JetDeepCSV_prob_b + JetDeepCSV_prob_bb;
    		}
			
			CHSjetPtSum15 += jet.pt();
			if (JetFromPV) {
				CHSjetPtSum15PV += jet.pt();
				if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_loose) {
    				nLooseBtagedCHSJetsPV++;
    				if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_medium) {
    					nMediumBtagedCHSJetsPV++;
    					if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_tight) {
    						nTightBtagedCHSJetsPV++;
    					}
    				}
    			}
			}
			if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_loose) {
    			nLooseBtagedCHSJets++;
    			if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_medium) {
    				nMediumBtagedCHSJets++;
    				if (JetDeepCSV_prob_b + JetDeepCSV_prob_bb > WPBTag_tight) {
    					nTightBtagedCHSJets++;
    				}
    			}
    		}

			if (jet.pt() > 20) {
				CHSjetPtSum20 += jet.pt();
				nCHSJets20++;
				if (JetFromPV) {
					nCHSJets20PV++;
					CHSjetPtSum20PV += jet.pt();
				}
			}
		}
	}

	return nJets20 > 0;
};

void TreeMaker::CountTracks(const edm::Event& event) {
	nTrks = 0;
	edm::Handle<reco::TrackCollection> tracks;
	event.getByToken(TrackToken_, tracks);
	if (!tracks.isValid()) return;
	nTrks = tracks->size();
};

void TreeMaker::AddWT(const edm::Event& iEvent) {

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
	edm::Handle<bool>   PhotonEmission1;
	iEvent.getByToken(Photon1Token_, PhotonEmission1);
	edm::Handle<bool>   PhotonEmission2;
	iEvent.getByToken(Photon2Token_, PhotonEmission2);

	if (!WTisValidHandle.isValid()) {
		//if (monitoring) std::cout << "WTisValidHandle is not valid" << std::endl;
		return;
	}
	WTisValid = *WTisValidHandle;
	if (WTisValid) {
		WT        = *WTHandle;
		WTFlip    = *WTFlipHandle;
		WThminus  = *WThminusHandle;
		WThplus   = *WThplusHandle;
		TauSpinnerMother = *TauMotherHandle;
		PhotonFlag1 = *PhotonEmission1;
		PhotonFlag2 = *PhotonEmission2;
	} else {
		//if (monitoring) std::cout << "WT Collections are not valid" << std::endl;
		WT        = null;
		WTFlip    = null;
		WThminus  = null;
		WThplus   = null;
		TauSpinnerMother = null;
		PhotonFlag1 = null;
		PhotonFlag2 = null;
	}

};

void TreeMaker::AddWTData(const edm::Event& iEvent) {

	int One = 1;
	WTisValid = One;
	WT        = One;
	WTFlip    = One;
	WThminus  = One;
	WThplus   = One;
	TauSpinnerMother = null;
}

void TreeMaker::GenTauDM (const edm::Event& iEvent) {

  edm::Handle<reco::GenParticleCollection> GP;
  iEvent.getByToken(GenParticleToken_, GP);

  int channel = null;
  gentau_dm = null;
  gentau_nPi0 = null;

  if(GP.isValid()) {
    for(unsigned i = 0 ; i < GP->size() ; i++) {
      if (abs((*GP)[i].pdgId())==15) {
        channel = -1;
	int tau_pdgid = (*GP)[i].pdgId();
        if (monitoring) std::cout << "GenTauDM particle number " << i << std::endl;

        //std::vector<SimpleParticle> tau_daughters_simple;
	std::vector<MySimpleParticle> tau_daughters;
        for (unsigned k = 0; k < (*GP)[i].numberOfDaughters(); k++) {
          MySimpleParticle tp((*GP)[i].daughter(k)->p4().Px(), (*GP)[i].daughter(k)->p4().Py(), (*GP)[i].daughter(k)->p4().Pz(), (*GP)[i].daughter(k)->p4().E(), (*GP)[i].daughter(k)->pdgId());
          tau_daughters.push_back(tp);
          if (monitoring) std::cout << "Tau_daughter[" << k << "] = " << (*GP)[i].daughter(k)->pdgId() << std::endl;
        }

	if (tau_daughters.size() >= 2 && (abs(tau_daughters[0].pdgid()) == 15 || abs(tau_daughters[1].pdgid()) == 15)) continue;

        std::vector<int>  pdgid;
        for (unsigned int l = 0; l < tau_daughters.size(); l++) {
          pdgid.push_back(tau_daughters[l].pdgid());
        }
        if( pdgid.size() == 2 &&
            (
              ( tau_pdgid == 15 && DecaychannelMatch(tau_daughters, 16,-211) ) ||
              ( tau_pdgid ==-15 && DecaychannelMatch(tau_daughters,-16, 211) ) ||
              ( tau_pdgid == 15 && DecaychannelMatch(tau_daughters, 16,-321) ) ||
              ( tau_pdgid ==-15 && DecaychannelMatch(tau_daughters,-16, 321) )
            )
          ) {
          channel = 3;
	  gentau_dm = 0;
	  gentau_nPi0 = 0;
          if (abs(pdgid[1]) == 321) channel = 6;
        }
        else if ( pdgid.size() == 3 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211, 111) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 111) ) ||
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211, 130) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 130) ) ||
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211, 310) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 310) ) ||
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-321, 111) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 321, 111) )
            )
          ) {
 
          channel = 4;
	  gentau_dm = 1;
	  gentau_nPi0 = 1;
          if(abs(pdgid[1])==321 || abs(pdgid[2])==130 || abs(pdgid[2])==310) channel = 7;
        }
        else if( pdgid.size() == 4 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111,-211) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 211) )
            )
          ) {
          channel = 5;
	  gentau_dm = 2;
	  gentau_nPi0 = 2;
        }
        else if( pdgid.size() == 4 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211,-211, 211) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 211,-211) )
            )
          ) {
          channel = 6;
	  gentau_dm = 10;
	  gentau_nPi0 = 0;
        }
        else if( pdgid.size()==5 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16,-211,-211, 211, 111) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 211, 211,-211, 111) )
            )
          ) {
          channel = 8;
	  gentau_dm = 11;
	  gentau_nPi0 = 1;
        }
        else if( pdgid.size()==5 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 111, 111, 111,-211) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 111, 211) )
            )
          ) {
          channel = 9;
	  gentau_nPi0 = 3;
        }
        else if( pdgid.size()==3 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 11,-12) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-11, 12) )
            )
          ) {
          channel = 1;
        }
        else if( pdgid.size()==4 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 11,-12, 22) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-11, 12, 22) )
            )
          ) {
          channel = 1;
        }
        else if( pdgid.size()==3 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 13,-14) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-13, 14) )
            )
          ) {
          channel = 2;
        }
        else if( pdgid.size()==4 &&
            (
              ( tau_pdgid== 15 && DecaychannelMatch(tau_daughters, 16, 13,-14, 22) ) ||
              ( tau_pdgid==-15 && DecaychannelMatch(tau_daughters,-16,-13, 14, 22) )
            )
          ) {
          channel = 2;
        }
      }
      if (channel == -1 && abs((*GP)[i].pdgId())==15 && monitoring) {
        std::cout << "Unidentified decay channel of tau was found" << std::endl;
	    std::cout << "List of daughters:" << std::endl;
        for (unsigned t = 0; t < (*GP)[i].numberOfDaughters(); t++) {
          std::cout << "TauDaughter[" << t << "] pdgId = " << (*GP)[i].daughter(t)->pdgId() << std::endl;
        }
      }
    }
  }
  GenHadronDecayChannel = channel;
  if (GP.isValid() && GenHadronDecayChannel == null && tau_found == 1 && monitoring) {
    std::cout << "There are no tau leptons among Gen particles" << std::endl;
    std::cout << "tau_eta    = " << tau_eta << "   " << "tau_phi    = " << tau_phi << std::endl;
    for(unsigned j = 0 ; j < GP->size() ; j++) {
      if (abs(tau_eta - (*GP)[j].eta()) < 0.1 && abs(tau_phi - (*GP)[j].phi()) < 0.1) {
        std::cout << "GP[" << j << "] pdgId = " << (*GP)[j].pdgId() << std::endl;
      }
    }
    std::cout << "piChar_eta = " << piChar_eta << "   " << "piChar_phi = " << piChar_phi << std::endl;
    for(unsigned j = 0 ; j < GP->size() ; j++) {
      if (abs(piChar_eta - (*GP)[j].eta()) < 0.1 && abs(piChar_phi - (*GP)[j].phi()) < 0.1) {
        std::cout << "GP[" << j << "] pdgId = " << (*GP)[j].pdgId() << std::endl;
      }
    }
    std::cout << "piZero_eta = " << piZero_eta << "   " << "piZero_phi = " << piZero_phi << std::endl;
    for(unsigned j = 0 ; j < GP->size() ; j++) {
      if (abs(piZero_eta - (*GP)[j].eta()) < 0.1 && abs(piZero_phi - (*GP)[j].phi()) < 0.1) {
        std::cout << "GP[" << j << "] pdgId = " << (*GP)[j].pdgId() << std::endl;
      }
    }
  }
}

bool TreeMaker::DecaychannelMatch(std::vector<MySimpleParticle> &particles, int p1, int p2, int p3, int p4, int p5, int p6) {
  // Copy pdgid-s of all particles
  std::vector<int> list;
  
  // Fill vector list with PGDId-s of tau daughters
  for (unsigned int i = 0; i < particles.size(); i++) list.push_back(particles[i].pdgid());
  
  // Create array out of pdgid-s
  int p[6] = { p1, p2, p3, p4, p5, p6 };

  // 1) Check if 'particles' contain all pdgid-s on the list 'p'
  for (int i = 0; i < 6; i++) {
    // if the pdgid is zero - finish (the list of searched particles ends)
    if(p[i] == 0) break;
    bool found = false;
    // Loop over tau daughters
    for (unsigned int j = 0;j < list.size(); j++) {
      // if pdgid is found - erese it from the list and search for the next one
      if(list[j] == p[i]) {
        found = true;
        list.erase(list.begin()+j);
        break;
      }
    }
    
    if(!found) return false;
  }
  
  // if there are more particles on the list - there is no match
  if(list.size()!=0) return false;

  
  // 2) Rearrange particles to match the order of pdgid-s listed in array 'p'

  std::vector<MySimpleParticle> newList;
  
  for(int i = 0; i < 6; i++) {
    // if the pdgid is zero - finish
    if(p[i] == 0) break;
    for(unsigned int j = 0; j < particles.size(); j++) {
      // if pdgid is found - copy it to new list and erese from the old one
      if(particles[j].pdgid() == p[i]) {
        newList.push_back(particles[j]);
        particles.erase(particles.begin()+j);
        break;
      }
    }
  }
  
  particles = newList;

  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.add("treemaker", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeMaker);
