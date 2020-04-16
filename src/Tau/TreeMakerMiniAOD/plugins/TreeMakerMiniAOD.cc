// -*- C++ -*-
//
// Package:    Tau/TreeMakerMiniAOD
// Class:      TreeMakerMiniAOD
//
/**\class TreeMakerMiniAOD TreeMaker.cc Tau/TreeMaker/plugins/TreeMaker.cc

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

#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"


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

class TreeMakerMiniAOD : public edm::EDAnalyzer {
public:
	explicit TreeMakerMiniAOD(const edm::ParameterSet&);
	~TreeMakerMiniAOD();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;
	
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

	virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	
	// ----------member data ---------------------------

	edm::EDGetTokenT<pat::TauCollection> TauCollectionToken_;
	edm::EDGetTokenT<pat::MuonCollection> MuonCollectionToken_;
	edm::EDGetTokenT<pat::ElectronCollection> ElectronCollectionToken_;
	edm::EDGetTokenT<pat::JetCollection> JetCollectionToken_;
	edm::EDGetTokenT<pat::JetCollection> PuppiJetCollectionToken_;
	edm::EDGetTokenT<pat::METCollection> MetCollectionToken_;
	edm::EDGetTokenT<reco::VertexCollection> PVToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> GenParticleToken_;
	
	edm::EDGetTokenT<pat::IsolatedTrackCollection> TrackToken_;

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
	double tau_mediumElectronRejection;
	double tau_tightElectronRejection;
	

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

	double DiPhoton_pt;
	double DiPhoton_eta;
	double DiPhoton_phi;
	double DiPhoton_m;
	
	double pipiMass;
	double upsilon;

    double jetPtSum15, jetPtSum20, jetPtSum15PV, jetPtSum20PV;
	int nJets20, nJets20PV;
	double PuppijetPtSum15, PuppijetPtSum20, PuppijetPtSum15PV, PuppijetPtSum20PV;
	int nPuppiJets20, nPuppiJets20PV;
	int nLooseBtagedPuppiJets, nMediumBtagedPuppiJets, nTightBtagedPuppiJets;
	int nLooseBtagedPuppiJetsPV, nMediumBtagedPuppiJetsPV, nTightBtagedPuppiJetsPV;
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
	int nPhotons;
	double met;
	double met_phi;
	double met_eta;
	double met_significance;
	double met_mEtSig;
	
	
	double m_t;
	
	double dPhi;
	
	math::XYZPoint pv_position;

	double WT;
	double WTFlip;
	double WThminus;
	double WThplus;
	bool WTisValid;
	int  TauSpinnerMother;

	//////////////////////////////////////////////////////
	bool isMC;
	double tauPtMin;
	double piPtMin;
	double tauEtaMax;
	double tauDzMax;
	double null;
	bool monitoring;

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
TreeMakerMiniAOD::TreeMakerMiniAOD(const edm::ParameterSet& iConfig) {
	//now do what ever initialization is needed
	iT =0;

	isMC						= iConfig.getParameter<bool>("isMC");
	monitoring					= iConfig.getParameter<bool>("monitoring");
	tauPtMin 					= iConfig.getParameter<double>("tauPtMin");
	piPtMin 					= iConfig.getParameter<double>("piPtMin");
	tauEtaMax 					= iConfig.getParameter<double>("tauEtaMax");
	tauDzMax 					= iConfig.getParameter<double>("tauDzMax");
	null 						= iConfig.getParameter<double>("null");

	std::string tauCollection         = iConfig.getParameter<std::string>("tauCollection");
	std::string muonCollection        = iConfig.getParameter<std::string>("muonCollection");
	std::string electronCollection    = iConfig.getParameter<std::string>("electronCollection");
	std::string jetCollection         = iConfig.getParameter<std::string>("jetCollection");
	std::string PuppijetCollection    = iConfig.getParameter<std::string>("PuppijetCollection");
	std::string metCollection         = iConfig.getParameter<std::string>("metCollection");
	std::string vertexCollection      = iConfig.getParameter<std::string>("vertexCollection");
	std::string genParticleCollection = iConfig.getParameter<std::string>("genParticleCollection");
	std::string trackCollection       = iConfig.getParameter<std::string>("trackCollection");
	
	TauCollectionToken_ 		= consumes<pat::TauCollection>(edm::InputTag(tauCollection));
	MuonCollectionToken_ 		= consumes<pat::MuonCollection>(edm::InputTag(muonCollection));
	ElectronCollectionToken_	= consumes<pat::ElectronCollection>(edm::InputTag(electronCollection));
	JetCollectionToken_ 		= consumes<pat::JetCollection>(edm::InputTag(jetCollection));
	PuppiJetCollectionToken_        = consumes<pat::JetCollection>(edm::InputTag(PuppijetCollection));
	MetCollectionToken_ 		= consumes<pat::METCollection>(edm::InputTag(metCollection));
	PVToken_ 			= consumes<reco::VertexCollection>(edm::InputTag(vertexCollection));
	GenParticleToken_ 		= consumes<reco::GenParticleCollection>(edm::InputTag(genParticleCollection));
	TrackToken_			= consumes<pat::IsolatedTrackCollection>(edm::InputTag(trackCollection));

	if (isMC) {
		TauSpinnerWTToken_        = consumes<double>(iConfig.getParameter<edm::InputTag>("WTCollection"));
		TauSpinnerWTFlipToken_    = consumes<double>(iConfig.getParameter<edm::InputTag>("WTFlipCollection"));
		TauSpinnerWThminusToken_  = consumes<double>(iConfig.getParameter<edm::InputTag>("WThminusCollection"));
		TauSpinnerWThplusToken_   = consumes<double>(iConfig.getParameter<edm::InputTag>("WThplusCollection"));
		TauSpinnerWTisValidToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("WTisValidCollection"));
		TauSpinnerMotherToken_    = consumes<int>(iConfig.getParameter<edm::InputTag>("MotherCollection"));
	}

}


TreeMakerMiniAOD::~TreeMakerMiniAOD() {
	 // do anything here that needs to be done at desctruction time
	 // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void TreeMakerMiniAOD::analyze(const edm::Event& event, const edm::EventSetup&) {
	t_Run   = event.id().run();
	t_Event = event.id().event();

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

	TLorentzVector pi0, pi1;
	pi0.SetPtEtaPhiM(piZero_pt, piZero_eta, piZero_phi, piZero_m);
	pi1.SetPtEtaPhiM(piChar_pt, piChar_eta, piChar_phi, piChar_m);

	pipiMass = (pi0 + pi1).M();
//	upsilon  = (pi1.E() - pi0.E()) / tau_pt;
	upsilon  = 2 * piChar_pt / tau_pt - 1;
	dPhi	 = dphi(tau_phi, met_phi);
	m_t      = sqrt(2 * tau_pt * met * (1 - cos(dPhi)));

	tree->Fill();
};


// ------------ method called once each job just before starting event loop  ------------
void TreeMakerMiniAOD::beginJob() {
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
	tree->Branch("tau_mediumElectronRejection",&tau_mediumElectronRejection,"tau_mediumElectronRejection/D");
	tree->Branch("tau_tightElectronRejection",&tau_tightElectronRejection,"tau_tightElectronRejection/D");
	
	tree->Branch("piChar_pt", &piChar_pt, "piChar_pt/D");
	tree->Branch("piChar_eta", &piChar_eta, "piChar_eta/D");
	tree->Branch("piChar_phi", &piChar_phi, "piChar_phi/D");
	tree->Branch("piChar_q", &piChar_q, "piChar_q/D");
	tree->Branch("piChar_m", &piChar_m, "piChar_m/D");

	tree->Branch("piZero_pt", &piZero_pt, "piZero_pt/D");
	tree->Branch("piZero_eta", &piZero_eta, "piZero_eta/D");
	tree->Branch("piZero_phi", &piZero_phi, "piZero_phi/D");
	tree->Branch("piZero_m", &piZero_m, "piZero_m/D");

	tree->Branch("DiPhoton_pt", &DiPhoton_pt, "DiPhoton_pt/D");
	tree->Branch("DiPhoton_eta", &DiPhoton_eta, "DiPhoton_eta/D");
	tree->Branch("DiPhoton_phi", &DiPhoton_phi, "DiPhoton_phi/D");
	tree->Branch("DiPhoton_m", &DiPhoton_m, "DiPhoton_m/D");

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

	// PuppiJets
	tree->Branch("PuppijetPtSum15", &PuppijetPtSum15, "PuppijetPtSum15/D");
	tree->Branch("PuppijetPtSum20", &PuppijetPtSum20, "PuppijetPtSum20/D");
	tree->Branch("PuppijetPtSum15PV", &PuppijetPtSum15PV, "PuppijetPtSum15PV/D");
	tree->Branch("PuppijetPtSum20PV", &PuppijetPtSum20PV, "PuppijetPtSum20PV/D");
	tree->Branch("nPuppiJets20", &nPuppiJets20, "nPuppiJets20/I");
	tree->Branch("nPuppiJets20PV", &nPuppiJets20PV, "nPuppiJets20PV/I");
	tree->Branch("nLooseBtagedPuppiJets", &nLooseBtagedPuppiJets, "nLooseBtagedPuppiJets/I");
	tree->Branch("nMediumBtagedPuppiJets", &nMediumBtagedPuppiJets, "nMediumBtagedPuppiJets/I");
	tree->Branch("nTightBtagedPuppiJets", &nTightBtagedPuppiJets, "nTightBtagedPuppiJets/I");
	tree->Branch("nLooseBtagedPuppiJetsPV", &nLooseBtagedPuppiJetsPV, "nLooseBtagedPuppiJetsPV/I");
	tree->Branch("nMediumBtagedPuppiJetsPV", &nMediumBtagedPuppiJetsPV, "nMediumBtagedPuppiJetsPV/I");
	tree->Branch("nTightBtagedPuppiJetsPV", &nTightBtagedPuppiJetsPV, "nTightBtagedPuppiJetsPV/I");

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
	tree->Branch("nPhotons",&nPhotons,"nPhotons/I");

	tree->Branch("WT", &WT, "WT/D");
	tree->Branch("WTFlip", &WTFlip, "WTFlip/D");
	tree->Branch("WThminus", &WThminus, "WThminus/D");
	tree->Branch("WThplus", &WThplus, "WThplus/D");
	tree->Branch("TauSpinnerMother", &TauSpinnerMother, "TauSpinnerMother/I");
	//tree->Branch("WTisValid", &WTisValid, "WTisValid/D")
	// add more branches
	
	allTauPt = FS->make<TH1F>("allTauPt","allTauPt",300,0,300);
}

// ------------ method called once each job just after ending the event loop  ------------
void TreeMakerMiniAOD::endJob() {

}

// ------------ method called when starting to processes a run  ------------

void TreeMakerMiniAOD::beginRun(const edm::Run& Run, const edm::EventSetup& Setup) {
}

// ------------ method called when ending the processing of a run  ------------
/*
void
TreeMakerMiniAOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TreeMakerMiniAOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TreeMakerMiniAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

//-----------------------------------------------------------------------------------------

bool TreeMakerMiniAOD::AddTau(const edm::Event& event) {
	tau_found = 0;
	nTauC     = 0;
	math::XYZTLorentzVector DiPhoton_p4;
	int nPhotons_temp = 0;

	edm::Handle<pat::TauCollection> taus;
	event.getByToken(TauCollectionToken_, taus);
	if (!taus.isValid() || taus->size() == 0) return false;

	// search for the tau candidate with the minimum isolation and the
	// maximum transverse momentum
	size_t index  = taus->size();;
	for (size_t i = 0; i < taus->size(); ++i) {
#define cut(condition) if (!(condition)) continue;
		auto& tau = (*taus)[i];
		cut(tau.pt() > tauPtMin); // tau transverse momentum
		cut(TMath::Abs(tau.eta()) < tauEtaMax); // tau pseudorapidity
		//cut((*taus)[i].tauID("byLooseIsolationMVArun2v1DBoldDMwLT")); // at least loose Iso

		//auto& chargedHadrons = tau.signalPFChargedHadrCands();
		//std::cout << "NPi+ = " << chargedHadrons.size() << std::endl;
		//cut(chargedHadrons.size() == 1); // single charged hadron (pi+)
		//auto& pi0s = tau.signalPiZeroCandidates();
		//const std::vector<reco::RecoTauPiZero>& pi0s = tau.signalPiZeroCandidates();
		//std::cout << "NPi0 = " << pi0s.size() << std::endl;
		//cut(pi0s.size() > 0); // at least one pi0

		//cut(tau.decayMode() == 1);
		if (monitoring) {
		  std::cout << "##### Tau number " << i << " ######" << std::endl;
		  std::cout << "Discriminator = " << (*taus)[i].tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") << std::endl;
		  std::cout << "LooseMVA = " << (*taus)[i].tauID("byLooseIsolationMVArun2v1DBoldDMwLT") << std::endl;
		}
		//tau.embedSignalPFChargedHadrCands();
		//tau.embedSignalPFNeutralHadrCands();
		//tau.embedLeadPFNeutralCand();

		const reco::CandidatePtr leadPiCh = tau.leadChargedHadrCand();
		const reco::CandidatePtr leadPi0  = tau.leadNeutralCand();
		reco::CandidatePtrVector VectorPi0 = tau.signalNeutrHadrCands();
		reco::CandidatePtrVector VectorPiCh = tau.signalChargedHadrCands();
		reco::CandidatePtrVector VectorSignalCands = tau.signalCands(); 

		cut(VectorPiCh.size() == 1); // single charged hadron (pi+)

		const std::vector<reco::PFCandidatePtr>& VectorNeutralPF = tau.signalPFNeutrHadrCands();
		//const std::vector< reco::RecoTauPiZero>& VectorPi0PF = tau.signalPiZeroCandidates();

		//reco::CandidatePtrVector NeutralHadrCands =	tau.signalNeutralHadrCandPtrs_;

        //reco::CandidatePtrVector VectorPi0_v1 = tau.isolationNeutrHadrCands();
		//reco::CandidatePtrVector VectorPiCh_v1 = tau.isolationChargedHadrCands();

		//const std::vector<reco::PFCandidatePtr>& isolationPFNeutrHadrCands();

		//if (monitoring) std::cout << "isPFTau    = " << tau.isPFTau() << std::endl;
		if (monitoring) std::cout << "Decay mode = " << tau.decayMode() << std::endl;
		if (monitoring) std::cout << "tau Pt     = " << tau.pt() << std::endl;
		if (monitoring) std::cout << "Size of PiCh vector = " << VectorPiCh.size() << std::endl;
		//std::cout << "Size of PiCh_v1 vector = " << VectorPiCh_v1.size() << std::endl;
		if (monitoring) std::cout << "Size of Pi0 vector  = " << VectorPi0.size() << std::endl;
		//std::cout << "Size of Pi0_v1 vector  = " << VectorPi0_v1.size() << std::endl;
		if (monitoring) std::cout << "Size of SignalCands vector  = " << VectorSignalCands.size() << std::endl;
		//if (monitoring) std::cout << "ecalStripSumEOverPLead  = " << tau.ecalStripSumEOverPLead() << std::endl;
		//std::cout << "Size of SignaNeutralHadrCands vector  = " << VectorSignalCands.size() << std::endl;
		//std::cout << "Size of VectorNeutralPF  = " << VectorNeutralPF.size() << std::endl;
		//std::cout << "Size of VectorPi0PF      = " << VectorPi0PF.size() << std::endl;

        if (monitoring) {
          std::cout << "Pt_PiCh_lead   = " << leadPiCh->pt() << std::endl;
		  std::cout << "eta_PiCh_lead  = " << leadPiCh->eta() << std::endl;
		  std::cout << "phi_PiCh_lead  = " << leadPiCh->phi() << std::endl;
		  std::cout << "mass_PiCh_lead = " << leadPiCh->mass() << std::endl;
        }
		

		if (VectorPi0.size() > 0 && monitoring) {
		  std::cout << "Pt_Zero_lead   = " << VectorPi0[0]->pt() << std::endl;
		  std::cout << "eta_Zero_lead  = " << VectorPi0[0]->eta() << std::endl;
		  std::cout << "phi_Zero_lead  = " << VectorPi0[0]->phi() << std::endl;
		  std::cout << "mass_Zero_lead = " << VectorPi0[0]->mass() << std::endl;
		}

		if (VectorSignalCands.size() > 0 && monitoring) {
		  for (unsigned l = 0; l < VectorSignalCands.size(); l++) {
		  	std::cout << "Signal Candidate number " << l << std::endl;
		  	std::cout << "pdgId = " << VectorSignalCands[l]->pdgId() << std::endl;
			std::cout << "Pt    = " << VectorSignalCands[l]->pt() << std::endl;
			std::cout << "eta   = " << VectorSignalCands[l]->eta() << std::endl;
			std::cout << "phi   = " << VectorSignalCands[l]->phi() << std::endl;
			std::cout << "mass  = " << VectorSignalCands[l]->mass() << std::endl;
		  }
		}

		reco::CandidatePtrVector VectorSignalCands_Photons;
		VectorSignalCands_Photons.clear();

		if (VectorSignalCands.size() > 0) {
		  for (unsigned l = 0; l < VectorSignalCands.size(); l++) {
		  	if (VectorSignalCands[l]->pdgId() == 22) {
		  	  VectorSignalCands_Photons.push_back(VectorSignalCands[l]);
		  	}
		  }
		}

		//cut(VectorSignalCands_Photons.size() <= 2);

		math::XYZTLorentzVector Photons_p4;
		//math::XYZTLorentzVector DiPhoton_p4;
		double InvaMass_temp = null;

		if (VectorSignalCands_Photons.size() > 0) {
		  for (unsigned l = 0; l < VectorSignalCands_Photons.size(); l++) {
		  	Photons_p4 += VectorSignalCands_Photons[l]->p4();
		  	for (unsigned n = l; n < VectorSignalCands_Photons.size(); n++) {
		  	  if (n == l) continue;
		  	  double InvaMass_pi0 = (VectorSignalCands_Photons[l]->p4() + VectorSignalCands_Photons[n]->p4()).M();
		  	  double Pt_pi0       = (VectorSignalCands_Photons[l]->p4() + VectorSignalCands_Photons[n]->p4()).Pt();
		  	  if (monitoring) std::cout << "mass Photon [" << l <<"][" << n << "] = " << InvaMass_pi0 << std::endl;
		  	  if (l == 0 && n == 1) {
		  	  	InvaMass_temp = InvaMass_pi0;
		  	  	DiPhoton_p4 = VectorSignalCands_Photons[l]->p4() + VectorSignalCands_Photons[n]->p4();
		  	  } else if (abs(InvaMass_pi0 - 0.135) < abs(InvaMass_temp - 0.135)) {
		  	  	InvaMass_temp = InvaMass_pi0;
		  	  	DiPhoton_p4 = VectorSignalCands_Photons[l]->p4() + VectorSignalCands_Photons[n]->p4();
		  	  }
		  	}
		  }
		}

		nPhotons_temp = VectorSignalCands_Photons.size();

		if (VectorPiCh.size() > 0 && monitoring) {
		  for (unsigned l = 0; l < VectorPiCh.size(); l++) {
		  	std::cout << "PiCh number " << l << std::endl;
			std::cout << "Pt_PiCh   = " << VectorPiCh[l]->pt() << std::endl;
			std::cout << "eta_PiCh  = " << VectorPiCh[l]->eta() << std::endl;
			std::cout << "phi_PiCh  = " << VectorPiCh[l]->phi() << std::endl;
			std::cout << "mass_PiCh = " << VectorPiCh[l]->mass() << std::endl;
		  }
		}

		/*
		if (tau.decayMode() == 1) {
			std::cout << "SumPhotons_M   = " << Photons_p4.M() << std::endl;
			std::cout << "SumPhotons_Pt  = " << Photons_p4.Pt() << std::endl;
			std::cout << "SumPhotons_Eta = " << Photons_p4.Eta() << std::endl;
			std::cout << "SumPhotons_Phi = " << Photons_p4.Phi() << std::endl;
		}
		*/

		math::XYZTLorentzVector tau_p4_temp = tau.p4();
    	math::XYZTLorentzVector piChar_p4_temp = tau.leadChargedHadrCand()->p4();

		if (monitoring) {
          std::cout << "Calculated Pi0 parameters:" << std::endl;
          std::cout << "PiZero_m   = " << (tau_p4_temp - piChar_p4_temp).M() << std::endl;
	      std::cout << "PiZero_pt  = " << (tau_p4_temp - piChar_p4_temp).pt() << std::endl;
	      std::cout << "PiZero_eta = " << (tau_p4_temp - piChar_p4_temp).eta() << std::endl;
	      std::cout << "PiZero_phi = " << (tau_p4_temp - piChar_p4_temp).phi() << std::endl;
	      std::cout << "DiPhoton parameters:" << std::endl;
          std::cout << "DiPhoton_m   = " << DiPhoton_p4.M() << std::endl;
	      std::cout << "DiPhoton_pt  = " << DiPhoton_p4.pt() << std::endl;
	      std::cout << "DiPhoton_eta = " << DiPhoton_p4.eta() << std::endl;
	      std::cout << "DiPhoton_phi = " << DiPhoton_p4.phi() << std::endl;
    	}
		
		//auto& pi0 = pi0s.front();
		//auto& pi1 = *chargedHadrons.front();

		//cut(pi1.pt() > piPtMin); // pi+ transverse momentum

		++nTauC;
		allTauPt->Fill(tau.pt());

		cut((pv_position - tau.vertex()).R() < tauDzMax); // tau vertex displacement

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

	if (index == taus->size()) return false;

	auto& tau = (*taus)[index];
	const reco::CandidatePtr leadPiCh = tau.leadChargedHadrCand();
	const reco::CandidatePtr leadPi0  = tau.leadNeutralCand();
	//auto& pi0 = tau.signalPiZeroCandidates().front();
	//auto& pi1 = *tau.signalPFChargedHadrCands().front();

    //reco::CandidateCollection Tau_daughters = tau.daughters;


	tau_pt                     = tau.pt();
	tau_eta                    = tau.eta();
	tau_phi                    = tau.phi();
	tau_dm                     = tau.decayMode();
	tau_dz                     = (pv_position - tau.vertex()).R();
	tau_q                      = tau.charge();
	tau_m                      = tau.mass();
	
	if (monitoring) std::cout << "List of pat::tau discriminators (tauID):" << std::endl;
	if (monitoring) std::cout << "Selected tau discriminator (" << index << ") = " << tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") << std::endl;
	tau_absIso                  = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
	tau_looseCombinedIso        = tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
	tau_mediumCombinedIso       = tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
	tau_tightCombinedIso        = tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
	tau_looseMvaIso             = tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT"); // there are more possible discriminators 
	tau_mediumMvaIso            = tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
	tau_tightMvaIso             = tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
	tau_looseMuonRejection      = tau.tauID("againstMuonLoose3");
	tau_tightMuonRejection      = tau.tauID("againstMuonTight3");
	tau_looseElectronRejection  = tau.tauID("againstElectronLooseMVA6");
	tau_mediumElectronRejection = tau.tauID("againstElectronMediumMVA6");
	tau_tightElectronRejection  = tau.tauID("againstElectronTightMVA6");

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

	DiPhoton_pt  = DiPhoton_p4.pt();
	DiPhoton_eta = DiPhoton_p4.eta();
	DiPhoton_phi = DiPhoton_p4.phi();
	DiPhoton_m   = DiPhoton_p4.M();

	nPhotons = nPhotons_temp;

    /*
	piZero_pt  = leadPi0->pt();
	piZero_eta = leadPi0->eta();
	piZero_phi = leadPi0->phi();
	piZero_m   = leadPi0->mass();
	
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
	*/
	nPi0 = 1;
	nTau      = taus->size();
	tau_found = 1;

	return true;
};

bool TreeMakerMiniAOD::FindGenTau(const edm::Event& event) {
	gentau_found  = 0;
	genTauFromW   = null;
	genTauMother  = null;
	dR            = null;
	W_pt          = null;

	edm::Handle<reco::GenParticleCollection> genParticles;
	event.getByToken(GenParticleToken_, genParticles);
	if (!genParticles.isValid()) return false;

	const int pdg_tau = 15;
	const int pdg_pi0 = 111;
	const int pdg_pi1 = 211;
	const int pdg_W   = 24;

	const reco::Candidate* tau = nullptr;
	double dRmin = null;
	for (auto& particle: *genParticles) {
		// look for the tau -> pi+ pi0 neutrino decay most oriented towards
		// the reconstructed tau (if present)
#define cut(condition) if (!(condition)) continue;
		cut(abs(particle.pdgId()) == pdg_tau);
		cut(particle.numberOfDaughters() == 3);

		const reco::Candidate* pi0 = nullptr;
		const reco::Candidate* pi1 = nullptr;
		for (int i = 0; i < 3; ++i) {
			const reco::Candidate* daughter = particle.daughter(i);
			int id = abs(daughter->pdgId());
			if (id == pdg_pi0)
				pi0 = daughter;
			else if (id == pdg_pi1)
				pi1 = daughter;
		};
		cut(pi0 && pi1);

		cut(particle.pt() > tauPtMin);
		cut(TMath::Abs(particle.eta()) < tauEtaMax);
		cut((pv_position - particle.vertex()).R() < tauDzMax);
		cut(pi1->pt() > piPtMin);

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
#undef cut
	};
	if (!tau) return false;

	dR           = dRmin;
	gentau_found = 1;
	genTauFromW  = 0;
	for (auto p = tau->mother(); p; p = p->mother()) {
		genTauMother = p->pdgId();
		if (abs(p->pdgId()) == pdg_W) {
			genTauFromW = 1;
			W_pt        = p->pt();
			break;
		};
	};
	return true;
};

bool TreeMakerMiniAOD::CheckMuon(const edm::Event& event) {
	edm::Handle<pat::MuonCollection> muons;
	event.getByToken(MuonCollectionToken_, muons);
	if (!muons.isValid()) return false;
	for (auto& muon: *muons) if (muon.pt() > 15) return true;
	return false;
}

bool TreeMakerMiniAOD::CheckElectron(const edm::Event& event) {
	edm::Handle<pat::ElectronCollection> electrons;
	event.getByToken(ElectronCollectionToken_, electrons);
	if (!electrons.isValid()) return false;
	for (auto& electron: *electrons) if (electron.pt() > 15) return true;
	return false;
}

bool TreeMakerMiniAOD::AddMET(const edm::Event& event) {
	edm::Handle<pat::METCollection> mets;
	event.getByToken(MetCollectionToken_, mets);
	if (!mets.isValid() || !mets->size()) return false;

	auto& MET = mets->front();
	met = MET.pt();
	met_phi = MET.phi();
	met_eta = MET.eta();
	met_significance = MET.significance();
	met_mEtSig       = MET.mEtSig();
	return true;
}

bool TreeMakerMiniAOD::AddVertex(const edm::Event& event) {
	edm::Handle<reco::VertexCollection> vertices;
	event.getByToken(PVToken_, vertices);
	if (!vertices.isValid()) return false;

	nVtx = vertices->size();
	if (nVtx == 0) return false;
	pv_position = vertices->front().position();
	return true;
};

bool TreeMakerMiniAOD::JetPtSum(const edm::Event& event) {
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
	// Puppi Jets
	PuppijetPtSum15 = 0;
	PuppijetPtSum20 = 0;
	nPuppiJets20    = 0;
	PuppijetPtSum15PV = 0;
	PuppijetPtSum20PV = 0;
	nPuppiJets20PV    = 0;
	nLooseBtagedPuppiJets    = 0;
	nMediumBtagedPuppiJets   = 0;
	nTightBtagedPuppiJets    = 0;
	nLooseBtagedPuppiJetsPV  = 0;
	nMediumBtagedPuppiJetsPV = 0;
	nTightBtagedPuppiJetsPV  = 0;

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

    // Jets
	edm::Handle<pat::JetCollection> jets;
	event.getByToken(JetCollectionToken_, jets);

	if (!jets.isValid()) return false;

    /*
	for (auto& jet: *jets) {
		if (TMath::Abs(jet.eta()) > 3) continue;
		if (jet.pt() > 10) jetPtSum10 += jet.pt();
		if (jet.pt() > 15) jetPtSum15 += jet.pt();
		if (jet.pt() > 20) {
			jetPtSum20 += jet.pt();
            ijet++;
			std::cout << "Jet Number " << ijet << std::endl;
			std::cout << "Pt   = " << jet.pt() << std::endl;
			std::cout << "btag = " << jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") << std::endl;
			++nJets20;
		}
	}
	*/
	int i = 0;
	for (auto& jet: *jets) {
		if (TMath::Abs(jet.eta()) > 3) continue;
		if (jet.pt() > 15) {
			i++;
			if (monitoring) std::cout << std::endl << "Jet[" << i << "] Pt = " << jet.pt() << std::endl;
			reco::TrackRefVector JetTracksRef = jet.associatedTracks();
			std::map <int, int> VertexTracksMap;
			std::map <int, std::pair<int, double>> VertexTracksPtMap;
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

    		jetPtSum15 += jet.pt();
			if (JetFromPV) {
				jetPtSum15PV += jet.pt();
				if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_loose) {
    				nLooseBtagedJetsPV++;
    				if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_medium) {
    					nMediumBtagedJetsPV++;
    					if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_tight) {
    						nTightBtagedJetsPV++;
    					}
    				}
    			}
			}
			if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_loose) {
    			nLooseBtagedJets++;
    			if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_medium) {
    				nMediumBtagedJets++;
    				if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_tight) {
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

	edm::Handle<pat::JetCollection> Puppijets;
	event.getByToken(PuppiJetCollectionToken_, Puppijets);
	if (!Puppijets.isValid()) return false;

	/*
	ijet = 0;
	for (auto& jet: *Puppijets) {
		if (TMath::Abs(jet.eta()) > 3) continue;
		if (jet.pt() > 20) {
			PuppijetPtSum20 += jet.pt();
			ijet++;
			std::cout << "Puppi Jet Number " << ijet << std::endl;
			std::cout << "Pt   = " << jet.pt() << std::endl;
			std::cout << "btag = " << jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") << std::endl;
			++nPuppiJets20;
		}
	}
	*/
	int j = 0;
	// Lopp over Jets
	for (auto& jet: *Puppijets) {
		if (TMath::Abs(jet.eta()) > 3) continue;
		if (jet.pt() > 15) {
			j++;
			if (monitoring) {
				std::cout << std::endl;
				std::cout << "PuppiJet[" << j << "] Pt = " << jet.pt() << std::endl;
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

    		if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > BJet_btag) {
    			BJet_pt   = jet.pt();
	            BJet_eta  = jet.eta();
				BJet_phi  = jet.phi();
				BJet_m    = jet.mass();
				BJet_btag = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb");
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
			
			PuppijetPtSum15 += jet.pt();
			if (JetFromPV) {
				PuppijetPtSum15PV += jet.pt();
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

			if (jet.pt() > 20) {
				PuppijetPtSum20 += jet.pt();
				nPuppiJets20++;
				if (JetFromPV) {
					nPuppiJets20PV++;
					PuppijetPtSum20PV += jet.pt();
				}
			}
		}
	}

	return nJets20 > 0;
};

void TreeMakerMiniAOD::CountTracks(const edm::Event& event) {
	nTrks = 0;
	edm::Handle<pat::IsolatedTrackCollection> tracks;
	event.getByToken(TrackToken_, tracks);
	if (!tracks.isValid()) return;
	nTrks = tracks->size();
};

void TreeMakerMiniAOD::AddWT(const edm::Event& iEvent) {

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

void TreeMakerMiniAOD::AddWTData(const edm::Event& iEvent) {

	int One = 1;
	WTisValid = One;
	WT        = One;
	WTFlip    = One;
	WThminus  = One;
	WThplus   = One;
	TauSpinnerMother = null;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TreeMakerMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.add("TreeMakerMiniAOD", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TreeMakerMiniAOD);
