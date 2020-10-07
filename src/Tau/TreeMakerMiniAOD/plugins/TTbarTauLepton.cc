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

	virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	
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
	std::vector<std::string>  trigNamesEmpty;

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
	
	double tau_absIso;
	double tau_againstElectronRaw;
	double tau_IsoMVArun2v1DBdR03oldDMwLTraw;
	double tau_IsoMVArun2v1DBnewDMwLTraw;
	double tau_IsoMVArun2v1DBoldDMwLTraw;
	double tau_IsoMVArun2v1PWdR03oldDMwLTraw;
	double tau_IsoMVArun2v1PWnewDMwLTraw;
	double tau_IsoMVArun2v1PWoldDMwLTraw;
	double tau_Deep2017v2ElectronRejection;
	double tau_Deep2017v2MuonRejection;
	double tau_Deep2017v2JetRejection;

	double tau_VVLooseDeepTau2017v2VSjet;
	double tau_VLooseDeepTau2017v2VSjet;
	double tau_LooseDeepTau2017v2VSjet;
	double tau_MediumDeepTau2017v2VSjet;
	double tau_TightDeepTau2017v2VSjet;
	double tau_VTightDeepTau2017v2VSjet;
	double tau_VVTightDeepTau2017v2VSjet;

	double tau_LooseDeepTau2017v2VSmu;
	double tau_MediumDeepTau2017v2VSmu;
	double tau_TightDeepTau2017v2VSmu;

	double tau_VVLooseDeepTau2017v2VSe;
	double tau_VLooseDeepTau2017v2VSe;
	double tau_LooseDeepTau2017v2VSe;
	double tau_MediumDeepTau2017v2VSe;
	double tau_TightDeepTau2017v2VSe;
	double tau_VTightDeepTau2017v2VSe;
	double tau_VVTightDeepTau2017v2VSe;

	
	double tau_looseCombinedIso;
	double tau_mediumCombinedIso;
	double tau_tightCombinedIso;
	double tau_VVlooseMvaIso;
	double tau_VlooseMvaIso;
	double tau_looseMvaIso;
	double tau_mediumMvaIso;
	double tau_tightMvaIso;
	double tau_VtightMvaIso;
	double tau_VVtightMvaIso;
	double tau_looseMuonRejection;
	double tau_tightMuonRejection;
	double tau_looseElectronRejection;
	double tau_mediumElectronRejection;
	double tau_tightElectronRejection;
	double tau_VtightElectronRejection;
	double decayModeFindingNewDMs;
	double decayModeFinding;

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
	
	/*
	double tau2_pt;
	double tau2_eta;
	double tau2_phi;
	double tau2_dm;
	double tau2_dz;
	double tau2_q;
	double tau2_m;
	double tau2_absIso;
	*/

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

	double PuppijetPtSum15, PuppijetPtSum20, PuppijetPtSum15PV, PuppijetPtSum20PV;
	int nPuppiJets20, nPuppiJets20PV;
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
	double lepton2_tauAbsIso;
	double m_ll;
	int    nLeptonCandidates;

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
	bool useHLT;
	bool TauSpinnerOn;
	double tauPtMin;
	double piPtMin;
	double tauEtaMax;
	double tauDzMax;
	double METcut;
	double BJetPtMin;
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
TTbarTauLepton::TTbarTauLepton(const edm::ParameterSet& iConfig) {
	//now do what ever initialization is needed
	iT =0;

	isMC						= iConfig.getParameter<bool>("isMC");
	useHLT                      = iConfig.getParameter<bool>("useHLT");
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
	MuElePtMin                  = iConfig.getParameter<double>("MuElePtMin"); // 20
	EtaMax                      = iConfig.getParameter<double>("EtaMax"); // 2.4
	DeepTau  					= iConfig.getParameter<bool>("DeepTau");
	null 						= iConfig.getParameter<double>("null");

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
	trigNamesEmpty                    = iConfig.getParameter<std::vector<std::string>>("Triggers");
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

	// Operation with triggers (filtering, scaling)
	if (!TriggerOK(event)) {
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
	if (!AddTau(event)) {
		if (monitoring) std::cout << "Tau" << std::endl;
		return;
	}
	//AddPackedCandidates(event);
	//
	CountTracks(event);
	// Summary of jets parameters
	if (!JetPtSum(event)) {
		if (monitoring) std::cout << "Jet" << std::endl;
		return;
	}
	// Find second lepton for ttbar
	if (!AddLepton(event)) {
		if (monitoring) std::cout << "Lepton" << std::endl;
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

	TLorentzVector pi0, pi1;
	pi0.SetPtEtaPhiM(piZero_pt, piZero_eta, piZero_phi, piZero_m);
	pi1.SetPtEtaPhiM(piChar_pt, piChar_eta, piChar_phi, piChar_m);

	pipiMass = (pi0 + pi1).M();
//	upsilon  = (pi1.E() - pi0.E()) / tau_pt;
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
	tree->Branch("tau_looseCombinedIso",&tau_looseCombinedIso,"tau_looseCombinedIso/D");
	tree->Branch("tau_mediumCombinedIso",&tau_mediumCombinedIso,"tau_mediumCombinedIso/D");
	tree->Branch("tau_tightCombinedIso",&tau_tightCombinedIso,"tau_tightCombinedIso/D");
	tree->Branch("tau_VVlooseMvaIso",&tau_VVlooseMvaIso,"tau_VVlooseMvaIso/D");
	tree->Branch("tau_VlooseMvaIso",&tau_VlooseMvaIso,"tau_VlooseMvaIso/D");
	tree->Branch("tau_looseMvaIso",&tau_looseMvaIso,"tau_looseMvaIso/D");
	tree->Branch("tau_mediumMvaIso",&tau_mediumMvaIso,"tau_mediumMvaIso/D");
	tree->Branch("tau_tightMvaIso",&tau_tightMvaIso,"tau_tightMvaIso/D");
	tree->Branch("tau_VtightMvaIso",&tau_VtightMvaIso,"tau_VtightMvaIso/D");
	tree->Branch("tau_VVtightMvaIso",&tau_VVtightMvaIso,"tau_VVtightMvaIso/D");
	tree->Branch("tau_looseMuonRejection",&tau_looseMuonRejection,"tau_looseMuonRejection/D");
	tree->Branch("tau_tightMuonRejection",&tau_tightMuonRejection,"tau_tightMuonRejection/D");
	tree->Branch("tau_looseElectronRejection",&tau_looseElectronRejection,"tau_looseElectronRejection/D");
	tree->Branch("tau_mediumElectronRejection",&tau_mediumElectronRejection,"tau_mediumElectronRejection/D");
	tree->Branch("tau_tightElectronRejection",&tau_tightElectronRejection,"tau_tightElectronRejection/D");
	tree->Branch("tau_VtightElectronRejection",&tau_VtightElectronRejection,"tau_VtightElectronRejection/D");

	// Deep 2017v2
	tree->Branch("tau_Deep2017v2ElectronRejection",&tau_Deep2017v2ElectronRejection,"tau_Deep2017v2ElectronRejection/D");
	tree->Branch("tau_Deep2017v2MuonRejection",&tau_Deep2017v2MuonRejection,"tau_Deep2017v2MuonRejection/D");
	tree->Branch("tau_Deep2017v2JetRejection",&tau_Deep2017v2JetRejection,"tau_Deep2017v2JetRejection/D");

	tree->Branch("tau_VVLooseDeepTau2017v2VSjet",&tau_VVLooseDeepTau2017v2VSjet,"tau_VVLooseDeepTau2017v2VSjet/D");
	tree->Branch("tau_VLooseDeepTau2017v2VSjet",&tau_VLooseDeepTau2017v2VSjet,"tau_VLooseDeepTau2017v2VSjet/D");
	tree->Branch("tau_LooseDeepTau2017v2VSjet",&tau_LooseDeepTau2017v2VSjet,"tau_LooseDeepTau2017v2VSjet/D");
	tree->Branch("tau_MediumDeepTau2017v2VSjet",&tau_MediumDeepTau2017v2VSjet,"tau_MediumDeepTau2017v2VSjet/D");
	tree->Branch("tau_TightDeepTau2017v2VSjet",&tau_TightDeepTau2017v2VSjet,"tau_TightDeepTau2017v2VSjet/D");
	tree->Branch("tau_VTightDeepTau2017v2VSjet",&tau_VTightDeepTau2017v2VSjet,"tau_VTightDeepTau2017v2VSjet/D");
	tree->Branch("tau_VVTightDeepTau2017v2VSjet",&tau_VVTightDeepTau2017v2VSjet,"tau_VVTightDeepTau2017v2VSjet/D");
	tree->Branch("tau_LooseDeepTau2017v2VSmu",&tau_LooseDeepTau2017v2VSmu,"tau_LooseDeepTau2017v2VSmu/D");
	tree->Branch("tau_MediumDeepTau2017v2VSmu",&tau_MediumDeepTau2017v2VSmu,"tau_MediumDeepTau2017v2VSmu/D");
	tree->Branch("tau_TightDeepTau2017v2VSmu",&tau_TightDeepTau2017v2VSmu,"tau_TightDeepTau2017v2VSmu/D");

	tree->Branch("tau_VVLooseDeepTau2017v2VSe",&tau_VVLooseDeepTau2017v2VSe,"tau_VVLooseDeepTau2017v2VSe/D");
	tree->Branch("tau_VLooseDeepTau2017v2VSe",&tau_VLooseDeepTau2017v2VSe,"tau_VLooseDeepTau2017v2VSe/D");
	tree->Branch("tau_LooseDeepTau2017v2VSe",&tau_LooseDeepTau2017v2VSe,"tau_LooseDeepTau2017v2VSe/D");
	tree->Branch("tau_MediumDeepTau2017v2VSe",&tau_MediumDeepTau2017v2VSe,"tau_MediumDeepTau2017v2VSe/D");
	tree->Branch("tau_TightDeepTau2017v2VSe",&tau_TightDeepTau2017v2VSe,"tau_TightDeepTau2017v2VSe/D");
	tree->Branch("tau_VTightDeepTau2017v2VSe",&tau_VTightDeepTau2017v2VSe,"tau_VTightDeepTau2017v2VSe/D");
	tree->Branch("tau_VVTightDeepTau2017v2VSe",&tau_VVTightDeepTau2017v2VSe,"tau_VVTightDeepTau2017v2VSe/D");

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
	tree->Branch("lepton2_tauAbsIso", &lepton2_tauAbsIso, "lepton2_tauAbsIso/D");
	tree->Branch("m_ll", &m_ll, "m_ll/D");
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
	//tree->Branch("m_BJet1BJet2", &m_BJet1BJet2, "m_BJet1BJet2/D");

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
	//tree->Branch("WTisValid", &WTisValid, "WTisValid/D")
	// add more branches
	
	allTauPt = FS->make<TH1F>("allTauPt","allTauPt",300,0,300);
}

// ------------ method called once each job just after ending the event loop  ------------
void TTbarTauLepton::endJob() {

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

	bool triggerOK = false;
	nTauTriggers = 0;
	nJetHTTriggers = 0;
	nMETTriggers = 0;
	nBTagCSVTriggers = 0;
	nBTagMuTriggers = 0;
	nSingleMuonTriggers = 0;
	nSingleElectronTriggers = 0;
    /////////////////////////////TriggerResults////////////////////////////////////
	edm::Handle<edm::TriggerResults> triggerResults;
	iEvent.getByToken(tok_trigRes, triggerResults);
	if (triggerResults.isValid()) {
		const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
		const std::vector<std::string> & triggerNames_ = triggerNames.triggerNames();
		for ( unsigned int iHLT=0; iHLT<triggerResults->size(); iHLT++ ) {
			int hlt    = triggerResults->accept(iHLT);
			if ( hlt > 0 ) {
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
						if (monitoringHLT) std::cout << "SingleMuon Trigger" << std::endl;
					}
				}
				for ( unsigned int i=0; i<trigNamesSingleElectron.size(); ++i ) {
					if ( triggerNames_[iHLT].find(trigNamesSingleElectron[i].c_str())!= std::string::npos ) {
						nSingleElectronTriggers++;
						if (monitoringHLT) std::cout << "SingleElectron Trigger" << std::endl;
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
	if (nTauTriggers + nJetHTTriggers + nMETTriggers + nBTagCSVTriggers + nSingleMuonTriggers + nSingleElectronTriggers > 0) triggerOK = true;
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
	for (size_t i = 0; i < taus->size(); ++i) {
#define cut(condition) if (!(condition)) continue;
		auto& tau = (*taus)[i];
		if (monitoringTau) std::cout << "######### Tau candidate number " << i << " ################" << std::endl;
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
	tau_VtightElectronRejection  = tau.tauID("againstElectronVTightMVA6");
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

		cut(particle.pt() > tauPtMin);
		//if (monitoringGen) std::cout << "gentau pt" << std::endl;
		cut(TMath::Abs(particle.eta()) < tauEtaMax);
		//if (monitoringGen) std::cout << "gentau eta" << std::endl;
		cut((pv_position - particle.vertex()).R() < tauDzMax);
		//if (monitoringGen) std::cout << "gentau dz" << std::endl;
		//cut(pi1->pt() > piPtMin);

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
		};
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
		if (abs(particle.pdgId()) == pdg_mu || abs(particle.pdgId()) == pdg_electron || abs(particle.pdgId()) == pdg_tau) {
			double lepdR1 = sqrt(deltaR2(particle.eta(), particle.phi(), lepton1_eta, lepton1_phi));
			double lepdR2 = sqrt(deltaR2(particle.eta(), particle.phi(), lepton2_eta, lepton2_phi));
			double dRtau  = sqrt(deltaR2(particle.eta(), particle.phi(), gentau_eta, gentau_phi));
			double leptondR_ = TMath::Min(lepdR1, lepdR2); 
			if ( (!lepton || leptondR_ < leptondRmin) && dRtau > 0.4) {
				lepton = &particle;
				leptondRmin = leptondR_;
				if (monitoringGen) {
					std::cout << "---------- Lepton" << std::endl;
					GenRecoMonitor *Lepton1Gen = new GenRecoMonitor(particle, lepton1_flavor, lepton1_pt, lepton1_eta, lepton1_phi);
					Lepton1Gen->PrintComp(true, true);
					delete Lepton1Gen;
					if (lepton2_pt > 0) {
						GenRecoMonitor *Lepton2Gen = new GenRecoMonitor(particle, lepton2_flavor, lepton2_pt, lepton2_eta, lepton2_phi);
						Lepton2Gen->PrintComp(true, true);
						delete Lepton2Gen;
					}
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

	gentau_vis_pt = TauVisP4.pt();
	gentau_vis_eta = TauVisP4.eta();
	gentau_vis_phi = TauVisP4.phi();
	gentau_vis_energy = TauVisP4.energy();

	SumNu_pt = SumNuP4.pt();
	SumNu_eta = SumNuP4.eta();
	SumNu_phi = SumNuP4.phi();
	SumNu_energy = SumNuP4.energy();

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

/*

void TTbarTauLepton::AddGenMET(const edm::Event& event) {

	edm::Handle<GenMETCollection> genMetTrue;
	event.getByToken(tok_GenMetTrue_, genMetTrue);

	if(!isMC) {
		TrueMetPt = null;
		TrueMetEta = null;
		TrueMetPhi = null;
		TrueMetEnergy = null;
		return;
	}

    if (genMetTrue.isValid() && !(*genMetTrue).empty()) {
        if (monitoring) std::cout << "True MET vector size = " << genMetTrue->size() << std::endl;
        auto& METTrue = genMetTrue->front();
        if (monitoringMET) {
            std::cout << "True MET:" << std::endl;
            std::cout << "Pt      = " << METTrue.pt() << std::endl
            << "Eta     = " << METTrue.eta() << std::endl
            << "Phi     = " << METTrue.phi() << std::endl
            << "Energy  = " << METTrue.energy() << std::endl;
        }
        TrueMetPt = METTrue.pt();
        TrueMetEta = METTrue.eta();
        TrueMetPhi = METTrue.phi();
        TrueMetEnergy = METTrue.energy();
    }

}
*/

/*
void TTbarTauLepton::AddPackedCandidates(const edm::Event& event) {

	edm::Handle<pat::PackedCandidateCollection> PackedCandidates;
	event.getByToken(PackedCandidateCollectionToken_, PackedCandidates);
	if (!PackedCandidates.isValid() || PackedCandidates->size() == 0) {
		std::cout << "Packed candidates are not avalible" << std::endl;
	}

	std::cout << "Packed candidates close to tau lepton:" << std::endl;
	int nCand = 0;
	for (auto& PackedCandidate: *PackedCandidates) {
		double deltaRTau = sqrt(deltaR2(tau_eta, tau_phi, PackedCandidate.eta(), PackedCandidate.phi()));
		if (deltaRTau < 0.4) std::cout << "Candidate [" << nCand << "]: pdgId = " << PackedCandidate.pdgId() << "   Pt = " << PackedCandidate.pt() << std::endl;
		nCand++;
	}

};
*/

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

	/*

	edm::Handle<reco::VertexCompositePtrCandidateCollection> SecondaryVertices;
	event.getByToken(SVToken_, SecondaryVertices);

	if (SecondaryVertices.isValid() && SecondaryVertices->size() > 0) {
		std::cout << "Number of secondary vertices = " << SecondaryVertices->size() << std::endl;
		for (unsigned nSV = 0; nSV < SecondaryVertices->size(); nSV++) {
			std::cout << "SV [" << nSV << "] position = " << "(" << (*vertices)[nSV].position().x() << ", " << (*vertices)[nSV].position().y() << ", " << (*vertices)[nSV].position().z() << ")" << std::endl; 
		}
	} else {
		std::cout << "SV collection is not valid or number of SVs = " << SecondaryVertices->size() << std::endl;
	}
	*/

	nVtx = vertices->size();
	if (nVtx == 0) return false;
	pv_position = vertices->front().position();
	return true;
};

// Jets analysis

bool TTbarTauLepton::JetPtSum (const edm::Event& event) {

	if (monitoringBJets) std::cout << std::endl << "!!! Jets !!!" << std::endl;

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
		if (TMath::Abs(jet.eta()) > 3) continue;
		if (jet.pt() > 15) {
			j++;
			if (monitoringJets) {
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
			if (jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb") > WPBTag_loose) {
				// BJetFound = true;
				if (TMath::Abs(jet.eta()) > EtaMax) continue;
				if (jet.pt() < BJetPtMin) continue;
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

	if (monitoringBJets) std::cout << std::endl << "!!! bTagged Jets !!!" << std::endl;

	SortJets(looseBJets);

	if (monitoringBJets) std::cout << "List of loose bJets (" << looseBJets.size() << ", " << nLooseBJets << "):" << std::endl;
	for(unsigned ijet = 0; ijet < looseBJets.size(); ijet++) {
		if (monitoringBJets) {
			std::cout << "jet " << ijet << std::endl;
			looseBJets[ijet].Monitoring();
		}
	}

	if (looseBJets.size() < 2) return false;

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

	looseBJets.clear();

	return nPuppiJets20 > 0;
};

bool TTbarTauLepton::AddLepton (const edm::Event& event) {

	if (monitoringLeptons) std::cout << std::endl << "!!! Leptons !!!:" << std::endl;

	nLeptonCandidates = 0;

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
	lepton2_tauAbsIso   = null;
	m_ll                = null;

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
	if (electrons.isValid() && !(*electrons).empty()) {
		//
		if (monitoringLeptons) std::cout << "Isolated electrons:" << std::endl;
		for (auto& electron: *electrons) {
			nLeptons++;
			if (monitoringLeptons) std::cout << "Electron " << nLeptons << std::endl;
#define cut(condition) if (!(condition)) continue;
			cut(electron.pt() > MuElePtMin);
			cut(abs(electron.eta()) < EtaMax);
			cut((pv_position - electron.vertex()).R() < tauDzMax);
			if (monitoringLeptons) std::cout << "Passed cuts" << std::endl;
			double deltaRJet1 = sqrt(deltaR2(BJet1_eta, BJet1_phi, electron.eta(), electron.phi()));
			double deltaRJet2 = sqrt(deltaR2(BJet2_eta, BJet2_phi, electron.eta(), electron.phi()));
			double deltaRTau  = sqrt(deltaR2(tau_eta, tau_phi, electron.eta(), electron.phi()));

			if (deltaRTau < 0.3) {
				if (monitoringLeptons) std::cout << "In Tau cone: dR = " << deltaRTau << std::endl;
				return false;
			}
			if (deltaRJet1 > 0.4 && deltaRJet2 > 0.4) {
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
				if (monitoringLeptons) std::cout << "Added to vector" << std::endl;
	//delete Lepton;
			} else if (monitoringLeptons) {
				std::cout << "Close to b-jets: dR1 = " << deltaRJet1 << ", dR2 = " << deltaRJet2 << std::endl;
			}
#undef cut
		}
	}

	if (muons.isValid() && !(*muons).empty()) {
		if (monitoringLeptons) std::cout << "Isolated muons:" << std::endl;
		for (auto& muon: *muons) {
			nLeptons++;
			if (monitoringLeptons) std::cout << "Muon " << nLeptons << std::endl;
#define cut(condition) if (!(condition)) continue;
			cut(muon.pt() > MuElePtMin);
			cut(abs(muon.eta()) < EtaMax);
			cut((pv_position - muon.vertex()).R() < tauDzMax);
			if (monitoringLeptons) std::cout << "Passed cuts" << std::endl;
			double deltaRJet1 = sqrt(deltaR2(BJet1_eta, BJet1_phi, muon.eta(), muon.phi()));
			double deltaRJet2 = sqrt(deltaR2(BJet2_eta, BJet2_phi, muon.eta(), muon.phi()));
			double deltaRTau  = sqrt(deltaR2(tau_eta, tau_phi, muon.eta(), muon.phi()));
			//double dphiJet1 = dphi(Jet1_phi, muon.phi());
			//double dphiJet2 = dphi(Jet2_phi, muon.phi());
			if (deltaRTau < 0.3) {
				if (monitoringLeptons) std::cout << "In Tau cone: dR = " << deltaRTau << std::endl;
				return false;
			}
			if (deltaRJet1 > 0.4 && deltaRJet2 > 0.4) {
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
				if (monitoringLeptons) std::cout << "Added to vector" << std::endl;
	//delete Lepton;
			} else if (monitoringLeptons) {
				std::cout << "Close to b-jets: dR1 = " << deltaRJet1 << ", dR2 = " << deltaRJet2 << std::endl;
			}
#undef cut
		}
	}

	// tau as second lepton candidate
	if (taus.isValid() && !(*taus).empty()) {
		//
		if (monitoringLeptons) std::cout << "Isolated taus:" << std::endl;
		for (auto& tau: *taus) {
			nLeptons++;
			if (monitoringLeptons) std::cout << "Tau " << nLeptons << std::endl;
#define cut(condition) if (!(condition)) continue;
			cut(tau.pt() > MuElePtMin);
			cut(abs(tau.eta()) < EtaMax);
			cut((pv_position - tau.vertex()).R() < tauDzMax);
			if (monitoringLeptons) std::cout << "Passed cuts" << std::endl;
			double deltaRJet1 = sqrt(deltaR2(BJet1_eta, BJet1_phi, tau.eta(), tau.phi()));
			double deltaRJet2 = sqrt(deltaR2(BJet2_eta, BJet2_phi, tau.eta(), tau.phi()));
			double deltaRTau  = sqrt(deltaR2(tau_eta, tau_phi, tau.eta(), tau.phi()));
			cut(deltaRTau > 0.3);
			if (deltaRJet1 > 0.4 && deltaRJet2 > 0.4) {
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
				if (monitoringLeptons) std::cout << "Added to vector" << std::endl;
	//delete Lepton;
			} else if (monitoringLeptons) {
				std::cout << "Close to b-jets: dR1 = " << deltaRJet1 << ", dR2 = " << deltaRJet2 << std::endl;
			}
#undef cut
		}
	}

	if (nLeptonCandidates < 1) {
		if (monitoringLeptons) std::cout << "Number of lepton candiates = " << LepCandidates.size() << ", " << nLeptonCandidates << std::endl;
		return false;
	}

	SortLeptons(LepCandidates);

	if (monitoringLeptons) std::cout << "List of selected leptons:" << std::endl;
	for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
		if (monitoringLeptons) {
			std::cout << "lepton " << ilep << std::endl;
			LepCandidates[ilep].Monitoring();
			std::cout << "delta(BJet1) = " << sqrt(deltaR2(BJet1_eta, BJet1_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
			std::cout << "delta(BJet2) = " << sqrt(deltaR2(BJet2_eta, BJet2_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
			std::cout << "delta(Tau)   = " << sqrt(deltaR2(tau_eta, tau_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
		}
	}

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
		lepton2_tauAbsIso = LepCandidates[1].tauAbsIso;
		m_ll             = (LepCandidates[0].FourMomentum + LepCandidates[1].FourMomentum).M();
	}

	LepCandidates.clear();

	return (nLeptonCandidates > 0);
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