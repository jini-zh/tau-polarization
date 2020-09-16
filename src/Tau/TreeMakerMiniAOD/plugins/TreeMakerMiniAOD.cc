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
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"

#include "Tau/TreeMakerMiniAOD/plugins/ParticleMonitor.h"
#include "Tau/TreeMakerMiniAOD/plugins/BJetCandidate.h"
#include "Tau/TreeMakerMiniAOD/plugins/LeptonCandidate.h"
#include "Tau/TauAnalyzer/plugins/MySimpleParticle.h"


static inline double sqr(double x) {
	return x * x;
}

static double dphi(double phi1, double phi2) {
	double result = TMath::Abs(phi1 - phi2);
	if (result < TMath::Pi()) return result;
	return 2 * TMath::Pi() - result;
}

/*
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Definition of class BJetCandidate
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
class BJetCandidate {
public :
    //--- variables
    double Pt;
    double E;
    double Eta, Phi;
    math::XYZTLorentzVector FourMomentum;
    double Bprob;
    bool FromPV;

    //--- constructor & destructor
    BJetCandidate();
    virtual ~BJetCandidate();

    //--- functions

};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Constructor of BJetCandidate object
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
BJetCandidate::BJetCandidate() {

}

BJetCandidate::~BJetCandidate()
{
// Destructor
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Definition of class LeptonCandidate
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
class LeptonCandidate {
public :
    //--- variables
    double Pt;
    math::XYZTLorentzVector FourMomentum;
    double Eta, Phi;
    int Flavor;
    int Charge;
    double Dz;
    // Methods from pat::Muon
    float caloIso;
    float ecalIso;
    float hcalIso;
    float puppiChargedHadronIso;
    float puppiNeutralHadronIso;
    float puppiNoLeptonsChargedHadronIso;
    float puppiNoLeptonsNeutralHadronIso;
    float puppiNoLeptonsPhotonIso;
    float puppiPhotonIso;
    float trackIso;

    //--- constructor & destructor
    LeptonCandidate();
    virtual ~LeptonCandidate();

    //--- functions

};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Constructor of LeptonCandidate object
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
LeptonCandidate::LeptonCandidate() {

}

LeptonCandidate::~LeptonCandidate()
{
// Destructor
}

// Jets in order of decreasing btag

void SortJets(std::vector<BJetCandidate>& items) {
  bool swapped;
  do {
    swapped = false;
    for (unsigned i = 1; i < items.size(); i++) {
      if (items[i-1].Bprob < items[i].Bprob) {
        std::swap(items[i-1], items[i]);
        swapped = true;
      }
    }
  } while (swapped != false);
}

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
*/

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

class TreeMakerMiniAOD : public edm::EDAnalyzer {
public:
	explicit TreeMakerMiniAOD(const edm::ParameterSet&);
	~TreeMakerMiniAOD();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
	virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	virtual void endJob() override;
	
	bool AddTau            (const edm::Event&);
	bool FindGenTau        (const edm::Event&);
	bool AddLepton         (const edm::Event&);
	bool AddMET            (const edm::Event&);
	bool AddVertex         (const edm::Event&);
	void CountTracks       (const edm::Event&);
	bool JetPtSum          (const edm::Event&);
	void AddWT             (const edm::Event&);
	void AddWTData         (const edm::Event&);
	bool DecaychannelMatch (std::vector<MySimpleParticle> &particles, int p1, int p2=0, int p3=0, int p4=0, int p5=0, int p6=0);
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
	edm::EDGetTokenT<pat::JetCollection> JetCollectionToken_;
	edm::EDGetTokenT<pat::JetCollection> PuppiJetCollectionToken_;
	edm::EDGetTokenT<pat::METCollection> MetCollectionToken_;
	edm::EDGetTokenT<reco::VertexCollection> PVToken_;
	//edm::EDGetTokenT<pat::PackedCandidateCollection> GenParticleToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> GenParticleToken_;
	
	edm::EDGetTokenT<pat::IsolatedTrackCollection> TrackToken_;
	edm::EDGetTokenT<GenEventInfoProduct> GenEventInfoToken_;
	edm::EDGetTokenT<GenRunInfoProduct> GenRunInfoToken_;
	// HLT
	edm::InputTag theTriggerResultsLabel;
	edm::EDGetTokenT<edm::TriggerResults> tok_trigRes;
	std::vector<std::string>  trigNamesTau;
	std::vector<std::string>  trigNamesJetHT;
	std::vector<std::string>  trigNamesMET;
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
	
	double tau_absIso;
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
	
	double tau2_pt;
	double tau2_eta;
	double tau2_phi;
	double tau2_dm;
	double tau2_dz;
	double tau2_q;
	double tau2_m;
	double tau2_absIso;

	double tau_found;
	double gentau_found;
	double gentau_dm, gentau_nPi0;
	double dR;
	double bquark1dR;
	double bquark2dR;
	double leptondR;
	double genTauFromW;
	int genTauMother;
	double W_pt;

	// Generated tau parameters
	double gentau_pt, gentau_px, gentau_py, gentau_pz, gentau_eta, gentau_phi;
	double genPiChar_pt, genPiChar_px, genPiChar_py, genPiChar_pz, genPiChar_eta, genPiChar_phi;
	double genPi0_pt, genPi0_px, genPi0_py, genPi0_pz, genPi0_eta, genPi0_phi;
	double nuW_pt, nuW_px, nuW_py, nuW_pz, nuW_eta, nuW_phi;
	double nutau_pt, nutau_px, nutau_py, nutau_pz, nutau_eta, nutau_phi;
	double nunu_pt;

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
	double lepton1_dz;
	double lepton1_flavor;
	double lepton1_charge;
	double lepton1_E;
	double lepton1_trackIso;

	double lepton2_pt;
	double lepton2_eta;
	double lepton2_phi;
	double lepton2_dz;
	double lepton2_flavor;
	double lepton2_charge;
	double lepton2_E;
	double lepton2_trackIso;
	double m_ll;
	int nLeptonCandidates;

	int nPi0;
	int nPhotons;
	int nGamma;
	double met;
	double met_phi;
	double met_eta;
	double met_energy;
	double tauMET_mass;
	double met_significance;
	double met_mEtSig;
	
	double m_t;
	double dPhi;

	// Trigger matching
	int nTauTriggers;
	int nJetHTTriggers;
	int nMETTriggers;
	
	math::XYZPoint pv_position;

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
TreeMakerMiniAOD::TreeMakerMiniAOD(const edm::ParameterSet& iConfig) {
	//now do what ever initialization is needed
	iT =0;

	isMC						= iConfig.getParameter<bool>("isMC");
	monitoring					= iConfig.getParameter<bool>("monitoring");
	monitoringHLT				= iConfig.getParameter<bool>("monitoringHLT");
	monitoringTau				= iConfig.getParameter<bool>("monitoringTau");
	monitoringGen				= iConfig.getParameter<bool>("monitoringGen");
	monitoringJets				= iConfig.getParameter<bool>("monitoringJets");
	monitoringBJets				= iConfig.getParameter<bool>("monitoringBJets");
	monitoringLeptons			= iConfig.getParameter<bool>("monitoringLeptons");
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
	std::string jetCollection         = iConfig.getParameter<std::string>("jetCollection");
	std::string PuppijetCollection    = iConfig.getParameter<std::string>("PuppijetCollection");
	std::string metCollection         = iConfig.getParameter<std::string>("metCollection");
	std::string vertexCollection      = iConfig.getParameter<std::string>("vertexCollection");
	std::string genParticleCollection = iConfig.getParameter<std::string>("genParticleCollection");
	std::string trackCollection       = iConfig.getParameter<std::string>("trackCollection");
	theTriggerResultsLabel		      = edm::InputTag("TriggerResults","","HLT");
	trigNamesTau                      = iConfig.getParameter<std::vector<std::string>>("TauTriggers");
	trigNamesJetHT                    = iConfig.getParameter<std::vector<std::string>>("JetHTTriggers");
	trigNamesMET                      = iConfig.getParameter<std::vector<std::string>>("METTriggers");
	trigNamesEmpty                    = iConfig.getParameter<std::vector<std::string>>("Triggers");
	
	TauCollectionToken_ 		= consumes<pat::TauCollection>(edm::InputTag(tauCollection));
	MuonCollectionToken_ 		= consumes<pat::MuonCollection>(edm::InputTag(muonCollection));
	ElectronCollectionToken_	= consumes<pat::ElectronCollection>(edm::InputTag(electronCollection));
	JetCollectionToken_ 		= consumes<pat::JetCollection>(edm::InputTag(jetCollection));
	PuppiJetCollectionToken_    = consumes<pat::JetCollection>(edm::InputTag(PuppijetCollection));
	MetCollectionToken_ 		= consumes<pat::METCollection>(edm::InputTag(metCollection));
	PVToken_ 			        = consumes<reco::VertexCollection>(edm::InputTag(vertexCollection));
	//GenParticleToken_ 		= consumes<pat::PackedCandidateCollection>(edm::InputTag(genParticleCollection));
	GenParticleToken_           = consumes<reco::GenParticleCollection>(edm::InputTag(genParticleCollection));
	TrackToken_			        = consumes<pat::IsolatedTrackCollection>(edm::InputTag(trackCollection));
	tok_trigRes                 = consumes<edm::TriggerResults>(theTriggerResultsLabel);

	if (isMC && TauSpinnerOn) {
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
	GenTauDM(event);
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
	tree->Branch("tau_absIso",&tau_absIso,"tau_absIso/D");

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

	// Second tau parameters
	tree->Branch("tau2_pt",&tau2_pt,"tau2_pt/D");
	tree->Branch("tau2_eta",&tau2_eta,"tau2_eta/D");
	tree->Branch("tau2_phi",&tau2_phi,"tau2_phi/D");
	tree->Branch("tau2_q",&tau2_q,"tau2_q/D");
	tree->Branch("tau2_m",&tau2_m,"tau2_m/D");
	tree->Branch("tau2_dm",&tau2_dm,"tau2_dm/D");
	tree->Branch("tau2_dz",&tau2_dz,"tau2_dz/D");
	tree->Branch("tau2_absIso",&tau2_absIso,"tau2_absIso/D");

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

	tree->Branch("pipiMass", &pipiMass, "pipiMass/D");
	tree->Branch("upsilon", &upsilon, "upsilon/D");

	// Generated particles
	tree->Branch("tau_found",&tau_found, "tau_found/D");
	tree->Branch("gentau_found",&gentau_found, "gentau_found/D");
	tree->Branch("dR", &dR, "dR/D");
	tree->Branch("leptondR", &leptondR, "leptondR/D");
	tree->Branch("bquark1dR", &bquark1dR, "bquark1dR/D");
	tree->Branch("bquark2dR", &bquark2dR, "bquark2dR/D");
	tree->Branch("genTauFromW", &genTauFromW, "genTauFromW/D");
	tree->Branch("W_pt", &W_pt, "W_pt/D");
	tree->Branch("genTauMother", &genTauMother, "genTauMother/I");
	tree->Branch("gentau_dm",&gentau_dm,"gentau_dm/D");
	tree->Branch("gentau_nPi0",&gentau_nPi0,"gentau_nPi0/D");
	tree->Branch("gentau_pt",&gentau_pt,"gentau_pt/D");
	tree->Branch("genPiChar_pt",&genPiChar_pt,"genPiChar_pt/D");
	tree->Branch("genPi0_pt",&genPi0_pt,"genPi0_pt/D");
	tree->Branch("nutau_pt",&nutau_pt,"nutau_pt/D");
	tree->Branch("nuW_pt",&nuW_pt,"nuW_pt/D");
	tree->Branch("nunu_pt",&nunu_pt,"nunu_pt/D");

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
	tree->Branch("lepton1_trackIso", &lepton1_trackIso, "lepton1_trackIso/D");
	tree->Branch("lepton2_pt", &lepton2_pt, "lepton2_pt/D");
	tree->Branch("lepton2_E", &lepton2_E, "lepton2_E/D");
	tree->Branch("lepton2_eta", &lepton2_eta, "lepton2_eta/D");
	tree->Branch("lepton2_phi", &lepton2_phi, "lepton2_phi/D");
	tree->Branch("lepton2_dz", &lepton2_dz, "lepton2_dz/D");
	tree->Branch("lepton2_flavor", &lepton2_flavor, "lepton2_flavor/I");
	tree->Branch("lepton2_charge", &lepton2_charge, "lepton2_charge/I");
	tree->Branch("lepton2_trackIso", &lepton2_trackIso, "lepton2_trackIso/D");
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
	
	// Vertices, tracks, tau candidates
	tree->Branch("nVtx",&nVtx,"nVtx/I");
	tree->Branch("nTrks",&nTrks,"nTrks/I");
	tree->Branch("nTau",&nTau,"nTau/I");
	tree->Branch("nTauC",&nTauC,"nTauC/I");

	tree->Branch("nPi0",&nPi0,"nPi0/I");
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

void
TreeMakerMiniAOD::endRun(edm::Run const& Run, edm::EventSetup const&)
{

}

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

//---------------------------------TRIGGER-------------------------------------------------
bool TreeMakerMiniAOD::TriggerOK (const edm::Event& iEvent) {

	if (monitoringHLT) std::cout << std::endl << "!!! Triggers !!!" << std::endl;

	nTauTriggers = 0;
	nJetHTTriggers = 0;
	nMETTriggers = 0;
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
	return true;
}

//-----------------------------------------------------------------------------------------

bool TreeMakerMiniAOD::AddTau(const edm::Event& event) {

	if (monitoringTau) std::cout << std::endl << "!!! Tau !!!" << std::endl;

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
	int index_2 = null;
	std::map <int, std::pair <double, double>> TauIdMap;
	for (size_t i = 0; i < taus->size(); ++i) {
#define cut(condition) if (!(condition)) continue;
		auto& tau = (*taus)[i];
		std::cout << "######### Tau candidate number " << i << " ################" << std::endl;
		//if (monitoringTau) TauMonitor (tau, DeepTau);
		cut(tau.pt() > tauPtMin); // tau transverse momentum
		std::cout << "Pt OK" << std::endl;
		cut(TMath::Abs(tau.eta()) < tauEtaMax); // tau pseudorapidity
		std::cout << "Eta OK" << std::endl;
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
		std::cout << "Isolation OK" << std::endl;

		const reco::CandidatePtr leadPiCh = tau.leadChargedHadrCand();
		//reco::CandidatePtrVector VectorPiCh = tau.signalChargedHadrCands();
		reco::CandidatePtrVector VectorSignalCands = tau.signalCands(); 

		if (monitoringTau) TauMonitor (tau, DeepTau, pv_position);

		reco::CandidatePtrVector VectorSignalCands_Photons;
		VectorSignalCands_Photons.clear();

		if (VectorSignalCands.size() > 0) {
			for (unsigned l = 0; l < VectorSignalCands.size(); l++) {
				if (VectorSignalCands[l]->pdgId() == 22) {
					VectorSignalCands_Photons.push_back(VectorSignalCands[l]);
				}
			}
		}

		math::XYZTLorentzVector Photons_p4;
		double InvaMass_temp = null;

		if (VectorSignalCands_Photons.size() > 0) {
			for (unsigned l = 0; l < VectorSignalCands_Photons.size(); l++) {
				Photons_p4 += VectorSignalCands_Photons[l]->p4();
				for (unsigned n = l; n < VectorSignalCands_Photons.size(); n++) {
					if (n == l) continue;
					double InvaMass_pi0 = (VectorSignalCands_Photons[l]->p4() + VectorSignalCands_Photons[n]->p4()).M();
					double Pt_pi0       = (VectorSignalCands_Photons[l]->p4() + VectorSignalCands_Photons[n]->p4()).Pt();
					if (monitoringTau) std::cout << "mass Photon [" << l <<"][" << n << "] = " << InvaMass_pi0 << std::endl;
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

		//math::XYZTLorentzVector tau_p4_temp = tau.p4();
    	//math::XYZTLorentzVector piChar_p4_temp = tau.leadChargedHadrCand()->p4();

		++nTauC;
		allTauPt->Fill(tau.pt());
		std::pair <int, double> PairAbsIsoPt (tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"), tau.pt());
		TauIdMap.insert(std::pair<int, std::pair<double, double>> (i, PairAbsIsoPt));

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

    if (nTauC > 1) {
		//std::cout << "Tau candidates:" << std::endl;
		for (auto MapId_iter = TauIdMap.begin(); MapId_iter !=  TauIdMap.end(); ++MapId_iter) {
			if ((*MapId_iter).second.first < tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") ||
				(*MapId_iter).second.first == tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits") && (*MapId_iter).second.second < tau.pt()) {
				index_2 = (*MapId_iter).first;
			}
		}
	} else {
		index_2 = null;
	}

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
	} else {
		DiPhoton_pt  = null;
		DiPhoton_eta = null;
		DiPhoton_phi = null;
		DiPhoton_m   = null;
	}

	nPi0      = 1;
	nPhotons  = nPhotons_temp;
	nGamma    = tau.signalGammaCands().size();
	nTau      = taus->size();
	tau_found = 1;

	if (index_2 >= 0) {
		auto& tau2  = (*taus)[index_2];
		tau2_pt     = tau2.pt();
		tau2_eta    = tau2.eta();
		tau2_phi    = tau2.phi();
		tau2_dm     = tau2.decayMode();
		tau2_dz     = (pv_position - tau2.vertex()).R();
		tau2_q      = tau2.charge();
		tau2_m      = tau2.mass();
		tau2_absIso = tau2.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
	} else {
		tau2_pt     = null;
		tau2_eta    = null;
		tau2_phi    = null;
		tau2_dm     = null;
		tau2_dz     = null;
		tau2_q      = null;
		tau2_m      = null;
		tau2_absIso = null;
	}

	return true;
};

bool TreeMakerMiniAOD::FindGenTau(const edm::Event& event) {

	if (monitoringGen) std::cout << std::endl << "!!! Generated particles !!!" << std::endl;

	gentau_found   = 0;
	genTauFromW   = null;
	genTauMother  = null;
	dR            = null;
	bquark1dR     = null;
	bquark2dR     = null;
	leptondR      = null;
	W_pt          = null;
	gentau_pt     = null;
	genPiChar_pt  = null;
	genPi0_pt     = null;
	nutau_pt      = null;
	nuW_pt        = null;
	nunu_pt       = null;
	double nuW_m  = 0.;

	if (monitoringGen) std::cout << "gentau Flag 0" << std::endl;

	//edm::Handle<pat::PackedCandidateCollection> genParticles;
	edm::Handle<reco::GenParticleCollection> genParticles;
	event.getByToken(GenParticleToken_, genParticles);
	if (!genParticles.isValid()) return false;

	if (monitoringGen) std::cout << "gentau Flag 1" << std::endl;

	const int pdg_tau      = 15;
	const int pdg_pi0      = 111;
	const int pdg_pi1      = 211;
	const int pdg_W        = 24;
	const int pdg_nu_tau   = 16;
	const int pdg_electron = 11;
	const int pdg_nu_ele   = 12;
	const int pdg_mu       = 13;
	const int pdg_nu_mu    = 14;
	const int pdg_bquark   = 5;
	const int pdg_tquark   = 6;

	const reco::Candidate* tau    = nullptr;
	const reco::Candidate* nu_W   = nullptr;
	const reco::Candidate* nu_tau = nullptr;
	const reco::Candidate* pi0    = nullptr;
	const reco::Candidate* pi1    = nullptr;
	double dRmin = null;
	// new vectors of gen particles
	std::vector<const reco::Candidate*> leptons;
	std::vector<const reco::Candidate*> bquarks;
	std::vector<const reco::Candidate*> Wbosons;
	std::vector<const reco::Candidate*> tquarks;


	// Look for tau among generated particles
	for (auto& particle: *genParticles) {
#define cut(condition) if (!(condition)) continue;
		// look for the tau -> pi+ pi0 neutrino decay most oriented towards
		// the reconstructed tau (if present)
		cut(abs(particle.pdgId()) == pdg_tau);
		if (monitoringGen) std::cout << "Gentau in GP was found" << std::endl;

		for (unsigned i = 0; i < particle.numberOfDaughters(); ++i) {
			const reco::Candidate* daughter = particle.daughter(i);
			int id = abs(daughter->pdgId());
			if (id == pdg_pi0)
				pi0 = daughter;
			else if (id == pdg_pi1)
				pi1 = daughter;
			else if (id == pdg_nu_tau)
				nu_tau = daughter;
		};

		cut(particle.pt() > tauPtMin);
		if (monitoringGen) std::cout << "gentau pt" << std::endl;
		cut(TMath::Abs(particle.eta()) < tauEtaMax);
		if (monitoringGen) std::cout << "gentau eta" << std::endl;
		cut((pv_position - particle.vertex()).R() < tauDzMax);
		if (monitoringGen) std::cout << "gentau dz" << std::endl;
		//cut(pi1->pt() > piPtMin);

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
		if (monitoringGen) {
			std::cout << "gentau pt = " << particle.pt() << std::endl;
			std::cout << "gentau eta = " << particle.eta() << std::endl;
			std::cout << "gentau phi = " << particle.phi() << std::endl;
			std::cout << "gentau dRmin = " << dRmin << std::endl;
		}
	};

	if (!tau) return false;

	if (monitoringGen) std::cout << "gentau Flag 2" << std::endl;

	gentau_pt  = tau->pt();
	gentau_px  = tau->px();
	gentau_py  = tau->py();
	gentau_pz  = tau->pz();
	gentau_eta = tau->eta();
	gentau_phi = tau->phi();

	dR           = dRmin;
	gentau_found = 1;
	genTauFromW  = 0;

	// Look for W which is mother particle for tau
	for (auto p = tau->mother(); p; p = p->mother()) {
		genTauMother = p->pdgId();
		if (monitoringGen) std::cout << "gnetau mother pdg ID = " << genTauMother << std::endl;
		if (abs(p->pdgId()) == pdg_W) {
			genTauFromW = 1;
			W_pt        = p->pt();
			Wbosons.push_back(p);
			for (unsigned l = 0; l < p->numberOfDaughters(); l++) {
				const reco::Candidate* Wdaughter = p->daughter(l);
				if (monitoringGen) std::cout << "W_daughter(" << l << ") = " << Wdaughter->pdgId() << std::endl;
				if (abs(Wdaughter->pdgId()) == pdg_nu_tau) {
					nu_W = Wdaughter;
				} else continue;
			}
			for (auto p1 = p->mother(); p1; p1 = p1->mother()) {
				if (abs(p1->pdgId()) == pdg_tquark) {
					if (monitoringGen) {
						std::cout << "t quark from tau decay chain found" << std::endl;
						std::cout << "pt  = " << p1->pt() << std::endl;
						std::cout << "eta = " << p1->eta() << std::endl;
						std::cout << "phi = " << p1->phi() << std::endl;
					}
					//FirstTCharge = p1->charge();
					tquarks.push_back(p1);
					for (unsigned m = 0; m < p1->numberOfDaughters(); m++) {
						const reco::Candidate* tdaughter = p1->daughter(m);
						if (abs(tdaughter->pdgId()) == pdg_bquark) {
							bquarks.push_back(tdaughter);
							double bdR1 = sqrt(deltaR2(tdaughter->eta(), tdaughter->phi(), BJet1_eta, BJet1_phi));
							double bdR2 = sqrt(deltaR2(tdaughter->eta(), tdaughter->phi(), BJet2_eta, BJet2_phi));
							bquark1dR = TMath::Min(bdR1, bdR2);
							if (monitoringGen) {
								std::cout << "b quark from tau decay chain found" << std::endl;
								std::cout << "pt  = " << tdaughter->pt() << std::endl;
								std::cout << "eta = " << tdaughter->eta() << std::endl;
								std::cout << "phi = " << tdaughter->phi() << std::endl;
								std::cout << "dR1 = " << bdR1 << std::endl;
								std::cout << "dR2 = " << bdR2 << std::endl;
							}
						} else continue;
					}
					break;
				}
			}
			break;
		};
	};

	// Look for second t-quark and it's daughter particles among generated particles
	for (auto& particle: *genParticles) {
		if (abs(particle.pdgId()) == pdg_tquark) {
			//if (tquarks.size() > 0 && (*tquarks[0]).charge() * particle.charge() >= 0) continue;
			const reco::Candidate* tquark = &particle;
			GetLastDaughter(tquark);
			bool sametquark = false;
			if (tquarks.size() > 0) {
				for (unsigned Nt = 0; Nt < tquarks.size(); Nt++) {
					double ttdR = sqrt(deltaR2((*tquarks[Nt]).eta(), (*tquarks[Nt]).phi(), particle.eta(), particle.phi()));
					if (monitoringGen) std::cout << "ttdR = " << ttdR << std::endl;
					if (ttdR < 0.4) {
						sametquark = true;
						break;
					}
				}
			}
			if (sametquark) continue;
			if (monitoringGen) {
				std::cout << "Second t quark found ID = " << particle.pdgId() << std::endl;
				std::cout << "pt  = " << particle.pt() << std::endl;
				std::cout << "eta = " << particle.eta() << std::endl;
				std::cout << "phi = " << particle.phi() << std::endl;
			}
			tquarks.push_back(tquark);
			for (unsigned l = 0; l < tquark->numberOfDaughters(); l++) {
				const reco::Candidate* tdaughter = tquark->daughter(l);
				if (monitoringGen) {
					std::cout << "t_daughter(" << l << ") ID = " << tdaughter->pdgId() << std::endl;
				}
				if (abs(tdaughter->pdgId()) == pdg_bquark) {
					bquarks.push_back(tdaughter);
					double bdR1 = sqrt(deltaR2(tdaughter->eta(), tdaughter->phi(), BJet1_eta, BJet1_phi));
					double bdR2 = sqrt(deltaR2(tdaughter->eta(), tdaughter->phi(), BJet2_eta, BJet2_phi));
					bquark2dR = TMath::Min(bdR1, bdR2);
					if (monitoringGen) {
						std::cout << "Second b quark found" << std::endl;
						std::cout << "pt  = " << tdaughter->pt() << std::endl;
						std::cout << "eta = " << tdaughter->eta() << std::endl;
						std::cout << "phi = " << tdaughter->phi() << std::endl;
						std::cout << "dR1 = " << bdR1 << std::endl;
						std::cout << "dR2 = " << bdR2 << std::endl;
					}
				} else if (abs(tdaughter->pdgId()) == pdg_W) {
					Wbosons.push_back(tdaughter);
					if (monitoringGen) std::cout << "Second W boson found" << std::endl;
					GetLastDaughter(tdaughter);
					for (unsigned m = 0; m < tdaughter->numberOfDaughters(); m++) {
						const reco::Candidate* Wdaughter = tdaughter->daughter(m);
						if (monitoringGen) std::cout << "W_daughter(" << m << ") ID = " << Wdaughter->pdgId() << std::endl;
						if (abs(Wdaughter->pdgId()) == pdg_electron || abs(Wdaughter->pdgId()) == pdg_mu || abs(Wdaughter->pdgId()) == pdg_tau) {
							leptons.push_back(Wdaughter);
							double lepdR1 = sqrt(deltaR2(Wdaughter->eta(), Wdaughter->phi(), lepton1_eta, lepton1_phi));
							double lepdR2 = sqrt(deltaR2(Wdaughter->eta(), Wdaughter->phi(), lepton2_eta, lepton2_phi));
							double dRtau  = sqrt(deltaR2(Wdaughter->eta(), Wdaughter->phi(), gentau_eta, gentau_phi));
							leptondR = TMath::Min(lepdR1, lepdR2); 
							if (monitoringGen) {
								std::cout << "Second lepton from W found with ID = " << Wdaughter->pdgId() << std::endl;
								std::cout << "pt  = " << Wdaughter->pt() << std::endl;
								std::cout << "eta = " << Wdaughter->eta() << std::endl;
								std::cout << "phi = " << Wdaughter->phi() << std::endl;
								std::cout << "dR1 = " << lepdR1 << std::endl;
								std::cout << "dR2 = " << lepdR2 << std::endl;
								std::cout << "dRtau = " << dRtau << std::endl;
							}
						} else continue;
					}
				} else continue;
			}
			//if (sqrt(deltaR2((*tquarks[0]).eta(), (*tquarks[0]).phi(), particle.eta(), particle.phi())) > 0.5) break;
		}
	}

	if (monitoringGen) std::cout << "gentau Flag 3" << std::endl;

	if (pi1) {
		genPiChar_pt  = pi1->pt();
		genPiChar_px  = pi1->px();
		genPiChar_py  = pi1->py();
		genPiChar_pz  = pi1->pz();
		genPiChar_eta = pi1->eta();
		genPiChar_phi = pi1->phi();
	} else {
		genPiChar_pt  = null;
		genPiChar_px  = null;
		genPiChar_py  = null;
		genPiChar_pz  = null;
		genPiChar_eta = null;
		genPiChar_phi = null;
	}

	double genPi0_m;

	if (pi0) {
		genPi0_pt  = pi0->pt();
		genPi0_px  = pi0->px();
		genPi0_py  = pi0->py();
		genPi0_pz  = pi0->pz();
		genPi0_eta = pi0->eta();
		genPi0_phi = pi0->phi();
	} else {
		genPi0_pt  = null;
		genPi0_px  = null;
		genPi0_py  = null;
		genPi0_pz  = null;
		genPi0_eta = null;
		genPi0_phi = null;
	}

	if (monitoringGen) std::cout << "gentau Flag 4" << std::endl;

	if (!nu_W) return false; 

	if (monitoringGen) std::cout << "gentau Flag 5" << std::endl;

	if (genTauFromW > 0) {
		nuW_pt  = nu_W->pt();
		nuW_px  = nu_W->px();
		nuW_py  = nu_W->py();
		nuW_pz  = nu_W->pz();
		nuW_eta = nu_W->eta();
		nuW_phi = nu_W->phi();
		nutau_pt  = nu_tau->pt();
		nutau_px  = nu_tau->px();
		nutau_py  = nu_tau->py();
		nutau_pz  = nu_tau->pz();
		nutau_eta = nu_tau->eta();
		nutau_phi = nu_tau->phi();
		double nutau_m = nu_tau->mass();
		double nuW_m = nu_W->mass();

		TLorentzVector gennutau, gennuW;
		gennutau.SetPtEtaPhiM(nutau_pt, nutau_eta, nutau_phi, nutau_m);
		gennuW.SetPtEtaPhiM(nuW_pt, nuW_eta, nuW_phi, nuW_m);
		nunu_pt = (gennutau + gennuW).Pt();
	}

	return true;
};

void TreeMakerMiniAOD::GenTauDM (const edm::Event& iEvent) {

  //edm::Handle<pat::PackedCandidateCollection> GP;
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
        if (monitoringGen) std::cout << "GenTauDM particle number " << i << std::endl;

        //std::vector<SimpleParticle> tau_daughters_simple;
	std::vector<MySimpleParticle> tau_daughters;
        for (unsigned k = 0; k < (*GP)[i].numberOfDaughters(); k++) {
          MySimpleParticle tp((*GP)[i].daughter(k)->p4().Px(), (*GP)[i].daughter(k)->p4().Py(), (*GP)[i].daughter(k)->p4().Pz(), (*GP)[i].daughter(k)->p4().E(), (*GP)[i].daughter(k)->pdgId());
          tau_daughters.push_back(tp);
          if (monitoringGen) std::cout << "Tau_daughter[" << k << "] = " << (*GP)[i].daughter(k)->pdgId() << std::endl;
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
          	gentau_dm = 3;
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
      if (channel == -1 && abs((*GP)[i].pdgId())==15 && monitoringGen) {
        std::cout << "Unidentified decay channel of tau was found" << std::endl;
	    std::cout << "List of daughters:" << std::endl;
        for (unsigned t = 0; t < (*GP)[i].numberOfDaughters(); t++) {
          std::cout << "TauDaughter[" << t << "] pdgId = " << (*GP)[i].daughter(t)->pdgId() << std::endl;
        }
      }
    }
  }
  /*
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
  */
};

bool TreeMakerMiniAOD::AddMET(const edm::Event& event) {
	edm::Handle<pat::METCollection> mets;
	event.getByToken(MetCollectionToken_, mets);
	if (!mets.isValid() || !mets->size()) return false;

	auto& MET = mets->front();
	met = MET.pt();
	met_phi = MET.phi();
	met_eta = MET.eta();
	met_energy = MET.energy();
	met_significance = MET.significance();
	met_mEtSig       = MET.mEtSig();

	TLorentzVector METp4, tau_p4;
	METp4.SetPtEtaPhiE(met, met_eta, met_phi, met_energy);
	tau_p4.SetPtEtaPhiM(tau_pt, tau_eta, tau_phi, tau_m);
	tauMET_mass = (METp4 + tau_p4).M();
	//std::cout << "tau + MET inv mass = " << tauMET_mass << std::endl;

	if (met < METcut) return false;
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

// Jets analysis

bool TreeMakerMiniAOD::JetPtSum (const edm::Event& event) {

	if (monitoringBJets) std::cout << std::endl << "!!! Jets !!!" << std::endl;

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

    // Jets
	edm::Handle<pat::JetCollection> jets;
	event.getByToken(JetCollectionToken_, jets);

	if (!jets.isValid()) return false;

	int i = 0;
	for (auto& jet: *jets) {
		if (TMath::Abs(jet.eta()) > 3) continue;
		if (jet.pt() > 15) {
			i++;
			if (monitoringJets) std::cout << std::endl << "Jet[" << i << "] Pt = " << jet.pt() << std::endl;
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
				if (monitoringJets) std::cout << "Vertex number   = " << (*Map_iter1).first << "  entries = " << (*Map_iter1).second.first << "  SumPt = " << (*Map_iter1).second.second << std::endl;
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
				/*
				BJet.Eta = jet.eta();
				BJet.Phi = jet.phi();
				BJet.Pt = jet.pt();
				BJet.Bprob = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb");
				BJet.FourMomentum = jet.p4();
				BJet.E = jet.energy();
				BJet.FromPV = JetFromPV;
				*/
				looseBJets.push_back(BJet);
				//std::pair<double, bool> pairBprobFromPV(JetDeepCSV_prob_b + JetDeepCSV_prob_bb, JetFromPV);
				//looseBJetsProb.insert(std::pair<int, std::pair<double, bool>>(j, pairBprobFromPV));
				//
				// Attempted to obtain btag information (parameters of b-meson), not successful
				/*
				std::cout << "Tag Info Lables (" << jet.tagInfoLabels().size() << "):" << std::endl;
				std::cout << "Tag Info (CATop, " << jet.tagInfo("CATop")->hasTracks() << "):" << std::endl;
				for (auto TagInfoLabelsIter = jet.tagInfoLabels().begin(); TagInfoLabelsIter != jet.tagInfoLabels().end(); ++TagInfoLabelsIter) {	
					std::cout << *TagInfoLabelsIter << std::endl;
					std::cout << "Has tag info = " << jet.hasTagInfo(*TagInfoLabelsIter) << std::endl;
				}
				*/
				/*
				if ()
				reco::TaggingVariableList TaggingVars = jet.tagInfo().TaggingVariables();
				std::cout << "Tagging variables list:" << std::endl;
				for (unsigned k = 0; k < TaggingVars.size(); k++) {
					std::cout << TaggingVars[k].first << " = " << TaggingVars[k].second << std::endl;
				}
				*/
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

	return nJets20 > 0;
};

bool TreeMakerMiniAOD::AddLepton (const edm::Event& event) {

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
	lepton2_pt       = null;
	lepton2_eta      = null;
	lepton2_phi      = null;
	lepton2_dz       = null;
	lepton2_flavor   = null;
	lepton2_charge   = null;
	lepton2_E        = null;
	lepton2_trackIso = null;
	m_ll             = null;

	edm::Handle<pat::ElectronCollection> electrons;
	event.getByToken(ElectronCollectionToken_, electrons);
	if (!electrons.isValid()) return false;

	edm::Handle<pat::MuonCollection> muons;
	event.getByToken(MuonCollectionToken_, muons);
	if (!muons.isValid()) return false;

	std::vector<LeptonCandidate> LepCandidates;

  // Select electrons and muons which satisfy following criteria
	// Not in bJets, pt above threshold, eta below threshold
	if (electrons.isValid() && !(*electrons).empty()) {
		//
		//if (monitoringLeptons) std::cout << "Isolated electrons:" << std::endl;
		for (auto& electron: *electrons) {
#define cut(condition) if (!(condition)) continue;
			cut(electron.pt() > MuElePtMin);
			cut(abs(electron.eta()) < EtaMax);
			cut((pv_position - electron.vertex()).R() < tauDzMax);
			double deltaRJet1 = sqrt(deltaR2(BJet1_eta, BJet1_phi, electron.eta(), electron.phi()));
			double deltaRJet2 = sqrt(deltaR2(BJet2_eta, BJet2_phi, electron.eta(), electron.phi()));
			double deltaRTau  = sqrt(deltaR2(tau_eta, tau_phi, electron.eta(), electron.phi()));
			//double dphiJet1 = dphi(Jet1_phi, electron.phi());
			//double dphiJet2 = dphi(Jet2_phi, electron.phi());
			if (deltaRTau < 0.3) return false;
			if (deltaRJet1 > 0.5 && deltaRJet2 > 0.5) {
				bool SameLepton = false;
				if(!LepCandidates.empty()) {
					for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
						if (abs(LepCandidates[ilep].Phi - electron.phi()) < 0.3) SameLepton = true;
					}
				}
				cut(!SameLepton);
				nLeptonCandidates++;
				LeptonCandidate Lepton(electron, pv_position);
				/*
				Lepton.Pt = electron.pt();
				Lepton.FourMomentum = electron.p4();
				Lepton.Eta = electron.eta();
				Lepton.Phi = electron.phi();
				Lepton.Dz = (pv_position - electron.vertex()).R();
				Lepton.Flavor = 11;
				Lepton.Charge = electron.charge();
				// Isolation variables
				Lepton.caloIso = electron.caloIso();
				Lepton.ecalIso = electron.ecalIso();
				Lepton.hcalIso = electron.hcalIso();
				Lepton.puppiChargedHadronIso = electron.puppiChargedHadronIso();
				Lepton.puppiNeutralHadronIso = electron.puppiNeutralHadronIso();
				Lepton.puppiNoLeptonsChargedHadronIso = electron.puppiNoLeptonsChargedHadronIso();
				Lepton.puppiNoLeptonsNeutralHadronIso = electron.puppiNoLeptonsNeutralHadronIso();
				Lepton.puppiNoLeptonsPhotonIso = electron.puppiNoLeptonsPhotonIso();
				Lepton.puppiPhotonIso = electron.puppiPhotonIso();
				Lepton.trackIso = electron.trackIso();
				*/
				LepCandidates.push_back(Lepton);
	//delete Lepton;
			}
#undef cut
		}
	}

	if (muons.isValid() && !(*muons).empty()) {
		//if (monitoringLeptons) std::cout << "Isolated muons:" << std::endl;
		for (auto& muon: *muons) {
#define cut(condition) if (!(condition)) continue;
			cut(muon.pt() > MuElePtMin);
			cut(abs(muon.eta()) < EtaMax);
			cut((pv_position - muon.vertex()).R() < tauDzMax);
			double deltaRJet1 = sqrt(deltaR2(BJet1_eta, BJet1_phi, muon.eta(), muon.phi()));
			double deltaRJet2 = sqrt(deltaR2(BJet2_eta, BJet2_phi, muon.eta(), muon.phi()));
			double deltaRTau  = sqrt(deltaR2(tau_eta, tau_phi, muon.eta(), muon.phi()));
			//double dphiJet1 = dphi(Jet1_phi, muon.phi());
			//double dphiJet2 = dphi(Jet2_phi, muon.phi());
			if (deltaRTau < 0.3) return false;
			if (deltaRJet1 > 0.5 && deltaRJet2 > 0.5) {
				bool SameLepton = false;
				if(!LepCandidates.empty()) {
					for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
						if (abs(LepCandidates[ilep].Phi - muon.phi()) < 0.3) SameLepton = true;
					}
				}
				cut(!SameLepton);
				LeptonCandidate Lepton(muon, pv_position);
				nLeptonCandidates++;
				/*
				Lepton.Pt = muon.pt();
				Lepton.FourMomentum = muon.p4();
				Lepton.Eta = muon.eta();
				Lepton.Phi = muon.phi();
				Lepton.Dz = (pv_position - muon.vertex()).R();
				Lepton.Flavor = 13;
				Lepton.Charge = muon.charge();
				// Isolation variables
				Lepton.caloIso = muon.caloIso();
				Lepton.ecalIso = muon.ecalIso();
				Lepton.hcalIso = muon.hcalIso();
				Lepton.puppiChargedHadronIso = muon.puppiChargedHadronIso();
				Lepton.puppiNeutralHadronIso = muon.puppiNeutralHadronIso();
				Lepton.puppiNoLeptonsChargedHadronIso = muon.puppiNoLeptonsChargedHadronIso();
				Lepton.puppiNoLeptonsNeutralHadronIso = muon.puppiNoLeptonsNeutralHadronIso();
				Lepton.puppiNoLeptonsPhotonIso = muon.puppiNoLeptonsPhotonIso();
				Lepton.puppiPhotonIso = muon.puppiPhotonIso();
				Lepton.trackIso = muon.trackIso();
				*/
				LepCandidates.push_back(Lepton);
	//delete Lepton;
			}
#undef cut
		}
	}

	if (nLeptonCandidates < 1) return false;

	SortLeptons(LepCandidates);

	if (monitoringLeptons) std::cout << "List of selected leptons:" << std::endl;
	for(unsigned ilep = 0; ilep < LepCandidates.size(); ilep++) {
		if (monitoringLeptons) {
			std::cout << "lepton " << ilep << std::endl;
			LepCandidates[ilep].Monitoring();
			/*
			std::cout << "flavor = " << LepCandidates[ilep].Flavor << std::endl;
			std::cout << "charge = " << LepCandidates[ilep].Charge << std::endl;
			std::cout << "dz     = " << LepCandidates[ilep].Dz << std::endl;
			std::cout << "pt     = " << LepCandidates[ilep].Pt << std::endl;
			std::cout << "eta    = " << LepCandidates[ilep].Eta << std::endl;
			std::cout << "phi    = " << LepCandidates[ilep].Phi << std::endl;
			std::cout << "delta(BJet1) = " << sqrt(deltaR2(BJet1_eta, BJet1_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
			std::cout << "delta(BJet2) = " << sqrt(deltaR2(BJet2_eta, BJet2_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
			std::cout << "delta(Tau)   = " << sqrt(deltaR2(tau_eta, tau_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
			std::cout << "caloIso = " << LepCandidates[ilep].caloIso << std::endl;
			std::cout << "ecalIso = " << LepCandidates[ilep].ecalIso << std::endl;
			std::cout << "hcalIso = " << LepCandidates[ilep].hcalIso << std::endl;
			std::cout << "puppiChargedHadronIso          = " << LepCandidates[ilep].puppiChargedHadronIso << std::endl;
			std::cout << "puppiNeutralHadronIso          = " << LepCandidates[ilep].puppiNeutralHadronIso << std::endl;
			std::cout << "puppiNoLeptonsChargedHadronIso = " << LepCandidates[ilep].puppiNoLeptonsChargedHadronIso << std::endl;
			std::cout << "puppiNoLeptonsNeutralHadronIso = " << LepCandidates[ilep].puppiNoLeptonsNeutralHadronIso << std::endl;
			std::cout << "puppiNoLeptonsPhotonIso        = " << LepCandidates[ilep].puppiNoLeptonsPhotonIso << std::endl;
			std::cout << "puppiPhotonIso                 = " << LepCandidates[ilep].puppiPhotonIso << std::endl;
			std::cout << "trackIso                       = " << LepCandidates[ilep].trackIso << std::endl;
			*/
			std::cout << "delta(BJet1) = " << sqrt(deltaR2(BJet1_eta, BJet1_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
			std::cout << "delta(BJet2) = " << sqrt(deltaR2(BJet2_eta, BJet2_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
			std::cout << "delta(Tau)   = " << sqrt(deltaR2(tau_eta, tau_phi, LepCandidates[ilep].Eta, LepCandidates[ilep].Phi)) << std::endl;
		}
	}

	lepton1_pt       = LepCandidates[0].Pt;
	lepton1_eta      = LepCandidates[0].Eta;
	lepton1_phi      = LepCandidates[0].Phi;
	//lepton1_dz       = LepCandidates[0].Dz;
	lepton1_flavor   = LepCandidates[0].Flavor;
	lepton1_charge   = LepCandidates[0].Charge;
	lepton1_E        = LepCandidates[0].FourMomentum.E();
	lepton1_trackIso = LepCandidates[0].trackIso;
	if (LepCandidates.size() > 1) {
		lepton2_pt       = LepCandidates[1].Pt;
		lepton2_eta      = LepCandidates[1].Eta;
		lepton2_phi      = LepCandidates[1].Phi;
		//lepton2_dz       = LepCandidates[1].Dz;
		lepton2_flavor   = LepCandidates[1].Flavor;
		lepton2_charge   = LepCandidates[1].Charge;
		lepton2_E        = LepCandidates[1].FourMomentum.E();
		lepton2_trackIso = LepCandidates[1].trackIso;
		m_ll             = (LepCandidates[0].FourMomentum + LepCandidates[1].FourMomentum).M();
	}

	LepCandidates.clear();

	return (nLeptonCandidates > 0);
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
};

bool TreeMakerMiniAOD::DecaychannelMatch(std::vector<MySimpleParticle> &particles, int p1, int p2, int p3, int p4, int p5, int p6) {
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
};

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
