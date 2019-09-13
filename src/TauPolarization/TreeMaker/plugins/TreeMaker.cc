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
#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"


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
	
	bool TriggerOK(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	bool AddTau(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	bool FindGenTau(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	bool CheckMuon(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	bool CheckElectron(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	bool AddMET(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	bool AddVertex(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	void CountTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	bool JetPtSum(const edm::Event& iEvent, const edm::EventSetup& iSetup);

	double dPhiFrom2P(double Px1, double Py1, double Px2, double Py2) {
		double prod = Px1*Px2 + Py1*Py2;
		double mod1 = TMath::Sqrt(Px1*Px1+Py1*Py1);
		double mod2 = TMath::Sqrt(Px2*Px2+Py2*Py2);
		double cosDPhi = prod/(mod1*mod2);
		return TMath::ACos(cosDPhi);
	}


	double deltaR(double Eta1, double Phi1, double Eta2, double Phi2) {
		double Px1 = TMath::Cos(Phi1);
		double Py1 = TMath::Sin(Phi1);
		double Px2 = TMath::Cos(Phi2);
		double Py2 = TMath::Sin(Phi2);
		double dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2);
		double dEta = Eta1 - Eta2;
		double d_R = TMath::Sqrt(dPhi*dPhi+dEta*dEta);
		return d_R;
	}

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
	edm::EDGetTokenT<reco::PFMETCollection> MetCollectionToken_;
	edm::EDGetTokenT<reco::VertexCollection> PVToken_;
	edm::EDGetTokenT<reco::GenParticleCollection> GenParticleToken_;
	
	edm::EDGetTokenT<reco::TrackCollection> TrackToken_;
	
	edm::EDGetTokenT<reco::PFTauDiscriminator> Token_absIso;

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
	

	double tau_found;
	double gentau_found;
	double dR;
	double genTauFromW;
	int genTauMother;

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

	double jetPtSum10, jetPtSum15, jetPtSum20;
	double RHT10, RHT15, RHT20;
	int nJets20;
	int nPi0;
	double met;
	double met_phi;
	
	double m_t;
	
	double dPhi;
	
	math::XYZPoint pv_position;

	//////////////////////////////////////////////////////
	bool isMC;
	double tauPtMin;
	double piPtMin;
	double tauEtaMax;
	double tauDzMax;
	double null;

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
	tauPtMin 					= iConfig.getParameter<double>("tauPtMin");
	piPtMin 					= iConfig.getParameter<double>("piPtMin");
	tauEtaMax 					= iConfig.getParameter<double>("tauEtaMax");
	tauDzMax 					= iConfig.getParameter<double>("tauDzMax");
	null 						= iConfig.getParameter<double>("null");
	
	trigNames 					= iConfig.getParameter<std::vector<std::string>>("Triggers");
	triggerEvent_				= edm::InputTag("hltTriggerSummaryAOD","","HLT");
	theTriggerResultsLabel		= edm::InputTag("TriggerResults","","HLT");

	std::string tauCollection         = iConfig.getParameter<std::string>("tauCollection");
	std::string muonCollection        = iConfig.getParameter<std::string>("muonCollection");
	std::string electronCollection    = iConfig.getParameter<std::string>("electronCollection");
	std::string jetCollection         = iConfig.getParameter<std::string>("jetCollection");
	std::string metCollection         = iConfig.getParameter<std::string>("metCollection");
	std::string vertexCollection      = iConfig.getParameter<std::string>("vertexCollection");
	std::string genParticleCollection = iConfig.getParameter<std::string>("genParticleCollection");
	std::string trackCollection       = iConfig.getParameter<std::string>("trackCollection");
	
	std::string absIsoDiscriminator					= iConfig.getParameter<std::string>("absIsoDiscriminator");
	std::string looseCombinedIsoDiscriminator		= iConfig.getParameter<std::string>("looseCombinedIsoDiscriminator");
	std::string mediumCombinedIsoDiscriminator		= iConfig.getParameter<std::string>("mediumCombinedIsoDiscriminator");
	std::string tightCombinedIsoDiscriminator		= iConfig.getParameter<std::string>("tightCombinedIsoDiscriminator");
	std::string looseMvaIsoDiscriminator			= iConfig.getParameter<std::string>("looseMvaIsoDiscriminator");
	std::string mediumMvaIsoDiscriminator			= iConfig.getParameter<std::string>("mediumMvaIsoDiscriminator");
	std::string tightMvaIsoDiscriminator			= iConfig.getParameter<std::string>("tightMvaIsoDiscriminator");
	std::string looseMuonRejectionDiscriminator		= iConfig.getParameter<std::string>("looseMuonRejectionDiscriminator");
	std::string tightMuonRejectionDiscriminator		= iConfig.getParameter<std::string>("tightMuonRejectionDiscriminator");
	std::string looseElectronRejectionDiscriminator	= iConfig.getParameter<std::string>("looseElectronRejectionDiscriminator");
	std::string tightElectronRejectionDiscriminator	= iConfig.getParameter<std::string>("tightElectronRejectionDiscriminator");
	
	
	tok_trigEvt					= consumes<trigger::TriggerEvent>(triggerEvent_);
	tok_trigRes					= consumes<edm::TriggerResults>(theTriggerResultsLabel);
	TauCollectionToken_ 		= consumes<reco::PFTauCollection>(edm::InputTag(tauCollection));
	MuonCollectionToken_ 		= consumes<reco::MuonCollection>(edm::InputTag(muonCollection));
	ElectronCollectionToken_	= consumes<reco::GsfElectronCollection>(edm::InputTag(electronCollection));
	JetCollectionToken_ 		= consumes<reco::PFJetCollection>(edm::InputTag(jetCollection));
	MetCollectionToken_ 		= consumes<reco::PFMETCollection>(edm::InputTag(metCollection));
	PVToken_ 					= consumes<reco::VertexCollection>(edm::InputTag(vertexCollection));
	GenParticleToken_ 			= consumes<reco::GenParticleCollection>(edm::InputTag(genParticleCollection));
	TrackToken_					= consumes<reco::TrackCollection>(edm::InputTag(trackCollection));

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

}


TreeMaker::~TreeMaker() {
	 // do anything here that needs to be done at desctruction time
	 // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace edm;
	t_Run   = iEvent.id().run();
	t_Event = iEvent.id().event();
	
	bool triggerOK = TriggerOK(iEvent,iSetup);
	if(!triggerOK) {
		return;
	}

	bool vertexFound = AddVertex(iEvent,iSetup);
	if (!vertexFound) { //cout << "V" << endl;
		return;
	}
	 
	bool muonsFound = CheckMuon(iEvent, iSetup);
	if (muonsFound) { //cout << "M" << endl;
		return;
	}

	bool electronFound = CheckElectron(iEvent, iSetup);
	if (electronFound) { //cout << "E" << endl;
		return;
	}
	 
	bool metFound = AddMET(iEvent,iSetup);
	if (!metFound) { //cout << "met" << endl;
		return;
	}

 	bool tauFound = AddTau(iEvent,iSetup);
	bool genTauFound = FindGenTau(iEvent, iSetup);
	if(!tauFound&&!genTauFound) { //cout << "T" << endl;
		return;
	}
	//cout << "Tau found" << endl;
	CountTracks(iEvent, iSetup);

	bool jetsFound = JetPtSum(iEvent, iSetup);
	if(!jetsFound) { //cout << "J" << endl;
		return;
	}

	TLorentzVector tauLV, piChLV, pi0LV;
	tauLV.SetPtEtaPhiM(tau_pt,tau_eta,tau_phi,tau_m);
	piChLV.SetPtEtaPhiM(piChar_pt,piChar_eta,piChar_phi,piChar_m);
	pi0LV.SetPtEtaPhiM(piZero_pt,piZero_eta,piZero_phi,piZero_m);


   	if(tau_pt>null) {
   		pipiMass = (piChLV+pi0LV).M();
   		//upsilon = (piChLV.E() - pi0LV.E())/tau_pt;
   		upsilon = 2*piChar_pt/tau_pt - 1;
	   	RHT10 = tau_pt/jetPtSum10;
	   	RHT15 = tau_pt/jetPtSum15;
	   	RHT20 = tau_pt/jetPtSum20;
   		dPhi = dPhiFrom2P(tauLV.Px(),tauLV.Py(),cos(met_phi),sin(met_phi));
		m_t = sqrt(2*tau_pt*met*(1-cos(dPhiFrom2P(tauLV.Px(),tauLV.Py(),cos(met_phi),sin(met_phi)))));
   	} else {
   		pipiMass = null;
		upsilon = null;
		RHT10 = null;
		RHT15 = null;
		RHT20 = null;
		dPhi = null;
		m_t = null;
   	}

  	tree->Fill();
}


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

	tree->Branch("genTauMother", &genTauMother, "genTauMother/I");
	

	tree->Branch("jetPtSum10", &jetPtSum10, "jetPtSum10/D");
	tree->Branch("RHT10", &RHT10, "RHT10/D");
	tree->Branch("jetPtSum15", &jetPtSum15, "jetPtSum15/D");
	tree->Branch("RHT15", &RHT15, "RHT15/D");
	tree->Branch("jetPtSum20", &jetPtSum20, "jetPtSum20/D");
	tree->Branch("RHT20", &RHT20, "RHT20/D");
	tree->Branch("nJets20", &nJets20, "nJets20/I");
	
	tree->Branch("met", &met, "met/D");
	tree->Branch("met_phi", &met_phi, "met_phi/D");
	tree->Branch("m_t", &m_t, "m_t/D");
	tree->Branch("dPhi", &dPhi, "dPhi/D");
	
	tree->Branch("nVtx",&nVtx,"nVtx/I");
	tree->Branch("nTrks",&nTrks,"nTrks/I");
	tree->Branch("nTau",&nTau,"nTau/I");
	tree->Branch("nTauC",&nTauC,"nTauC/I");

	tree->Branch("nPi0",&nPi0,"nPi0/I");
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
bool TreeMaker::TriggerOK(const edm::Event& iEvent, const edm::EventSetup& iSetup){
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

bool TreeMaker::AddTau(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	//tau collection
	using namespace reco;

	bool tauFound = false;
	nTau = 0;
	nTauC = 0;
	tau_found = 0;
	tau_pt = null;
	tau_eta = null;
	tau_phi = null;
	tau_dm = null;
	tau_dz = null;
	tau_q = null;
	tau_m = null;
	tau_absIso = null;
	tau_looseCombinedIso = null;
	tau_mediumCombinedIso = null;
	tau_tightCombinedIso = null;
	tau_looseMvaIso = null;
	tau_mediumMvaIso = null;
	tau_tightMvaIso = null;
	tau_looseMuonRejection = null;
	tau_tightMuonRejection = null;
	tau_looseElectronRejection = null;
	tau_tightElectronRejection = null;
	dR = null;

	piChar_pt = null;
	piChar_eta = null;
	piChar_phi = null;
	piChar_q = null;
	piChar_m = null;
	piZero_pt = null;
	piZero_eta = null;
	piZero_phi = null;
	piZero_m = null;

	nPi0 = null;

	edm::Handle<reco::PFTauCollection> Taus;
	iEvent.getByToken(TauCollectionToken_, Taus);
	
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

	iEvent.getByToken(Token_absIso, absIso);
	iEvent.getByToken(Token_looseCombinedIso, looseCombinedIso);
	iEvent.getByToken(Token_mediumCombinedIso, mediumCombinedIso);
	iEvent.getByToken(Token_tightCombinedIso, tightCombinedIso);
	iEvent.getByToken(Token_looseMvaIso, looseMvaIso);
	iEvent.getByToken(Token_mediumMvaIso, mediumMvaIso);
	iEvent.getByToken(Token_tightMvaIso, tightMvaIso);
	iEvent.getByToken(Token_looseMuonRejection, looseMuonRejection);
	iEvent.getByToken(Token_tightMuonRejection, tightMuonRejection);
	iEvent.getByToken(Token_looseElectronRejection, looseElectronRejection);
	iEvent.getByToken(Token_tightElectronRejection, tightElectronRejection);
	

	if(Taus.isValid()) {
		nTau = Taus->size();
		double minIso = 1e+10;
		double ptMax = 0;
		for(unsigned i = 0 ; i < Taus->size() ; i++) {
		  
			double _tau_pt = (*Taus)[i].pt();
			double _tau_eta = (*Taus)[i].eta();
			double _tau_phi = (*Taus)[i].phi();
			double _tau_m = (*Taus)[i].mass();
			double _tau_q = (*Taus)[i].charge();
			double _tau_dm = (*Taus)[i].decayMode();

			double xx = pv_position.x() - (*Taus)[i].vx();
			double yy = pv_position.y() - (*Taus)[i].vy();
			double zz = pv_position.z() - (*Taus)[i].vz();
			double _tau_dz = sqrt(xx*xx+yy*yy+zz*zz);
				
			double _tau_absIso = absIso->value(i);
			double _tau_looseCombinedIso = looseCombinedIso->value(i);
			double _tau_mediumCombinedIso = mediumCombinedIso->value(i);
			double _tau_tightCombinedIso = tightCombinedIso->value(i);
			double _tau_looseMvaIso = looseMvaIso->value(i);
			double _tau_mediumMvaIso = mediumMvaIso->value(i);
			double _tau_tightMvaIso = tightMvaIso->value(i);
			double _tau_looseMuonRejection = looseMuonRejection->value(i);
			double _tau_tightMuonRejection = tightMuonRejection->value(i);
			double _tau_looseElectronRejection = looseElectronRejection->value(i);
			double _tau_tightElectronRejection = tightElectronRejection->value(i);

			double _piChar_pt = null;
			double _piChar_eta = null;
			double _piChar_phi = null;
			double _piChar_q = null;
			double _piChar_m = null;

			double _piZero_pt = null;
			double _piZero_eta = null;
			double _piZero_phi = null;
			double _piZero_m = null;

			double _nPi0 = null;
			
			bool tauFlag_dm = false;
			//if((_tau_dm==1)&&((*Taus)[i].signalPFChargedHadrCands().size()==1)&&((*Taus)[i].signalPiZeroCandidates().size()==1)) {
			if(((*Taus)[i].signalPFChargedHadrCands().size()==1)&&((*Taus)[i].signalPiZeroCandidates().size()>0)) {
				_piChar_pt = (*Taus)[i].signalPFChargedHadrCands()[0]->pt();
				_piChar_eta = (*Taus)[i].signalPFChargedHadrCands()[0]->eta();
				_piChar_phi = (*Taus)[i].signalPFChargedHadrCands()[0]->phi();
				_piChar_q = (*Taus)[i].signalPFChargedHadrCands()[0]->charge();
				_piChar_m = (*Taus)[i].signalPFChargedHadrCands()[0]->mass();

				_piZero_pt = (*Taus)[i].signalPiZeroCandidates()[0].pt();
				_piZero_eta = (*Taus)[i].signalPiZeroCandidates()[0].eta();
				_piZero_phi = (*Taus)[i].signalPiZeroCandidates()[0].phi();
				_piZero_m = (*Taus)[i].signalPiZeroCandidates()[0].mass();
        	
        		_nPi0 = (*Taus)[i].signalPiZeroCandidates().size();

        		tauFlag_dm = true;
			}
			
			if((_tau_pt > tauPtMin)&&(_piChar_pt > piPtMin)&&(TMath::Abs(_tau_eta)<tauEtaMax)){
		    nTauC++;
		    allTauPt->Fill(_tau_pt);
		  }
			bool tauFlag = (_tau_pt > tauPtMin)&&(_piChar_pt > piPtMin)&&(TMath::Abs(_tau_eta)<tauEtaMax)&&(TMath::Abs(_tau_dz)<tauDzMax)&&((_tau_absIso<minIso)||((_tau_absIso==minIso)&&(_tau_pt>ptMax)));
			
			if(tauFlag&&tauFlag_dm) {
				tau_found = 1;
				tau_pt = _tau_pt;
	        	tau_eta = _tau_eta;
	        	tau_phi = _tau_phi;
	        	tau_dm = _tau_dm;
	        	tau_dz = _tau_dz;
	        	tau_q = _tau_q;
	        	tau_m = _tau_m;
	        	tau_absIso = _tau_absIso;
	        	tau_looseCombinedIso = _tau_looseCombinedIso;
	        	tau_mediumCombinedIso = _tau_mediumCombinedIso;
	        	tau_tightCombinedIso = _tau_tightCombinedIso;
	        	tau_looseMvaIso = _tau_looseMvaIso;
	        	tau_mediumMvaIso = _tau_mediumMvaIso;
	        	tau_tightMvaIso = _tau_tightMvaIso;
	        	tau_looseMuonRejection = _tau_looseMuonRejection;
	        	tau_tightMuonRejection = _tau_tightMuonRejection;
	        	tau_looseElectronRejection = _tau_looseElectronRejection;
	        	tau_tightElectronRejection = _tau_tightElectronRejection;
				  
	        	piChar_pt = _piChar_pt;
				piChar_eta = _piChar_eta;
				piChar_phi = _piChar_phi;
				piChar_q = _piChar_q;
				piChar_m = _piChar_m;

				piZero_pt = _piZero_pt;
				piZero_eta = _piZero_eta;
				piZero_phi = _piZero_phi;
				piZero_m = _piZero_m;

				nPi0 = _nPi0;

				minIso = tau_absIso;
				ptMax = _tau_pt;

				tauFound = true;
			}
		}
	}

	return tauFound;

}

bool TreeMaker::FindGenTau(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	bool tauFound = false;
	gentau_found = 0;
	genTauFromW = null;
	genTauMother = null;
	edm::Handle<reco::GenParticleCollection> GP;
	iEvent.getByToken(GenParticleToken_, GP);
	double dRmin = 10;
	
	if(GP.isValid()) {
		for(unsigned i = 0 ; i < GP->size() ; i++) {
			if (abs((*GP)[i].pdgId())==15) {
				if((*GP)[i].numberOfDaughters()==3) {
					if(
	        			((abs((*GP)[i].daughter(1)->pdgId())==111)||(abs((*GP)[i].daughter(1)->pdgId())==211))&&
	        			((abs((*GP)[i].daughter(2)->pdgId())==111)||(abs((*GP)[i].daughter(2)->pdgId())==211))&&
	        			(abs((*GP)[i].daughter(1)->pdgId())+abs((*GP)[i].daughter(2)->pdgId())==322)
	        		) {
	      	
	      			//std::cout<<"::::::::"<<(*GP)[i].daughter(0)->pdgId()<<" "<<(*GP)[i].daughter(1)->pdgId()<<" "<<(*GP)[i].daughter(2)->pdgId()<<" "<<abs((*GP)[i].daughter(1)->pdgId())+abs((*GP)[i].daughter(2)->pdgId())<<std::endl;

						double gentau_pt = (*GP)[i].pt();
						double gentau_phi = (*GP)[i].phi();
						double gentau_eta = (*GP)[i].eta();
						double xxx = pv_position.x() - (*GP)[i].vx();
						double yyy = pv_position.y() - (*GP)[i].vy();
						double zzz = pv_position.z() - (*GP)[i].vz();
						double gentau_dz = sqrt(xxx*xxx+yyy*yyy+zzz*zzz);
						double _dR;
						if(tau_found) {
							_dR = deltaR(gentau_eta,gentau_phi,tau_eta,tau_phi);
						} else {
							_dR = null;
						}
						bool genTauFlag = (gentau_pt > tauPtMin)&&(TMath::Abs(gentau_eta)<tauEtaMax)&&(TMath::Abs(gentau_dz)<tauDzMax)&&(_dR<dRmin);
						bool genDFlag1 = (abs((*GP)[i].daughter(1)->pdgId())==211)&&((*GP)[i].daughter(1)->pt() > piPtMin);
						bool genDFlag2 = (abs((*GP)[i].daughter(2)->pdgId())==211)&&((*GP)[i].daughter(2)->pt() > piPtMin);
						if(genTauFlag&&(genDFlag1 || genDFlag2)) {
							dR = _dR;
							gentau_found = 1;
							genTauMother = (*GP)[i].mother()->pdgId();
							if(abs((*GP)[i].mother()->pdgId())==24) {
								genTauFromW = 1;
								genTauMother = (*GP)[i].mother()->pdgId();
							} else {
								if((abs((*GP)[i].mother()->pdgId())==15)&&(abs((*GP)[i].mother()->mother()->pdgId())==24)){
									genTauFromW = 1;
									genTauMother = (*GP)[i].mother()->mother()->pdgId();
								} else {
									if((abs((*GP)[i].mother()->mother()->pdgId())==15)&&(abs((*GP)[i].mother()->mother()->mother()->pdgId())==24)){
										genTauFromW = 1;
										genTauMother = (*GP)[i].mother()->mother()->mother()->pdgId();
									} else {
										genTauFromW = 0;
									}
								}
							}
							tauFound = true;
						}
					}
				}
			}
		}
	}
	return tauFound;
}


bool TreeMaker::CheckMuon(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	
	using namespace reco;
	bool muonFound = false;
	edm::Handle<reco::MuonCollection> Muons;
	iEvent.getByToken(MuonCollectionToken_, Muons);

	if(Muons.isValid()) {
		for(unsigned i = 0 ; i < Muons->size() ; i++) {
			if ((*Muons)[i].pt()>15) muonFound = true;
		}
	}
	return muonFound;
}


bool TreeMaker::CheckElectron(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	
	using namespace reco;
	bool electronFound = false;
	edm::Handle<reco::GsfElectronCollection> Electrons;
	iEvent.getByToken(ElectronCollectionToken_, Electrons);

	if(Electrons.isValid()) {
		for(unsigned i = 0 ; i < Electrons->size() ; i++) {
			if ((*Electrons)[i].pt()>15) electronFound = true;
		}
	}
	return electronFound;
}

bool TreeMaker::AddMET(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	bool foundMET = false;
	edm::Handle<reco::PFMETCollection> recoMet;
	iEvent.getByToken(MetCollectionToken_, recoMet);

	if(recoMet->size() > 0) {
		met = (*recoMet)[0].pt();
		met_phi = (*recoMet)[0].phi();
		foundMET = true;
	}
	return foundMET;
}


bool TreeMaker::AddVertex(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

	bool vertexfound = false;
	nVtx = 0;
	pv_position = math::XYZPoint(0.,0.,0.);
	edm::Handle<reco::VertexCollection> Vertex;
	iEvent.getByToken(PVToken_, Vertex);
	
	if(Vertex.isValid()) {
	  nVtx = Vertex->size();
		if (Vertex->size()>0) {
			pv_position = (*Vertex)[0].position();
			vertexfound = true;
		}
	}
	return vertexfound;
}


bool TreeMaker::JetPtSum(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	bool jetsFound = false;
	double sum10 = 0;
	double sum15 = 0;
	double sum20 = 0;
	int nJets = 0;
	edm::Handle<reco::PFJetCollection> Jets;
	iEvent.getByToken(JetCollectionToken_, Jets);

	if(Jets.isValid()) {
		for(unsigned i = 0 ; i < Jets->size() ; i++) {
			if(	((*Jets)[i].pt()>10)	&&	(abs((*Jets)[i].eta())<3)	){
				sum10 = sum10 + (*Jets)[i].pt();
			}
			if(	((*Jets)[i].pt()>15)	&&	(abs((*Jets)[i].eta())<3)	){
				sum15 = sum15 + (*Jets)[i].pt();
			}
			if(	((*Jets)[i].pt()>20)	&&	(abs((*Jets)[i].eta())<3)	){
				sum20 = sum20 + (*Jets)[i].pt();
				nJets++;
			}
		}
		jetPtSum10 = sum10;
		jetPtSum15 = sum15;
		jetPtSum20 = sum20;
		nJets20 = nJets;
		if(sum10*sum15*sum20!=0){
			jetsFound = true;
		}
	}
	return jetsFound;
}

void TreeMaker::CountTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  nTrks = 0;
	edm::Handle<reco::TrackCollection> Tracks;
	iEvent.getByToken(TrackToken_, Tracks);
	if(Tracks.isValid()) {
	  nTrks = Tracks->size();
	}
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
