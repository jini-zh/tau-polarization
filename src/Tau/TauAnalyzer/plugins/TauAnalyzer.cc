// -*- C++ -*-
//
// Package:    Tau/TauAnalyzer
// Class:      TauAnalyzer
// 
/**\class TauAnalyzer TauAnalyzer.cc Tau/TauAnalyzer/plugins/TauAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrei Oskin
//         Created:  Sat, 14 Sep 2019 01:15:10 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <stdio.h>

// user include files
//#include "GeneratorInterface/TauolaInterface/interface/TauSpinnerCMS.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

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
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Tau/TauAnalyzer/plugins/MySimpleParticle.h"

using namespace std;
using namespace edm;
using namespace reco;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TauAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TauAnalyzer(const edm::ParameterSet&);
      ~TauAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

      bool FindGenTau(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      bool DecaychannelMatch(vector<MySimpleParticle> &particles, int p1, int p2=0, int p3=0, int p4=0, int p5=0, int p6=0);

      // ----------member data ---------------------------
};


edm::EDGetTokenT<double> TauSpinnerWTToken_;
edm::EDGetTokenT<double> TauSpinnerWTFlipToken_;
edm::EDGetTokenT<double> TauSpinnerWThminusToken_;
edm::EDGetTokenT<double> TauSpinnerWThplusToken_;
edm::EDGetTokenT<bool>   TauSpinnerWTisValidToken_;
edm::EDGetTokenT<std::vector<double>> PolarimetricVectorToken_;

edm::EDGetTokenT<reco::GenParticleCollection> GenParticleToken_;

TTree * tree;

int t_Run;
int t_Event;
int null;

double WT;
double WTFlip;
double WThminus;
double WThplus;
bool   WTisValid;
int channel;

double gentau_pt;
double gentau_e;
int tau_charge;
double hx, hy, hz, ht;

int PichargedNumber = 0;
double picharged_pt;
double picharged_e;
double x_MC;
std::vector<double> PolarimetricVector;

TH1D * WTHist;
TH1D * WTFlipHist;
TH1D * WThminusHist;
TH1D * WThplusHist;
TH1D * xHist;
TH1D * xHist_Weight;
TH1D * xHist_unWeight;
TH1D * xHist_reWeight;

TH1D * x_tauP_polMinus;
TH1D * x_tauM_polMinus;
TH1D * x_tauP_polPlus;
TH1D * x_tauM_polPlus;

// Poalrinetric vector components
TH1D * hx_Hist;
TH1D * hy_Hist;
TH1D * hz_Hist;
TH1D * ht_Hist;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TauAnalyzer::TauAnalyzer(const edm::ParameterSet& iConfig)

{

  /*
  std::string WTCollection        = iConfig.getParameter<std::string>("WTCollection");
  std::string WTFlipCollection    = iConfig.getParameter<std::string>("WTFlipCollection");
  std::string WThminusCollection  = iConfig.getParameter<std::string>("WThminusCollection");
  std::string WThplusCollection   = iConfig.getParameter<std::string>("WThplusCollection");
  std::string WTisValidCollection = iConfig.getParameter<std::string>("WTisValidCollection");

  TauSpinnerWTToken_        = consumes<double>(edm::InputTag(WTCollection));
  TauSpinnerWTFlipToken_    = consumes<double>(edm::InputTag(WTFlipCollection));
  TauSpinnerWThminusToken_  = consumes<double>(edm::InputTag(WThminusCollection));
  TauSpinnerWThplusToken_   = consumes<double>(edm::InputTag(WThplusCollection));
  TauSpinnerWTisValidToken_ = consumes<bool>(edm::InputTag(WTisValidCollection));
  */

  null            = iConfig.getParameter<int>("null");

  TauSpinnerWTToken_        = consumes<double>(iConfig.getParameter<edm::InputTag>("WTCollection"));
  TauSpinnerWTFlipToken_    = consumes<double>(iConfig.getParameter<edm::InputTag>("WTFlipCollection"));
  TauSpinnerWThminusToken_  = consumes<double>(iConfig.getParameter<edm::InputTag>("WThminusCollection"));
  TauSpinnerWThplusToken_   = consumes<double>(iConfig.getParameter<edm::InputTag>("WThplusCollection"));
  TauSpinnerWTisValidToken_ = consumes<bool>(iConfig.getParameter<edm::InputTag>("WTisValidCollection"));
  PolarimetricVectorToken_  = consumes<std::vector<double>>(iConfig.getParameter<edm::InputTag>("PolarimetricVector"));

  GenParticleToken_         = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));

   //now do what ever initialization is needed
   //usesResource("TFileService");

}


TauAnalyzer::~TauAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  t_Run   = iEvent.id().run();
  t_Event = iEvent.id().event();

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
  edm::Handle<std::vector<double>> PolarimetricVectorHandle;
  iEvent.getByToken(PolarimetricVectorToken_, PolarimetricVectorHandle);

  WT        = *WTHandle;
  WTFlip    = *WTFlipHandle;
  WThminus  = *WThminusHandle;
  WThplus   = *WThplusHandle;
  WTisValid = *WTisValidHandle;
  PolarimetricVector = *PolarimetricVectorHandle;

  hx = PolarimetricVector[0];
  hy = PolarimetricVector[1];
  hz = PolarimetricVector[2];
  ht = PolarimetricVector[3];

  bool genTauFound = FindGenTau(iEvent, iSetup);
  if(!genTauFound) {
    return;
  }
  if((WT >= 0) && (WT <= 2)) {
    WTHist->Fill(WT);
    WTFlipHist->Fill(WTFlip);
  }

  std::cout << "N(tau -> pi + tau_nu) = " << PichargedNumber << endl;

  tree->Fill();

}

bool TauAnalyzer::FindGenTau(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool tauFound = false;
  edm::Handle<reco::GenParticleCollection> GP;
  iEvent.getByToken(GenParticleToken_, GP);
  
  if(GP.isValid()) {
    std::cout << "GenParticles collection is valid" << std::endl;
    for(unsigned i = 0 ; i < GP->size() ; i++) {
      if (abs((*GP)[i].pdgId())==15) {
        channel = -1;
	int tau_pdgid = (*GP)[i].pdgId();

        //std::vector<SimpleParticle> tau_daughters_simple;
	std::vector<MySimpleParticle> tau_daughters;
        for (unsigned k = 0; k < (*GP)[i].numberOfDaughters(); k++) {
          MySimpleParticle tp((*GP)[i].daughter(k)->p4().Px(), (*GP)[i].daughter(k)->p4().Py(), (*GP)[i].daughter(k)->p4().Pz(), (*GP)[i].daughter(k)->p4().E(), (*GP)[i].daughter(k)->pdgId());
          tau_daughters.push_back(tp);
        }

        vector<int>  pdgid;
        double *HH = new double[4];

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
	  PichargedNumber++;
          channel = 0;
	  std::cout << "tau -> pi + tau_nu = " << endl;
          if (abs(pdgid[1]) == 321) channel = 4;
	  double Tau_e = (*GP)[i].energy();
	  double PiCh_e = tau_daughters[1].e();
	  gentau_pt = (*GP)[i].pt();
	  gentau_e = (*GP)[i].energy();
	  picharged_pt = sqrt(TMath::Power(tau_daughters[1].px(), 2) + TMath::Power(tau_daughters[1].py(), 2));
	  picharged_e = tau_daughters[1].e();
	  x_MC = PiCh_e/Tau_e;
	  std::cout << "Polarimetric vector of tau: " << std::endl;
	  for (unsigned l = 0; l < PolarimetricVector.size(); l++) {
	    std::cout << "HH[" << l << "] = " << PolarimetricVector[l] << std::endl;
          }
	  if ((x_MC >= 0) && (x_MC <= 1)) {
              std::cout << "WT = " << WT << std::endl;
              std::cout << "WTFlip = " << WTFlip << std::endl;
              xHist->Fill(x_MC);
              xHist_Weight->Fill(x_MC, WT);
              xHist_unWeight->Fill(x_MC, 1/WT);
              xHist_reWeight->Fill(x_MC, WTFlip);
	      hx_Hist->Fill(hx);
	      hy_Hist->Fill(hy);
	      hz_Hist->Fill(hz);
	      ht_Hist->Fill(ht);
            }
          tauFound = true;
        }
        /*           
        else if ( pdgid.size() == 3  &&
            (
              ( tau_pdgid == 15 && DecaychannelMatch(tau_daughters, 16,-211, 111) ) ||
              ( tau_pdgid ==-15 && DecaychannelMatch(tau_daughters,-16, 211, 111) )
            )
          ) {
          channel = 1;
        } else if( pdgid.size() == 4 &&
            (
              ( tau_pdgid == 15 && DecaychannelMatch(tau_daughters, 16, 111, 111,-211) ) ||
              ( tau_pdgid ==-15 && DecaychannelMatch(tau_daughters,-16, 111, 111, 211) )
            )
          ) {
          channel = 2;
        } else if( pdgid.size() == 4 &&
            (
              ( tau_pdgid == 15 && DecaychannelMatch(tau_daughters, 16,-211,-211, 211) ) ||
              ( tau_pdgid ==-15 && DecaychannelMatch(tau_daughters,-16, 211, 211,-211) )
            )DecaychannelMatch
          ) {
          channel = 3;
        } else if( pdgid.size() == 5 &&
            (
              ( tau_pdgid == 15 && DecaychannelMatch(tau_daughters, 16,-211,-211, 211, 111) ) ||
              ( tau_pdgid ==-15 && DecaychannelMatch(tau_daughters,-16, 211, 211,-211, 111) )
            )
          ) {
          channel = 7;
        }
	*/
      }
    }
  }
  return tauFound;
}


bool TauAnalyzer::DecaychannelMatch(vector<MySimpleParticle> &particles, int p1, int p2, int p3, int p4, int p5, int p6) {
  // Copy pdgid-s of all particles
  vector<int> list;
  
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

  vector<MySimpleParticle> newList;
  
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

// ------------ method called once each job just before starting event loop  ------------
void 
TauAnalyzer::beginJob()
{
  edm::Service<TFileService> FS;
  tree = FS->make<TTree>("tree", "tree", 1);

  tree->Branch("WT",&WT,"WT/D");
  tree->Branch("WTFlip",&WTFlip,"WTFlip/D");
  tree->Branch("WThminus",&WThminus,"WThminus/D");
  tree->Branch("WThplus",&WThplus,"WThplus/D");
  tree->Branch("WTisValid",&WTisValid,"WTisValid/D");
  tree->Branch("gentau_pt",&gentau_pt,"gentau_pt/D");
  tree->Branch("gentau_e",&gentau_e,"gentau_e/D");
  tree->Branch("picharged_pt",&picharged_pt,"picharged_pt/D");
  tree->Branch("picharged_e",&picharged_e,"picharged_e/D");
  tree->Branch("x_MC",&x_MC,"x_MC/D");
  tree->Branch("hx",&hx,"hx");
  tree->Branch("hy",&hy,"hy");
  tree->Branch("hz",&hz,"hz");
  tree->Branch("ht",&ht,"ht");

  WTHist        = FS->make<TH1D>("WTHist","WTHist",100,0,2);
  WTFlipHist    = FS->make<TH1D>("WTFlipHist","WTFlipHist",100,0,2);
  WThminusHist  = FS->make<TH1D>("WThminusHist","WThminusHist",100,0,2);
  WThplusHist   = FS->make<TH1D>("WThplusHist","WThplusHist",100,0,2);
  xHist           = FS->make<TH1D>("xHist","x",100,0,1);
  xHist_unWeight  = FS->make<TH1D>("xHist_unWeight","x/WT",100,0,1);
  xHist_Weight    = FS->make<TH1D>("xHist_Weight","x*WT",100,0,1);
  xHist_reWeight  = FS->make<TH1D>("xHist_reWeight","x*WTFlip",100,0,1);
  hx_Hist         = FS->make<TH1D>("hx_Hist","hx",100,-1,1);
  hy_Hist         = FS->make<TH1D>("hy_Hist","hy",100,-1,1);
  hz_Hist         = FS->make<TH1D>("hz_Hist","hz",100,-1,1);
  ht_Hist         = FS->make<TH1D>("ht_Hist","ht",100,-1,1);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void TauAnalyzer::beginRun(const edm::Run& Run, const edm::EventSetup& Setup) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  //descriptions.addDefault(desc);
  descriptions.add("tauanalyzer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyzer);
