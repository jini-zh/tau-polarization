#include "GeneratorInterface/TauolaInterface/interface/TauSpinnerCMS.h"

//MC-TESTER header files
#include "Tauola/Tauola.h"
#include "TauSpinner/tau_reweight_lib.h"
#include "TauSpinner/Tauola_wrapper.h"
#include "GeneratorInterface/TauolaInterface/interface/read_particles_from_HepMC.h"
#include "TLorentzVector.h"

#include "CLHEP/Random/RandomEngine.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ServiceRegistry/interface/RandomEngineSentry.h"

using namespace edm;
using namespace TauSpinner;

CLHEP::HepRandomEngine *TauSpinnerCMS::fRandomEngine = nullptr;
bool TauSpinnerCMS::isTauSpinnerConfigure = false;
bool TauSpinnerCMS::fInitialized = false;

unsigned EventNumber = 0;
unsigned TauCounter = 0;

TauSpinnerCMS::TauSpinnerCMS(const ParameterSet &pset)
    : EDProducer(),
      isReco_(pset.getParameter<bool>("isReco")),
      isTauolaConfigured_(pset.getParameter<bool>("isTauolaConfigured")),
      isLHPDFConfigured_(pset.getParameter<bool>("isLHPDFConfigured")),
      LHAPDFname_(pset.getUntrackedParameter("LHAPDFname", (string)("MSTW2008nnlo90cl.LHgrid"))),
      CMSEnergy_(pset.getParameter<double>("CMSEnergy"))  //GeV
      ,
      gensrc_(pset.getParameter<edm::InputTag>("gensrc")),
      MotherPDGID_(pset.getUntrackedParameter("MotherPDGID", (int)(24))), // -1 - all particles
      Ipol_(pset.getUntrackedParameter("Ipol", (int)(0))), // For Z/gamma case only???
      nonSM2_(pset.getUntrackedParameter("nonSM2", (int)(1))),
      nonSMN_(pset.getUntrackedParameter("nonSMN", (int)(0))),
      monitoring(pset.getParameter<bool>("monitoring")),
      roundOff_(pset.getUntrackedParameter("roundOff", (double)(0.01))) {
  //usesResource(edm::uniqueSharedResourceName());
  usesResource(edm::SharedResourceNames::kTauola);

  produces<bool>("TauSpinnerWTisValid").setBranchAlias("TauSpinnerWTisValid");
  produces<bool>("PhotonEmisson1").setBranchAlias("PhotonEmisson1");
  produces<bool>("PhotonEmisson2").setBranchAlias("PhotonEmisson2");
  produces<double>("TauSpinnerWT").setBranchAlias("TauSpinnerWT");
  produces<double>("TauSpinnerWTFlip").setBranchAlias("TauSpinnerWTFlip");
  produces<double>("TauSpinnerWThplus").setBranchAlias("TauSpinnerWThplus");
  produces<double>("TauSpinnerWThminus").setBranchAlias("TauSpinnerWThminus");
  produces<int>("TauSpinnerMother").setBranchAlias("TauSpinnerMother");
  produces<double>("TauPolarisation").setBranchAlias("TauPolarisation");
  produces<std::vector<double>>("PolarimetricVector").setBranchAlias("PolarimetricVector");

  if (isReco_) {
    GenParticleCollectionToken_ = consumes<reco::GenParticleCollection>(gensrc_);
  } else {
    hepmcCollectionToken_ = consumes<HepMCProduct>(gensrc_);
  }
}

void TauSpinnerCMS::beginJob(){};
void TauSpinnerCMS::endJob(){};

void TauSpinnerCMS::initialize() {
  // Now for Tauola and TauSpinner
  Tauolapp::Tauola::setRandomGenerator(TauSpinnerCMS::flat);
  if (!isTauolaConfigured_) {
    Tauolapp::Tauola::initialize();
  }
  if (!isLHPDFConfigured_) {
    LHAPDF::initPDFSetByName(LHAPDFname_);
  }
  if (!isTauSpinnerConfigure) {
    isTauSpinnerConfigure = true;
    bool Ipp = true;  // for pp collisions
    // Initialize TauSpinner
    //Ipol - polarization of input sample
    //nonSM2 - nonstandard model calculations
    //nonSMN
    TauSpinner::initialize_spinner(Ipp, Ipol_, nonSM2_, nonSMN_, CMSEnergy_);
  }
  fInitialized = true;
}

void TauSpinnerCMS::produce(edm::Event &e, const edm::EventSetup &iSetup) {
  RandomEngineSentry<TauSpinnerCMS> randomEngineSentry(this, e.streamID());
  if (!fInitialized)
    initialize();

  Tauolapp::Tauola::setRandomGenerator(
      TauSpinnerCMS::flat);  // rest tauola++ random number in case other modules use tauola++
  Tauolapp::jaki_.ktom = 1;  // rest for when you run after tauola
  double WT = 1.0;
  double WTFlip = 1.0;
  double polSM = -999;  //range [-1,1]
  SimpleParticle X, tau, tau2;
  std::vector<SimpleParticle> tau_daughters, tau_daughters2;
  int stat(0);
  int MotherPDGId = 0;
  std::vector<double> HHvector;
  bool Photon1 = false;
  bool Photon2 = false;

  EventNumber++;
  if (monitoring) {
    std::cout << endl << "### Event Number = " << EventNumber << endl;
    std::cout << "isReco_ = " << isReco_ << endl;
  }

  if (isReco_) {
    // readParticlesfromReco returns 0 if there is decay with tau
    // returns 1 if there are no taus
    stat = readParticlesfromReco(e, X, tau, tau2, tau_daughters, tau_daughters2);

    if (stat == 23) {
      Photon1 = true;
      Photon2 = true;
    } else if (stat == 2) {
      Photon1 = true;
      Photon2 = false;
    } else if (stat == 3) {
      Photon1 = false;
      Photon2 = true;
    }

    if (monitoring) {
      std::cout << "stat = " << stat << endl;
      std::cout << "Total number of tau leptons = " << TauCounter << endl;
    }

  } else {
    edm::Handle<HepMCProduct> evt;
    e.getByToken(hepmcCollectionToken_, evt);
    //Get EVENT
    HepMC::GenEvent *Evt = new HepMC::GenEvent(*(evt->GetEvent()));
    // см. документацию TauSpinner
    stat = readParticlesFromHepMC(Evt, X, tau, tau2, tau_daughters, tau_daughters2);
  }
  MotherPDGId = X.pdgid();
  if (MotherPDGID_ < 0 || abs(X.pdgid()) == MotherPDGID_) {
    if (stat != 1) {
      // Determine the weight
      // Wheights for W+ and H+
      if (abs(X.pdgid()) == 24 || abs(X.pdgid()) == 37 || abs(X.pdgid()) == 9900024) {
      	if (abs(X.pdgid()) == 24) {
      	  if (monitoring) std::cout << "W+ is mother" << std::endl;
      	} else if (abs(X.pdgid()) == 37) {
      	  if (monitoring) std::cout << "H+ is mother" << std::endl;
      	} else if (abs(X.pdgid()) == 9900024) {
      	  if (monitoring) std::cout << "W_R is mother (" << X.pdgid() << ")" << std::endl;
      	  // set PDGId of H+ instead of W_R
      	  X.setPdgid(SignumPDG(X.pdgid()) * 37);
	  //std::cout << "New PDGId = " << X.pdgid() << std::endl;
      	}
        TLorentzVector tau_1r(0, 0, 0, 0);
        TLorentzVector tau_1(tau.px(), tau.py(), tau.pz(), tau.e());
        for (unsigned int i = 0; i < tau_daughters.size(); i++) {
          tau_1r += TLorentzVector(
              tau_daughters.at(i).px(), tau_daughters.at(i).py(), tau_daughters.at(i).pz(), tau_daughters.at(i).e());
        }
        if (fabs(tau_1r.M() - tau_1.M()) < roundOff_) {
          WT = TauSpinner::calculateWeightFromParticlesWorHpn(
              X, tau, tau2, tau_daughters);  // note that tau2 is tau neutrino
          polSM = getTauSpin();
          WTFlip = (2.0 - WT) / WT;

          //-----------------------------------------------------------------------------------------------------
	  // Prepare kinematic fo HH and calculate HH for Tau daughters
	  double phi2 = 0.0, theta2 = 0.0;
	  // Create Particles from SimpleParticles because method calculateWeightFromParticlesWorHpn takes SimpleParticles
	  // while prepareKinematicForHH and calculateHH take Particles
 
	  Particle X_new     (   X.px(),    X.py(),    X.pz(),    X.e(),    X.pdgid() );
	  Particle tau_new   ( tau.px(),  tau.py(),  tau.pz(),  tau.e(),  tau.pdgid() );
	  Particle nu_tau_new(tau2.px(), tau2.py(), tau2.pz(), tau2.e(), tau2.pdgid() );
 
	  vector<Particle> tau_daughters_new;
 
	  // tau pdgid
	  int tau_pdgid = tau.pdgid();
 
	  // Create list of tau daughters
	  for(unsigned int i = 0; i < tau_daughters.size(); i++) {
	  Particle pp(tau_daughters[i].px(),
          tau_daughters[i].py(),
	  tau_daughters[i].pz(),
	  tau_daughters[i].e(),
	  tau_daughters[i].pdgid() );
 
	  tau_daughters_new.push_back(pp);
	  }
	  prepareKinematicForHH   (tau_new, nu_tau_new, tau_daughters_new, &phi2, &theta2);
	  double *HH = calculateHH(tau_pdgid, tau_daughters_new, phi2, theta2);
	  HHvector.clear();
	  for (unsigned k = 0; k < 4; k++) {
	    HHvector.push_back(HH[k]);
	    if (monitoring) std::cout << "HH[" << k << "] = " << HH[k] << std::endl;
          }
          //-----------------------------------------------------------------------------------------------------

          // Monitoring of values
          if (monitoring) {
            std::cout << endl << "WT     = " << WT << endl;
            std::cout << endl << "polSM  = " << polSM << endl;
            std::cout << endl << "WTFlip = " << WTFlip << endl;
          }
        }
      // Wheights for Z0, gamma, H0, A0
      } else if (X.pdgid() == 25 || X.pdgid() == 36 || X.pdgid() == 22 || X.pdgid() == 23) {
      	if (abs(X.pdgid()) == 25) {
      	  if (monitoring) std::cout << "h0 is mother" << std::endl;
      	} else if (abs(X.pdgid()) == 36) {
      	  if (monitoring) std::cout << "A0 is mother" << std::endl;
      	} else if (abs(X.pdgid()) == 22) {
      	  if (monitoring) std::cout << "gamma is mother" << std::endl;
      	} else if (abs(X.pdgid()) == 23) {
      	  if (monitoring) std::cout << "Z0 is mother" << std::endl;
      	}
        TLorentzVector tau_1r(0, 0, 0, 0), tau_2r(0, 0, 0, 0);
        TLorentzVector tau_1(tau.px(), tau.py(), tau.pz(), tau.e()), tau_2(tau2.px(), tau2.py(), tau2.pz(), tau2.e());
        for (unsigned int i = 0; i < tau_daughters.size(); i++) {
          tau_1r += TLorentzVector(
              tau_daughters.at(i).px(), tau_daughters.at(i).py(), tau_daughters.at(i).pz(), tau_daughters.at(i).e());
        }
        for (unsigned int i = 0; i < tau_daughters2.size(); i++) {
          tau_2r += TLorentzVector(tau_daughters2.at(i).px(),
                                   tau_daughters2.at(i).py(),
                                   tau_daughters2.at(i).pz(),
                                   tau_daughters2.at(i).e());
        }

        if (fabs(tau_1r.M() - tau_1.M()) < roundOff_ && fabs(tau_2r.M() - tau_2.M()) < roundOff_) {
          WT = TauSpinner::calculateWeightFromParticlesH(X, tau, tau2, tau_daughters, tau_daughters2);
          polSM = getTauSpin();
          if (monitoring) {
            std::cout << "polSM   = " << polSM << std::endl;
            std::cout << "This is weight and spin for initial mother particle" << std::endl;
            std::cout << "WT      = " << WT << std::endl;
          }
          if (X.pdgid() == 25 || X.pdgid() == 22 || X.pdgid() == 23) {
            if (X.pdgid() == 25)
              X.setPdgid(23);
            if (X.pdgid() == 22 || X.pdgid() == 23)
              X.setPdgid(25);

            // Attribute kinematics of taus to another mother particle decay
            //       h0 ---> Z0
            // Z0/gamma ---> h0
            double WTother = TauSpinner::calculateWeightFromParticlesH(X, tau, tau2, tau_daughters, tau_daughters2);
            WTFlip = WTother / WT;
            if (monitoring) {
              std::cout << "Parameters for mother particle with PDGID = " << X.pdgid() << std::endl;
              std::cout << "WTother = " << WTother << std::endl;
              std::cout << "WTFlip  = " << WTFlip << std::endl;
            }
          }
        }
      } else {
        cout << "TauSpinner: WARNING: Unexpected PDG for tau mother: " << X.pdgid() << endl;
      }
    }
  }
  bool isValid = true;
  if (!(0 <= WT && WT < 10)) {
    isValid = false;
    WT = 1.0;
    WTFlip = 1.0;
  }
  std::unique_ptr<bool> TauSpinnerWeightisValid(new bool);
  *TauSpinnerWeightisValid = isValid;
  e.put(std::move(TauSpinnerWeightisValid), "TauSpinnerWTisValid");

  // PDGId of mother particle
  std::unique_ptr<int> TauSpinnerMother(new int);
  *TauSpinnerMother = MotherPDGId;
  e.put(std::move(TauSpinnerMother), "TauSpinnerMother");

  // Polarimetric vector for Tau from W -> tau + nu
  std::unique_ptr<std::vector<double>> PolarimetricVector(new std::vector<double>);
  *PolarimetricVector = HHvector;
  e.put(std::move(PolarimetricVector), "PolarimetricVector");

  std::unique_ptr<bool> PhotonEmisson1(new bool);
  *PhotonEmisson1 = Photon1;
  e.put(std::move(PhotonEmisson1), "PhotonEmisson1");

  std::unique_ptr<bool> PhotonEmisson2(new bool);
  *PhotonEmisson2 = Photon2;
  e.put(std::move(PhotonEmisson2), "PhotonEmisson2");

  // Type of polariztion polSM
  std::unique_ptr<double> TauPolarisation(new double);
  *TauPolarisation = polSM;
  e.put(std::move(TauPolarisation), "TauPolarisation");

  // regular weight
  std::unique_ptr<double> TauSpinnerWeight(new double);
  *TauSpinnerWeight = WT;
  e.put(std::move(TauSpinnerWeight), "TauSpinnerWT");

  // flipped weight (ie Z->H or H->Z)
  std::unique_ptr<double> TauSpinnerWeightFlip(new double);
  *TauSpinnerWeightFlip = WTFlip;
  e.put(std::move(TauSpinnerWeightFlip), "TauSpinnerWTFlip");

  // h+ polarization
  double WThplus = WT;
  if (polSM < 0.0 && polSM != -999 && isValid)
    WThplus = 0;
  std::unique_ptr<double> TauSpinnerWeighthplus(new double);
  *TauSpinnerWeighthplus = WThplus;
  e.put(std::move(TauSpinnerWeighthplus), "TauSpinnerWThplus");

  // h- polarization
  double WThminus = WT;
  if (polSM > 0.0 && polSM != -999 && isValid)
    WThminus = 0;
  std::unique_ptr<double> TauSpinnerWeighthminus(new double);
  *TauSpinnerWeighthminus = WThminus;
  e.put(std::move(TauSpinnerWeighthminus), "TauSpinnerWThminus");
  return;
}

int TauSpinnerCMS::readParticlesfromReco(edm::Event &e,
                                         SimpleParticle &X,
                                         SimpleParticle &tau,
                                         SimpleParticle &tau2,
                                         std::vector<SimpleParticle> &tau_daughters,
                                         std::vector<SimpleParticle> &tau2_daughters) {
  if (monitoring) std::cout << "start readParticlesfromReco method" << std::endl;
  edm::Handle<reco::GenParticleCollection> genParticles;
  e.getByToken(GenParticleCollectionToken_, genParticles);
  // Loop over genParticles
  for (reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr != genParticles->end(); ++itr) {
    int pdgid = abs(itr->pdgId());
    //std::cout << "PDGID of genparticle = " << pdgid << std::endl;
    // 24      - W+
    // 9900024 - W_R
    // 37      - H+
    // 25      - H0
    // 22      - gamma
    // 23      - Z0
    // genParticles must be W, H, Z, gamma, etc
    // replace initial list of mother particles to avoid gamma checking
    //if (pdgid == 24 || pdgid == 37 || pdgid == 25 || pdgid == 36 || pdgid == 22 || pdgid == 23 || pdgid == 9900024) {
    if (pdgid == MotherPDGID_) {
      if (monitoring) std::cout << "Found mother particle with PDGId = " << pdgid << std::endl;
      const reco::GenParticle *hx = &(*itr);
      if (monitoring) std::cout << "Cheking initial mother daughters list (PDGId = " << itr->pdgId() << ")" << std::endl;
      for (unsigned int i = 0; i < itr->numberOfDaughters(); i++) {
	const reco::GenParticle *dauinit = static_cast<const reco::GenParticle *>(itr->daughter(i));
        if (monitoring) std::cout << "Daughter " << i << " pdgid = " << dauinit->pdgId() << std::endl;
      }
      const reco::GenParticle *LastMotherParticle;
      // If hx is same with it's mother particles returns false
      if (!isFirst(hx))
        continue;
      // recurrent function checks if PDGID of daughter of particle hx is the same with hx's PDGID
      // the replace argument with daughter particle and check it and so on

      GetLastSelfNew(hx, LastMotherParticle);
      if (monitoring) std::cout << "Cheking last mother daughters list (PDGId = " << LastMotherParticle->pdgId() << ")" << std::endl;
      for (unsigned int i = 0; i < LastMotherParticle->numberOfDaughters(); i++) {
	const reco::GenParticle *daulast = static_cast<const reco::GenParticle *>(LastMotherParticle->daughter(i));
        if (monitoring) std::cout << "Daughter " << i << " pdgid = " << daulast->pdgId() << std::endl;
      }
      // Now check the value LastMotherParticle and it's daughter particles
      /*
      if (pdgid == 24) {
        std::cout << "Last Mother PDGID   = " << LastMotherParticle->pdgId() << endl;
        for (unsigned int i = 0; i < LastMotherParticle->numberOfDaughters(); i++) {
      	  const reco::GenParticle *dau = static_cast<const reco::GenParticle *>(LastMotherParticle->daughter(i));
          std::cout << "Daughter number " << i << " PDGID = " << dau->pdgId() << endl;
        }
      }
      */

      const reco::GenParticle *recotau1 = nullptr;
      const reco::GenParticle *recotau2 = nullptr;
      const reco::GenParticle *LastRecotau1 = nullptr;
      const reco::GenParticle *LastRecotau2 = nullptr;
      bool Photon1_found = false;
      bool Photon2_found = false;
      unsigned int ntau(0), ntauornu(0);

      for (unsigned int i = 0; i < LastMotherParticle->numberOfDaughters(); i++) {
        const reco::Candidate *dau = LastMotherParticle->daughter(i);
        if (abs(dau->pdgId()) != pdgid) {
	  // Add Tau_nu_R pdgID = 9900016
          if (abs(dau->pdgId()) == 15 || abs(dau->pdgId()) == 16 || abs(dau->pdgId()) == 9900016) {
            if (ntau == 0 && abs(dau->pdgId()) == 15) {
              recotau1 = static_cast<const reco::GenParticle *>(dau);
              if(monitoring) {
                std::cout << "recotau1 pdgId = " << recotau1->pdgId() << endl;
                std::cout << "Simple check of recotau1 daughters" << endl;
              }
              for (unsigned int j = 0; j < recotau1->numberOfDaughters(); j++) {
                const reco::GenParticle *dau1 = static_cast<const reco::GenParticle *>(recotau1->daughter(j));
		if (abs(dau1->pdgId()) == 22) {
		  Photon1_found = true;
		  for (unsigned int k = 0; k < recotau1->numberOfDaughters(); k++) {
		    const reco::GenParticle *dau11 = static_cast<const reco::GenParticle *>(recotau1->daughter(k));
		    if (dau11->pdgId() == recotau1->pdgId()) LastRecotau1 = dau11;
		  }
		}
                if(monitoring) std::cout << "Daughter of recotau1 number " << j << "   " << dau1->pdgId() << endl;
              }

              //GetLastSelf(recotau1);
              GetLastSelfNew(recotau1, LastRecotau1);
              ntau++;
            } else if ((ntau == 1 && abs(dau->pdgId()) == 15) || abs(dau->pdgId()) == 16 || abs(dau->pdgId()) == 9900016) {
              recotau2 = static_cast<const reco::GenParticle *>(dau);
	      // recotau2 was const
	      // recotau2 = const_cast<reco::GenParticle *>(static_cast<const reco::GenParticle *>(dau));
	      //if (abs(recotau2->pdgId()) == 9900016) {
	      //  recotau2->setPdgId(SignumPDG(recotau2->pdgId()) * 16);
	      //}
              if(monitoring) {
                std::cout << "recotau2 pdgId = " << recotau2->pdgId() << endl;
                std::cout << "Simple check of recotau2 daughters" << endl;
		if (recotau2->numberOfDaughters() == 0) std::cout << "recotau2 has no daughters" << std::endl;
              }
              for (unsigned int j = 0; j < recotau2->numberOfDaughters(); j++) {
                const reco::GenParticle *dau2 = static_cast<const reco::GenParticle *>(recotau2->daughter(j));
		if (abs(dau2->pdgId()) == 22) {
		  Photon1_found = true;
		  for (unsigned int k = 0; k < recotau2->numberOfDaughters(); k++) {
		    const reco::GenParticle *dau22 = static_cast<const reco::GenParticle *>(recotau2->daughter(k));
		    if (dau22->pdgId() == recotau2->pdgId()) LastRecotau2 = dau22;
		  }
		}
                //reco::GenParticle *dau2 = static_cast<reco::GenParticle *>(recotau2->daughter(j));
                if(monitoring) std::cout << "Daughter of recotau2 number " << j << "   " << dau2->pdgId() << endl;
              }
              if (abs(dau->pdgId()) == 15) {
                ntau++;
                //GetLastSelf(recotau2);
                GetLastSelfNew(recotau2, LastRecotau2);
              }
            }
            ntauornu++;
            if(monitoring) {
              std::cout << "recotau1 pdgId = " << recotau1->pdgId() << endl;
              std::cout << "First generation daughters of recotau1" << endl;
	      for (unsigned int j = 0; j < recotau1->numberOfDaughters(); j++) {
                std::cout << "Daughter of recotau1 number " << j << "   " << recotau1->daughter(j)->pdgId() << endl;
	      }
	      /*
              std::cout << "recotau2 pdgId = " << recotau2->pdgId() << endl;
              std::cout << "First generation daughters of recotau2" << endl;
	      for (unsigned int j = 0; j < recotau2->numberOfDaughters(); j++) {
                std::cout << "Daughter of recotau2 number " << j << "   " << recotau2->daughter(j)->pdgId() << endl;
	      }
	      */
            }
          }
        }
      }
      TauCounter = TauCounter + ntau;

      if ((ntau == 2 && ntauornu == 2) || (ntau == 1 && ntauornu == 2)) {
        X.setPx(itr->p4().Px());
        X.setPy(itr->p4().Py());
        X.setPz(itr->p4().Pz());
        X.setE(itr->p4().E());
        X.setPdgid(itr->pdgId());
	if (!Photon1_found) {
          tau.setPx(recotau1->p4().Px());
          tau.setPy(recotau1->p4().Py());
          tau.setPz(recotau1->p4().Pz());
          tau.setE(recotau1->p4().E());
          tau.setPdgid(recotau1->pdgId());
          if (monitoring) std::cout << "GetRecoDaughters method for recotau1" << endl;
          GetRecoDaughters(recotau1, tau_daughters, recotau1->pdgId());
	} else {
          tau.setPx(LastRecotau1->p4().Px());
          tau.setPy(LastRecotau1->p4().Py());
          tau.setPz(LastRecotau1->p4().Pz());
          tau.setE(LastRecotau1->p4().E());
          tau.setPdgid(LastRecotau1->pdgId());
          if (monitoring) std::cout << "GetRecoDaughters method for LastRecotau1" << endl;
          GetRecoDaughters(LastRecotau1, tau_daughters, LastRecotau1->pdgId());
	}
	/*
        if (monitoring) {
          for (unsigned int l = 0; l < tau_daughters.size(); l++) {
            std::cout << "tau_daughters vector element " << l << " with pdgID = " << tau_daughters[l].pdgid() << endl; 
          }
	}
        */
	if (!Photon2_found) {
          tau2.setPx(recotau2->p4().Px());
          tau2.setPy(recotau2->p4().Py());
          tau2.setPz(recotau2->p4().Pz());
          tau2.setE(recotau2->p4().E());
	  // Change PDGId from nu_tau_R to nu_tau
          tau2.setPdgid(SignumPDG(recotau2->pdgId()) * 16);
          if (ntau == 2) {
	    if (monitoring) std::cout << "GetRecoDaughters method for recotau2" << endl;
            GetRecoDaughters(recotau2, tau2_daughters, SignumPDG(recotau2->pdgId()) * 16);
	  }
	} else {
          tau2.setPx(LastRecotau2->p4().Px());
          tau2.setPy(LastRecotau2->p4().Py());
          tau2.setPz(LastRecotau2->p4().Pz());
          tau2.setE(LastRecotau2->p4().E());
	  // Change PDGId from nu_tau_R to nu_tau
          tau2.setPdgid(SignumPDG(LastRecotau2->pdgId()) * 16);
          if (ntau == 2) {
	    if (monitoring) std::cout << "GetRecoDaughters method for LastRecotau2" << endl;
            GetRecoDaughters(LastRecotau2, tau2_daughters, SignumPDG(LastRecotau2->pdgId()) * 16);
	  }
	}
	/*
        if (monitoring) {
          for (unsigned int l = 0; l < tau2_daughters.size(); l++) {
            std::cout << "tau2_daughters vector element " << l << " with pdgID = " << tau2_daughters[l].pdgid() << endl; 
          }
	}
	*/
	if (Photon1_found || Photon2_found) {
	  if (Photon1_found && Photon2_found) return 23;
	  if (Photon1_found) return 2;
	  if (Photon2_found) return 3;
	}
        else return 0;
      }
      //std::cout << "Number of tau or nu_tau = " << ntauornu << endl;
    }
  }
  return 1;
}

// Get pointer to constant value GenParticle
// 
void TauSpinnerCMS::GetLastSelf(const reco::GenParticle *Particle) {
  
  if (Particle->pdgId() == 15) {
    std::cout << "GetLastSelf for tau" << endl;
  }
  if (Particle->pdgId() == 16) {
    std::cout << "GetLastSelf for tau_nu" << endl;
  }
  
  for (unsigned int i = 0; i < Particle->numberOfDaughters(); i++) {
    const reco::GenParticle *dau = static_cast<const reco::GenParticle *>(Particle->daughter(i));
    
    if (Particle->pdgId() == 15 || Particle->pdgId() == 16) {
      std::cout << "Daughter of tau (tau_nu) number " << i << "   " << dau->pdgId() << endl;
    }
    
    if (Particle->pdgId() == dau->pdgId()) {
      
      if (Particle->pdgId() == 15 || Particle->pdgId() == 16) {
        std::cout << "   |   " << endl;
      }
      
      Particle = dau;
      GetLastSelf(Particle);
    }
  }
}

void TauSpinnerCMS::GetLastSelfNew(const reco::GenParticle *Particle, const reco::GenParticle* &LastMother) {
  for (unsigned int i = 0; i < Particle->numberOfDaughters(); i++) {
    const reco::GenParticle *dau = static_cast<const reco::GenParticle *>(Particle->daughter(i));
    if (Particle->pdgId() == dau->pdgId()) {
      Particle = dau;
      if (monitoring) std::cout << "daughter and mother are same" << endl;
      GetLastSelfNew(Particle, LastMother);
    } else {
      LastMother = static_cast<const reco::GenParticle *>(Particle);
      if (monitoring) std::cout << "Last Mother PDGID in function = " << LastMother->pdgId() << endl;
      //LastMother = const_cast<reco::GenParticle *>(static_cast<const reco::GenParticle *>(Particle));
    }
  }
}

int TauSpinnerCMS::SignumPDG(int particlePDGId) {
  int Signum = 1;
  if (particlePDGId < 0) {
    Signum = -1;
  } else if (particlePDGId > 0) {
	Signum = 1;
  }
  return Signum;
}

// Get pointer to constant value GenParticle
bool TauSpinnerCMS::isFirst(const reco::GenParticle *Particle) {
  if (monitoring) {
    std::cout << "isFirst method called for particle with PDGId = " << Particle->pdgId() << std::endl;
    std::cout << "Number of mothers for this particle = " << Particle->numberOfMothers() << std::endl;
  }
  for (unsigned int i = 0; i < Particle->numberOfMothers(); i++) {
    const reco::GenParticle *moth = static_cast<const reco::GenParticle *>(Particle->mother(i));
    if (monitoring) std::cout << "Mother number " << i << " PDGId = " << moth->pdgId() << std::endl;
    if (Particle->pdgId() == moth->pdgId()) {
      return false;
    }
  }
  return true;
}

// Fills the vector of tau's daughters particles
// If firs order daughters of tau have their own daughter we look at last ones
// and so on untill we algorithm have reach particles without daughters or Pi_0
void TauSpinnerCMS::GetRecoDaughters(const reco::GenParticle *Particle,
                                     std::vector<SimpleParticle> &daughters,
                                     int parentpdgid) {
  // replace PDGId of tau_nu_R to PDGId of tau_nu
  // When algo found particle without daughters or pi+, it returning this particle
  if (monitoring) std::cout << "Mother PDGId = " << Particle->pdgId() << std::endl;
  if (Particle->numberOfDaughters() == 0 || abs(Particle->pdgId()) == 111) {
    if (monitoring) std::cout << "Pushed to daughters vector" << std::endl;
    SimpleParticle tp(
        Particle->p4().Px(), Particle->p4().Py(), Particle->p4().Pz(), Particle->p4().E(), Particle->pdgId());
    daughters.push_back(tp);
    return;
  }
  for (unsigned int i = 0; i < Particle->numberOfDaughters(); i++) {
    const reco::Candidate *dau = Particle->daughter(i);
    if (monitoring) std::cout << "daughter PDGId = " <<  dau->pdgId() << std::endl;
    GetRecoDaughters(static_cast<const reco::GenParticle *>(dau), daughters, Particle->pdgId());
  }
}

double TauSpinnerCMS::flat() {
  if (!fRandomEngine) {
    throw cms::Exception("LogicError")
        << "TauSpinnerCMS::flat: Attempt to generate random number when engine pointer is null\n"
        << "This might mean that the code was modified to generate a random number outside the\n"
        << "event and beginLuminosityBlock methods, which is not allowed.\n";
  }
  return fRandomEngine->flat();
}

DEFINE_FWK_MODULE(TauSpinnerCMS);
