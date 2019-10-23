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
      MotherPDGID_(pset.getUntrackedParameter("MotherPDGID", (int)(-1))), // -1 - all particles
      Ipol_(pset.getUntrackedParameter("Ipol", (int)(0))), // For Z/gamma case only???
      nonSM2_(pset.getUntrackedParameter("nonSM2", (int)(0))),
      nonSMN_(pset.getUntrackedParameter("nonSMN", (int)(0))),
      monitoring(pset.getParameter<bool>("monitoring")),
      roundOff_(pset.getUntrackedParameter("roundOff", (double)(0.01))) {
  //usesResource(edm::uniqueSharedResourceName());
  usesResource(edm::SharedResourceNames::kTauola);

  produces<bool>("TauSpinnerWTisValid").setBranchAlias("TauSpinnerWTisValid");
  produces<double>("TauSpinnerWT").setBranchAlias("TauSpinnerWT");
  produces<double>("TauSpinnerWTFlip").setBranchAlias("TauSpinnerWTFlip");
  produces<double>("TauSpinnerWThplus").setBranchAlias("TauSpinnerWThplus");
  produces<double>("TauSpinnerWThminus").setBranchAlias("TauSpinnerWThminus");
  produces<int>("TauSpinnerMother").setBranchAlias("TauSpinnerMother");

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

  EventNumber++;
  if (monitoring) {
    std::cout << endl << "### Event Number = " << EventNumber << endl;
    std::cout << "isReco_ = " << isReco_ << endl;
  }

  if (isReco_) {
  	// readParticlesfromReco returns 0 if there is decay with tau
  	// returns 1 if there are no taus
    stat = readParticlesfromReco(e, X, tau, tau2, tau_daughters, tau_daughters2);
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
      if (abs(X.pdgid()) == 24 || abs(X.pdgid()) == 37) {
      	if (abs(X.pdgid()) == 24) {
      	  if (monitoring) std::cout << "W+ is mother" << std::endl;
      	} else if (abs(X.pdgid()) == 37) {
      	  if (monitoring) std::cout << "H+ is mother" << std::endl;
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
  edm::Handle<reco::GenParticleCollection> genParticles;
  e.getByToken(GenParticleCollectionToken_, genParticles);
  // Loop over genParticles
  for (reco::GenParticleCollection::const_iterator itr = genParticles->begin(); itr != genParticles->end(); ++itr) {
    int pdgid = abs(itr->pdgId());
    // 24 - W+
    // 37 - H+
    // 25 - H0
    // 22 - gamma
    // 23 - Z0
    // genParticles must be W, H, Z, gamma, etc
    if (pdgid == 24 || pdgid == 37 || pdgid == 25 || pdgid == 36 || pdgid == 22 || pdgid == 23) {
      const reco::GenParticle *hx = &(*itr);
      const reco::GenParticle *LastMotherParticle;
      // If hx is same with it's mother particles returns false
      if (!isFirst(hx))
        continue;
      // recurrent function checks if PDGID of daughter of particle hx is the same with hx's PDGID
      // the replace argument with daughter particle and check it and so on
      // std::cout << "Calling GetLastSelf method for MOTHER particle" << endl;
      //if (pdgid == 24) {
      //  std::cout << "GetLastSelf method for W-boson " << endl;
      //}
      GetLastSelfNew(hx, LastMotherParticle);

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
      unsigned int ntau(0), ntauornu(0);
      /*
      // Loop over itr daughters
      for (unsigned int i = 0; i < itr->numberOfDaughters(); i++) {
        const reco::Candidate *dau = itr->daughter(i);
        if (abs(dau->pdgId()) != pdgid) {
          if (abs(dau->pdgId()) == 15 || abs(dau->pdgId()) == 16) {
            if (ntau == 0 && abs(dau->pdgId()) == 15) {
              std::cout << "Found tau in 1st step" << endl;
              recotau1 = static_cast<const reco::GenParticle *>(dau);
              GetLastSelf(recotau1);
              ntau++;
            } else if ((ntau == 1 && abs(dau->pdgId()) == 15) || abs(dau->pdgId()) == 16) {
              recotau2 = static_cast<const reco::GenParticle *>(dau);
              if (abs(dau->pdgId()) == 15) {
                ntau++;
                GetLastSelf(recotau2);
              }
            }
            ntauornu++;
          }
        }
      }
      */
      for (unsigned int i = 0; i < LastMotherParticle->numberOfDaughters(); i++) {
        const reco::Candidate *dau = LastMotherParticle->daughter(i);
        if (abs(dau->pdgId()) != pdgid) {
          if (abs(dau->pdgId()) == 15 || abs(dau->pdgId()) == 16) {
            if (ntau == 0 && abs(dau->pdgId()) == 15) {
              //std::cout << "Found tau 1!" << endl;
              recotau1 = static_cast<const reco::GenParticle *>(dau);
              if(monitoring) {
                std::cout << "recotau1 pdgId = " << recotau1->pdgId() << endl;
                std::cout << "Simple check of recotau1 daughters" << endl;
              }
              for (unsigned int j = 0; j < recotau1->numberOfDaughters(); j++) {
                const reco::GenParticle *dau1 = static_cast<const reco::GenParticle *>(recotau1->daughter(j));
                if(monitoring) std::cout << "Daughter of recotau1 number " << j << "   " << dau1->pdgId() << endl;
              }

              GetLastSelf(recotau1);
              ntau++;
            } else if ((ntau == 1 && abs(dau->pdgId()) == 15) || abs(dau->pdgId()) == 16) {
              recotau2 = static_cast<const reco::GenParticle *>(dau);
              if(monitoring) {
                std::cout << "recotau2 pdgId = " << recotau1->pdgId() << endl;
                std::cout << "Simple check of recotau1 daughters" << endl;
              }
              for (unsigned int j = 0; j < recotau1->numberOfDaughters(); j++) {
                const reco::GenParticle *dau1 = static_cast<const reco::GenParticle *>(recotau1->daughter(j));
                if(monitoring) std::cout << "Daughter of recotau2 number " << j << "   " << dau1->pdgId() << endl;
              }
              if (abs(dau->pdgId()) == 15) {
                ntau++;
                GetLastSelf(recotau2);
              }
            }
            ntauornu++;
          }
        }
      }
      TauCounter = TauCounter + ntau;
      // Заполняет четырёхвекторы материнской частицы (например, W) и двух дочерних (tau и nu_tau)
      // А также заполняет std::vector<SimpleParticle> четырёхвекторами дочерних по отношению к tau частиц
      if ((ntau == 2 && ntauornu == 2) || (ntau == 1 && ntauornu == 2)) {
        X.setPx(itr->p4().Px());
        X.setPy(itr->p4().Py());
        X.setPz(itr->p4().Pz());
        X.setE(itr->p4().E());
        X.setPdgid(itr->pdgId());
        tau.setPx(recotau1->p4().Px());
        tau.setPy(recotau1->p4().Py());
        tau.setPz(recotau1->p4().Pz());
        tau.setE(recotau1->p4().E());
        tau.setPdgid(recotau1->pdgId());
        // std::cout << "GetRecoDaughters method for recotau1" << endl;
        GetRecoDaughters(recotau1, tau_daughters, recotau1->pdgId());
        /*
        for (unsigned int l = 0; l < tau_daughters.size(); l++) {
          std::cout << "tau_daughters vector element " << l << " with pdgID = " << tau_daughters[l].pdgid() << endl; 
        }
        */
        tau2.setPx(recotau2->p4().Px());
        tau2.setPy(recotau2->p4().Py());
        tau2.setPz(recotau2->p4().Pz());
        tau2.setE(recotau2->p4().E());
        tau2.setPdgid(recotau2->pdgId());
        if (ntau == 2)
          GetRecoDaughters(recotau2, tau2_daughters, recotau2->pdgId());
        return 0;
      }
      //std::cout << "Number of tau or nu_tau = " << ntauornu << endl;
    }
  }
  return 1;
}

// Get pointer to constant value GenParticle
// 
void TauSpinnerCMS::GetLastSelf(const reco::GenParticle *Particle) {
  /*
  if (Particle->pdgId() == 15) {
    std::cout << "GetLastSelf for tau" << endl;
  }
  if (Particle->pdgId() == 16) {
    std::cout << "GetLastSelf for tau_nu" << endl;
  }
  */
  for (unsigned int i = 0; i < Particle->numberOfDaughters(); i++) {
    const reco::GenParticle *dau = static_cast<const reco::GenParticle *>(Particle->daughter(i));
    /*
    if (Particle->pdgId() == 15 || Particle->pdgId() == 16) {
      std::cout << "Daughter of tau (tau_nu) number " << i << "   " << dau->pdgId() << endl;
    }
    */
    if (Particle->pdgId() == dau->pdgId()) {
      /*
      if (Particle->pdgId() == 15 || Particle->pdgId() == 16) {
        std::cout << "   |   " << endl;
      }
      */
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
      //std::cout << "daughter and mother are same" << endl;
      GetLastSelfNew(Particle, LastMother);
    } else {
      LastMother = static_cast<const reco::GenParticle *>(Particle);
      //std::cout << "Last Mother PDGID in function = " << LastMother->pdgId() << endl;
      //LastMother = const_cast<reco::GenParticle *>(static_cast<const reco::GenParticle *>(Particle));
    }
  }
}

// Get pointer to constant value GenParticle
bool TauSpinnerCMS::isFirst(const reco::GenParticle *Particle) {
  for (unsigned int i = 0; i < Particle->numberOfMothers(); i++) {
    const reco::GenParticle *moth = static_cast<const reco::GenParticle *>(Particle->mother(i));
    if (Particle->pdgId() == moth->pdgId()) {
      return false;
    }
  }
  return true;
}

// Fills the vector of tau's daughters particles
// If firs order daughters of tau have ther ows daughter we look at last ones
// and so on untill we algorithm have reach particles without daughters or Pi_0
void TauSpinnerCMS::GetRecoDaughters(const reco::GenParticle *Particle,
                                     std::vector<SimpleParticle> &daughters,
                                     int parentpdgid) {
  if (Particle->numberOfDaughters() == 0 || abs(Particle->pdgId()) == 111) {
    SimpleParticle tp(
        Particle->p4().Px(), Particle->p4().Py(), Particle->p4().Pz(), Particle->p4().E(), Particle->pdgId());
    daughters.push_back(tp);
    return;
  }
  for (unsigned int i = 0; i < Particle->numberOfDaughters(); i++) {
    const reco::Candidate *dau = Particle->daughter(i);
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
