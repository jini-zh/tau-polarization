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
    BJetCandidate(const pat::Jet &jet, bool IsFromPV)
        : Pt(jet.pt()), E(jet.energy()), Eta(jet.eta()), Phi(jet.phi()), FourMomentum(jet.p4()),
        Bprob(jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb")), FromPV(IsFromPV)
    {
    }
    void Monitoring();
    virtual ~BJetCandidate();

    //--- functions

};

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// Constructor of BJetCandidate object
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
BJetCandidate::BJetCandidate() {

};

BJetCandidate::~BJetCandidate()
{
// Destructor
};

void BJetCandidate::Monitoring() {

    std::cout << "pt     = " << Pt << std::endl
    << "eta    = " << Eta << std::endl
    << "phi    = " << Phi << std::endl
    << "BProb  = " << Bprob << std::endl
    << "FromPV = " << FromPV << std::endl;

};

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
};