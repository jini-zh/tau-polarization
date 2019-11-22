
#include <memory>
#include <string>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
//
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"


#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

//
// class decleration
//

class MySimpleParticle {
  public:
    MySimpleParticle() { _px=_py=_pz=_e=0; _pdgid=0; }
    MySimpleParticle(double x, double y, double z, double e, int id) { _px=x; _py=y; _pz=z; _e=e; _pdgid=id; }
  
    double px()    { return _px;    }
    double py()    { return _py;    }
    double pz()    { return _pz;    }
    double e ()    { return _e;     }
    int    pdgid() { return _pdgid; }
  
    void setPx(double x ) { _px = x;     }
    void setPy(double y ) { _py = y;     }
    void setPz(double z ) { _pz = z;     }
    void setE (double e ) { _e  = e;     }
    void setPdgid(int id) { _pdgid = id; }
 
  private:
    double _px,_py,_pz,_e;
    int    _pdgid;
};
