#ifndef EgammaIsolationProducers_ModifiedRecHitIsolation_h
#define EgammaIsolationProducers_ModifiedRecHitIsolation_h
//*****************************************************************************
// File:      ModifiedRecHitIsolation.h
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer, adapted from EgammaHcalIsolation by S. Harper
// Institute: IIHE-VUB, RAL
//=============================================================================
//*****************************************************************************

//C++ includes
#include <vector>
#include <functional>

//CMSSW includes
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

class ModifiedRecHitIsolation {
 public:

  //constructors
  ModifiedRecHitIsolation (double extRadius,
                         double intRadius,
                         double etaSlice,
                         double etLow,
                         double eLow,
                         edm::ESHandle<CaloGeometry> ,
                         const EcalRecHitCollection&,
                         const EcalSeverityLevelAlgo*,
                         DetId::Detector detector,
                         std::vector<int> recHitFlags,
                         std::vector<int> recHitSeverity,
                         double etaSlice2nd,
                         double intRadius2nd);

  double getEtSum(const reco::GsfElectron * emObject, const reco::TrackBase& addTrk) const {return getSum_(emObject, addTrk, true);}
  double getEnergySum(const reco::GsfElectron * emObject, const reco::TrackBase& addTrk) const{ return  getSum_(emObject, addTrk, false);}

  double getEtSum(const reco::SuperCluster* emObject, const reco::TrackBase& addTrk) const {return getSum_(emObject, addTrk, true);}
  double getEnergySum(const reco::SuperCluster * emObject, const reco::TrackBase& addTrk) const{ return  getSum_(emObject, addTrk, false);}

  void setUseNumCrystals(bool b=true) { useNumCrystals_ = b; }
  void setVetoClustered(bool b=true) { vetoClustered_ = b; }
  void doSeverityChecks(const EcalRecHitCollection *const recHits,
			const std::vector<int>& v) {
    ecalBarHits_ = recHits;
    severitiesexcl_.clear();
    severitiesexcl_.insert(severitiesexcl_.begin(), v.begin(), v.end());
    std::sort(severitiesexcl_.begin(), severitiesexcl_.end());
  }

  void doFlagChecks(const std::vector<int>& v) {
    flags_.clear();
    flags_.insert(flags_.begin(), v.begin(), v.end());
    std::sort(flags_.begin(), flags_.end() );
  }

  //destructor
  ~ModifiedRecHitIsolation() ;

 private:
  double getSum_(const reco::GsfElectron *, const reco::TrackBase& addTrk, bool returnEt ) const;
  double getSum_(const reco::SuperCluster *, const reco::TrackBase& addTrk, bool returnEt ) const;

  double extRadius_;
  double intRadius_;
  double etaSlice_;
  double etLow_;
  double eLow_;
  double etaSlice2nd_;
  double intRadius2nd_;

  edm::ESHandle<CaloGeometry>  theCaloGeom_;
  const EcalRecHitCollection&  caloHits_;
  const EcalSeverityLevelAlgo* sevLevel_;

  bool useNumCrystals_;
  bool vetoClustered_;
  const EcalRecHitCollection *ecalBarHits_;
  std::vector<int> severitiesexcl_;
  std::vector<int> flags_;

  const CaloSubdetectorGeometry* subdet_[2]; // barrel+endcap
};

#endif
