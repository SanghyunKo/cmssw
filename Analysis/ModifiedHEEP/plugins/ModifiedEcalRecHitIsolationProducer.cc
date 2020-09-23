//*****************************************************************************
// File:      ModifiedEcalRecHitIsolationProducer.cc
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer
// Institute: IIHE-VUB
//=============================================================================
//*****************************************************************************


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Analysis/ModifiedHEEP/interface/ModifiedRecHitIsolation.h"

// Framework
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"

#include "CommonTools/Utils/interface/StringToEnumValue.h"

#include "Analysis/ModifiedHEEP/interface/ModifiedEleTkIsolFromCands.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

class ModifiedEcalRecHitIsolationProducer : public edm::stream::EDProducer<> {
public:
  explicit ModifiedEcalRecHitIsolationProducer(const edm::ParameterSet&);
  ~ModifiedEcalRecHitIsolationProducer() override;


  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  template <typename T> void setToken(edm::EDGetTokenT<T>& token,edm::InputTag tag){token=consumes<T>(tag);}
  template <typename T> void setToken(edm::EDGetTokenT<T>& token,const edm::ParameterSet& iPara,const std::string& tag){token=consumes<T>(iPara.getParameter<edm::InputTag>(tag));}
  template <typename T> void setToken(std::vector<edm::EDGetTokenT<T> >& tokens,const edm::ParameterSet& iPara,const std::string& tagName) {
    auto tags =iPara.getParameter<std::vector<edm::InputTag> >(tagName);
    for(auto& tag : tags) {
      edm::EDGetTokenT<T> token;
      setToken(token,tag);
      tokens.push_back(token);
    }
  }
  // ----------member data ---------------------------

  edm::EDGetTokenT<edm::View<reco::GsfElectron>> emObjectToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalBarrelRecHitToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalEndcapRecHitToken_;
  edm::EDGetTokenT<reco::GsfTrackCollection> gsfTrkToken_;

  double egIsoPtMinBarrel_; //minimum Et noise cut
  double egIsoEMinBarrel_;  //minimum E noise cut
  double egIsoPtMinEndcap_; //minimum Et noise cut
  double egIsoEMinEndcap_;  //minimum E noise cut
  double egIsoConeSizeOut_; //outer cone size
  double egIsoConeSizeInBarrel_; //inner cone size
  double egIsoConeSizeInEndcap_; //inner cone size
  double egIsoJurassicWidth_ ; // exclusion strip width for jurassic veto

  double egIsoJurassicWidth2nd_;
  double egIsoConeSizeIn2nd_;

  bool useIsolEt_; //switch for isolEt rather than isolE
  bool tryBoth_ ; // use rechits from barrel + endcap
  bool subtract_ ; // subtract SC energy (allows veto cone of zero size)

  bool useNumCrystals_ ; // veto on number of crystals
  bool vetoClustered_ ;  // veto all clusterd rechits

  std::vector<std::string> recHitFlagsEB_;
  std::vector<std::string> recHitFlagsEE_;
  std::vector<std::string> recHitSeverityEB_;
  std::vector<std::string> recHitSeverityEE_;

  std::vector<int> recHitFlagsEnumsEB_;
  std::vector<int> recHitFlagsEnumsEE_;
  std::vector<int> recHitSeverityEnumsEB_;
  std::vector<int> recHitSeverityEnumsEE_;

  ModifiedEleTkIsolFromCands trkIsoCalc_;

  edm::ParameterSet conf_;

};

ModifiedEcalRecHitIsolationProducer::ModifiedEcalRecHitIsolationProducer(const edm::ParameterSet& config):
  trkIsoCalc_(config.getParameter<edm::ParameterSet>("trkIsoConfig")),
  conf_(config)
{
  setToken(emObjectToken_,config,"emObjectProducer");
  setToken(ecalBarrelRecHitToken_,config,"ecalBarrelRecHitCollection");
  setToken(ecalEndcapRecHitToken_,config,"ecalEndcapRecHitCollection");
  setToken(gsfTrkToken_,config,"gsfTrks");
  // use configuration file to setup input/output collection names
  //inputs
  // emObjectProducer_               = conf_.getParameter<edm::InputTag>("emObjectProducer");
  // ecalBarrelRecHitProducer_       = conf_.getParameter<edm::InputTag>("ecalBarrelRecHitProducer");
  // ecalBarrelRecHitCollection_     = conf_.getParameter<edm::InputTag>("ecalBarrelRecHitCollection");
  // ecalEndcapRecHitProducer_       = conf_.getParameter<edm::InputTag>("ecalEndcapRecHitProducer");
  // ecalEndcapRecHitCollection_     = conf_.getParameter<edm::InputTag>("ecalEndcapRecHitCollection");

  //vetos
  egIsoPtMinBarrel_               = conf_.getParameter<double>("etMinBarrel");
  egIsoEMinBarrel_                = conf_.getParameter<double>("eMinBarrel");
  egIsoPtMinEndcap_               = conf_.getParameter<double>("etMinEndcap");
  egIsoEMinEndcap_                = conf_.getParameter<double>("eMinEndcap");
  egIsoConeSizeInBarrel_          = conf_.getParameter<double>("intRadiusBarrel");
  egIsoConeSizeInEndcap_          = conf_.getParameter<double>("intRadiusEndcap");
  egIsoConeSizeOut_               = conf_.getParameter<double>("extRadius");
  egIsoJurassicWidth_             = conf_.getParameter<double>("jurassicWidth");
  egIsoJurassicWidth2nd_          = conf_.getParameter<double>("jurassicWidth2nd");
  egIsoConeSizeIn2nd_             = conf_.getParameter<double>("intRadius2nd");

  // options
  useIsolEt_      = conf_.getParameter<bool>("useIsolEt");
  tryBoth_        = conf_.getParameter<bool>("tryBoth");
  subtract_       = conf_.getParameter<bool>("subtract");
  useNumCrystals_ = conf_.getParameter<bool>("useNumCrystals");
  vetoClustered_  = conf_.getParameter<bool>("vetoClustered");

  recHitFlagsEB_ = conf_.getParameter<std::vector<std::string>>("recHitFlagsExclBarrel");
  recHitFlagsEE_ = conf_.getParameter<std::vector<std::string>>("recHitFlagsExclEndcaps");
  recHitSeverityEB_ = conf_.getParameter<std::vector<std::string>>("recHitSeverityExclBarrel");
  recHitSeverityEE_ = conf_.getParameter<std::vector<std::string>>("recHitSeverityExclEndcaps");

  recHitFlagsEnumsEB_ = StringToEnumValue<EcalRecHit::Flags>(recHitFlagsEB_);
  recHitFlagsEnumsEE_ = StringToEnumValue<EcalRecHit::Flags>(recHitFlagsEE_);
  recHitSeverityEnumsEB_ = StringToEnumValue<EcalSeverityLevel::SeverityLevel>(recHitSeverityEB_);
  recHitSeverityEnumsEE_ = StringToEnumValue<EcalSeverityLevel::SeverityLevel>(recHitSeverityEE_);

  //register your products
  produces < edm::ValueMap<double> >("EcalRecHitIso");
}


ModifiedEcalRecHitIsolationProducer::~ModifiedEcalRecHitIsolationProducer(){}

//
// member functions
//

// ------------ method called to produce the data  ------------
void
ModifiedEcalRecHitIsolationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Get the  filtered objects
  edm::Handle< edm::View<reco::GsfElectron> > emObjectHandle;
  iEvent.getByToken(emObjectToken_,emObjectHandle);

  // Next get Ecal hits barrel
  edm::Handle<EcalRecHitCollection> ecalBarrelRecHitHandle; //EcalRecHitCollection is a typedef to
  iEvent.getByToken(ecalBarrelRecHitToken_, ecalBarrelRecHitHandle);

  // Next get Ecal hits endcap
  edm::Handle<EcalRecHitCollection> ecalEndcapRecHitHandle;
  iEvent.getByToken(ecalEndcapRecHitToken_, ecalEndcapRecHitHandle);

  edm::Handle<reco::GsfTrackCollection> gsfTrkHandle;
  iEvent.getByToken(gsfTrkToken_, gsfTrkHandle);

  edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
  iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);
  const EcalSeverityLevelAlgo* sevLevel = sevlv.product();

  //Get Calo Geometry
  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* caloGeom = pG.product();

  //reco::CandViewDoubleAssociations* isoMap = new reco::CandViewDoubleAssociations( reco::CandidateBaseRefProd( emObjectHandle ) );
  auto isoMap = std::make_unique<edm::ValueMap<double>>();
  edm::ValueMap<double>::Filler filler(*isoMap);
  std::vector<double> retV(emObjectHandle->size(),0);

  ModifiedRecHitIsolation ecalBarrelIsol(egIsoConeSizeOut_,egIsoConeSizeInBarrel_,egIsoJurassicWidth_,egIsoPtMinBarrel_,egIsoEMinBarrel_,caloGeom,*ecalBarrelRecHitHandle,sevLevel,DetId::Ecal,recHitFlagsEnumsEB_,recHitSeverityEnumsEB_,egIsoJurassicWidth2nd_,egIsoConeSizeIn2nd_);
  ecalBarrelIsol.setUseNumCrystals(useNumCrystals_);
  ecalBarrelIsol.setVetoClustered(vetoClustered_);

  ModifiedRecHitIsolation ecalEndcapIsol(egIsoConeSizeOut_,egIsoConeSizeInEndcap_,egIsoJurassicWidth_,egIsoPtMinEndcap_,egIsoEMinEndcap_,caloGeom,*ecalEndcapRecHitHandle,sevLevel,DetId::Ecal,recHitFlagsEnumsEE_,recHitSeverityEnumsEE_,egIsoJurassicWidth2nd_,egIsoConeSizeIn2nd_);
  ecalEndcapIsol.setUseNumCrystals(useNumCrystals_);
  ecalEndcapIsol.setVetoClustered(vetoClustered_);


  for( size_t i = 0 ; i < emObjectHandle->size(); ++i) {

    //i need to know if its in the barrel/endcap so I get the supercluster handle to find out the detector eta
    //this might not be the best way, are we guaranteed that eta<1.5 is barrel
    //this can be safely replaced by another method which determines where the emobject is
    //then we either get the isolation Et or isolation Energy depending on user selection
    double isoValue =0.;

    reco::SuperClusterRef superClus = emObjectHandle->at(i).get<reco::SuperClusterRef>();

    bool addGsfTrkSel = false;

    auto additionalGsfTrk = trkIsoCalc_.additionalGsfTrkSelector(*emObjectHandle->ptrAt(i),gsfTrkHandle, addGsfTrkSel);

    if(tryBoth_){ //barrel + endcap
      if(useIsolEt_) isoValue =  ecalBarrelIsol.getEtSum(&(emObjectHandle->at(i)),*(additionalGsfTrk.get())) + ecalEndcapIsol.getEtSum(&(emObjectHandle->at(i)),*(additionalGsfTrk.get()));
      else           isoValue =  ecalBarrelIsol.getEnergySum(&(emObjectHandle->at(i)),*(additionalGsfTrk.get())) + ecalEndcapIsol.getEnergySum(&(emObjectHandle->at(i)),*(additionalGsfTrk.get()));
    }
    else if( fabs(superClus->eta())<1.479) { //barrel
      if(useIsolEt_) isoValue =  ecalBarrelIsol.getEtSum(&(emObjectHandle->at(i)),*(additionalGsfTrk.get()));
      else           isoValue =  ecalBarrelIsol.getEnergySum(&(emObjectHandle->at(i)),*(additionalGsfTrk.get()));
    }
    else{ //endcap
      if(useIsolEt_) isoValue =  ecalEndcapIsol.getEtSum(&(emObjectHandle->at(i)),*(additionalGsfTrk.get()));
      else           isoValue =  ecalEndcapIsol.getEnergySum(&(emObjectHandle->at(i)),*(additionalGsfTrk.get()));
    }

    //we subtract off the electron energy here as well
    double subtractVal=0;

    if(useIsolEt_) subtractVal = superClus.get()->rawEnergy()*sin(2*atan(exp(-superClus.get()->eta())));
    else           subtractVal = superClus.get()->rawEnergy();

    if(subtract_) isoValue-= subtractVal;

    retV[i]=isoValue;
    //all done, isolation is now in the map

  }//end of loop over em objects

  filler.insert(emObjectHandle,retV.begin(),retV.end());
  filler.fill();

  iEvent.put(std::move(isoMap),"EcalRecHitIso");

}

// define this as a plug-in
DEFINE_FWK_MODULE(ModifiedEcalRecHitIsolationProducer);
