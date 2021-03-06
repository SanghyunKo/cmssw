//Package:    EgammaHFProdcers
// Class  :    HFEMClusterProducer
// Original Author:  Kevin Klapoetke (minnesota)
//
// $Id: HFEMClusterProducer.cc,v 1.2 2007/09/19 Kevin Klapoetke
//

#include <iostream>
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoEgamma/EgammaHFProducers/plugins/HFEMClusterProducer.h"
using namespace reco;
HFEMClusterProducer::HFEMClusterProducer(edm::ParameterSet const& conf)
    : hfreco_(consumes<HFRecHitCollection>(conf.getParameter<edm::InputTag>("hits"))) {
  produces<reco::HFEMClusterShapeCollection>();
  produces<reco::BasicClusterCollection>();
  produces<reco::SuperClusterCollection>();
  produces<reco::HFEMClusterShapeAssociationCollection>();
  algo_.setup(conf.getParameter<double>("minTowerEnergy"),
              conf.getParameter<double>("seedThresholdET"),
              conf.getParameter<double>("maximumSL"),
              conf.getParameter<double>("maximumRenergy"),
              conf.getParameter<bool>("usePMTFlag"),
              conf.getParameter<bool>("usePulseFlag"),
              conf.getParameter<bool>("forcePulseFlagMC"),
              conf.getParameter<int>("correctionType"));
}

void HFEMClusterProducer::produce(edm::Event& e, edm::EventSetup const& iSetup) {
  edm::Handle<HFRecHitCollection> hf_hits;

  e.getByToken(hfreco_, hf_hits);

  edm::ESHandle<CaloGeometry> geometry;
  iSetup.get<CaloGeometryRecord>().get(geometry);

  // create return data
  auto retdata1 = std::make_unique<HFEMClusterShapeCollection>();
  auto retdata2 = std::make_unique<SuperClusterCollection>();

  algo_.isMC(!e.isRealData());

  algo_.clusterize(*hf_hits, *geometry, *retdata1, *retdata2);
  edm::OrphanHandle<reco::SuperClusterCollection> SupHandle;
  edm::OrphanHandle<reco::HFEMClusterShapeCollection> ShapeHandle;

  // put the results
  ShapeHandle = e.put(std::move(retdata1));
  SupHandle = e.put(std::move(retdata2));

  auto retdata3 = std::make_unique<HFEMClusterShapeAssociationCollection>(SupHandle, ShapeHandle);

  for (unsigned int i = 0; i < ShapeHandle->size(); i++) {
    retdata3->insert(edm::Ref<reco::SuperClusterCollection>(SupHandle, i),
                     edm::Ref<reco::HFEMClusterShapeCollection>(ShapeHandle, i));
  }

  e.put(std::move(retdata3));
}
