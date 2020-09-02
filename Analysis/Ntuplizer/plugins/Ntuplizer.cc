// -*- C++ -*-
//
// Package:    Analysis/Ntuplizer
// Class:      Ntuplizer
//
/**\class Ntuplizer Ntuplizer.cc Analysis/Ntuplizer/plugins/Ntuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mr. Ko Sanghyun
//         Created:  Tue, 10 Apr 2018 07:28:14 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// CMSSW include files w/o mkedanlzr
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// ROOT
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "Analysis/Ntuplizer/interface/Ntuplizer.h"

using namespace std;
using namespace reco;
using namespace edm;
using namespace isodeposit;
using namespace pat;

namespace HEEPV70 {
  enum CutIndex {
    ET=0,ETA,DETAINSEED,DPHIIN,SIGMAIETAIETA,E2X5OVER5X5,HADEM,TRKISO,EMHADD1ISO,DXY,MISSHITS,ECALDRIVEN
  };
}

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Ntuplizer(const edm::ParameterSet&);
      ~Ntuplizer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void fillTriggers(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      virtual void fillMuons(const edm::Event& iEvent);
      virtual void fillElectrons(const edm::Event& iEvent);
      virtual void fillPhotons(const edm::Event& iEvent);
      virtual void fillJets(const edm::Event& iEvent);
      virtual void fillMETs(const edm::Event& iEvent);
      virtual void fillGenParticles(const edm::Event& iEvent);
      virtual void fillGenJets(const edm::Event& iEvent);
      virtual void fillGsfTracks(const edm::Event& iEvent);
      virtual void fillGeneralTracks(const edm::Event& iEvent);

      // ----------member data ---------------------------
      bool isMC;

      edm::EDGetTokenT<edm::TriggerResults> triggerToken;
      edm::EDGetTokenT<trigger::TriggerEvent> triggersummaryToken;
      edm::EDGetTokenT<reco::BeamSpot> beamspotToken;
      edm::EDGetTokenT<reco::VertexCollection> vertexToken;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUToken;
      edm::EDGetTokenT<std::vector<reco::Muon>> muonsToken;
      edm::EDGetTokenT<std::vector<reco::GsfElectron>> electronsToken;
      edm::EDGetTokenT<std::vector<reco::Conversion>> conversionsToken;
      edm::EDGetTokenT<std::vector<reco::Photon>> photonToken;
      edm::EDGetTokenT<std::vector<reco::PFJet>> jetToken;
      edm::EDGetTokenT<std::vector<reco::PFMET>> METToken;
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genparticlesToken;
      edm::EDGetTokenT<std::vector<reco::GenJet>> genjetToken;
      edm::EDGetTokenT<GenEventInfoProduct> generatorToken;
      edm::EDGetTokenT<std::vector<reco::GsfTrack>> GsfTrackToken;
      edm::EDGetTokenT<std::vector<reco::Track>> GeneralTrackToken;
      // edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>> HEEPVIDResultToken;
      edm::EDGetTokenT<edm::ValueMap<int>> nrSatCrysMapToken;
      edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken;
      edm::EDGetTokenT<edm::ValueMap<int>> nrMatchedTrkToken;
      edm::EDGetTokenT<edm::ValueMap<float>> rtMatchedTrkToken;
      edm::EDGetTokenT<edm::ValueMap<bool>> addGsfTrkSelToken;
      edm::EDGetTokenT<edm::ValueMap<double>> EcalRecHitIsoToken;

      edm::EDGetTokenT<double> rhoToken;

      edm::ESHandle<MagneticField> MagField;
      edm::ESHandle<TransientTrackBuilder> TTBuilder;
      edm::Handle<reco::BeamSpot> beamSpotHandle;
      edm::Handle<std::vector<reco::Vertex>> PVHandle;

      reco::BeamSpot beamSpot;
      reco::Vertex vtx;

      std::map<std::string, TTree*> tree_;
      NtupleEvent evt_;
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
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig):
isMC(iConfig.getUntrackedParameter<bool>("isMC")),
triggerToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
triggersummaryToken(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerSummary"))),
beamspotToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
vertexToken(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("Vertex"))),
PUToken(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PU"))),
muonsToken(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("Muons"))),
electronsToken(consumes<std::vector<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("Electrons"))),
conversionsToken(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("Conversions"))),
photonToken(consumes<std::vector<reco::Photon>>(iConfig.getParameter<edm::InputTag>("Photons"))),
jetToken(consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("Jets"))),
METToken(consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("MET"))),
genparticlesToken(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("GenParticles"))),
genjetToken(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("GenJets"))),
generatorToken(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("Generator"))),
GsfTrackToken(consumes<std::vector<reco::GsfTrack>>(iConfig.getParameter<edm::InputTag>("GsfTracks"))),
GeneralTrackToken(consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("GeneralTracks"))),
// HEEPVIDResultToken(consumes<edm::ValueMap<vid::CutFlowResult>>(iConfig.getParameter<edm::InputTag>("HEEPVID"))),
nrSatCrysMapToken(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("nrSatCrysMap"))),
trkIsoMapToken(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
nrMatchedTrkToken(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("nrMatchedTrkMap"))),
rtMatchedTrkToken(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("rtMatchedTrkMap"))),
addGsfTrkSelToken(consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("addGsfTrkSelMap"))),
EcalRecHitIsoToken(consumes<edm::ValueMap<double>>(iConfig.getParameter<edm::InputTag>("EcalRecHitIsoMap"))),

rhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("rho")))

{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


Ntuplizer::~Ntuplizer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  evt_.triggers.clear();
  evt_.muons.clear();
  evt_.electrons.clear();
  evt_.photons.clear();
  evt_.jets.clear();
  evt_.METs.clear();
  evt_.genparticles.clear();
  evt_.genjets.clear();
  evt_.GsfTracks.clear();
  evt_.tracks.clear();

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder);
  iSetup.get<IdealMagneticFieldRecord>().get(MagField);

  iEvent.getByToken(beamspotToken, beamSpotHandle); beamSpot = (*beamSpotHandle);
  iEvent.getByToken(vertexToken, PVHandle); vtx = PVHandle->front();

  evt_.run = iEvent.id().run();
  evt_.event = iEvent.id().event();
  evt_.lumi = iEvent.id().luminosityBlock();
  evt_.nVertices = PVHandle->size();

  if (isMC) {
    edm::Handle<std::vector<PileupSummaryInfo>> PUInfo;
    iEvent.getByToken(PUToken, PUInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVIdx;

    int npv = -1; int npvin = -1;
    for (PVIdx = PUInfo->begin(); PVIdx != PUInfo->end(); ++PVIdx) {
      int BX = PVIdx->getBunchCrossing();
      if (BX == 0) {
        npvin = PVIdx->getPU_NumInteractions();
        npv = PVIdx->getTrueNumInteractions();
        continue;
      }
    }
    evt_.nPU = npv;
    evt_.nPUin = npvin;
  }

  fillTriggers(iEvent, iSetup);
  fillMuons(iEvent);
  fillElectrons(iEvent);
  fillPhotons(iEvent);
  fillJets(iEvent);
  // fillMETs(iEvent);
  if(isMC) {
    fillGenParticles(iEvent);
    fillGenJets(iEvent);
  }
  fillGsfTracks(iEvent);
  fillGeneralTracks(iEvent);
  tree_["PhysicsTree"]->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
Ntuplizer::beginJob()
{
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> FS;
  tree_["PhysicsTree"] = FS->make<TTree>("PhysicsTree","PhysicsTree");
  tree_["PhysicsTree"]->Branch("Event", &evt_);
}

// ------------ method called once each job just after ending the event loop  ------------
void
Ntuplizer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void Ntuplizer::fillTriggers(const edm::Event &iEvent, const edm::EventSetup& iSetup) {
  NtupleTrigger Trig_;

  edm::Handle<TriggerResults> ResultHandle;
  iEvent.getByToken(triggerToken,ResultHandle);

  string Trigs[] = {
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*",
    "HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v*",
    "HLT_Ele27_WPTight_Gsf_v*",
    "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
    "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",
    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
    "HLT_Ele22_eta2p1_WP75_Gsf_v*",
    "HLT_Ele22_eta2p1_WPLoose_Gsf_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  };

  const unsigned int nTrigs = sizeof(Trigs)/sizeof(*Trigs);
  const unsigned int nTrig = ResultHandle.product()->size();
  std::vector<std::pair<std::string, int>> indices;
  edm::TriggerNames TrigList = iEvent.triggerNames(*ResultHandle);

  for (unsigned int iTrig = 0; iTrig != nTrig; iTrig++) {
    string TrigName_ = TrigList.triggerName(iTrig);
    for (unsigned int jTrig = 0; jTrig != nTrigs; jTrig++) {
      if (TrigName_.find(Trigs[jTrig].substr(0, Trigs[jTrig].find("*"))) != std::string::npos) {
        Trig_.name = TrigName_;
        if (ResultHandle.product()->accept(iTrig)) Trig_.isFired = true;
        else Trig_.isFired = false;
        evt_.triggers.push_back(Trig_);
      }
    }
  }
}

void Ntuplizer::fillMuons(const edm::Event& iEvent) {
  NtupleMuon mu_;

  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken,muons);

  evt_.nMuons = muons->size();
  if(muons->size() == 0) return;

  for (std::vector<reco::Muon>::const_iterator mu = muons->begin(); mu != muons->end(); ++mu) {
    mu_.pt = mu->pt();
    mu_.eta = mu->eta();
    mu_.phi = mu->phi();
    mu_.charge = mu->charge();
    mu_.isStandAloneMuon = mu->isStandAloneMuon();
    mu_.isGlobalMuon = mu->isGlobalMuon();
    mu_.isTrackerMuon = mu->isTrackerMuon();
    mu_.isPFMuon = mu->isPFMuon();
    mu_.px = mu->px();
    mu_.py = mu->py();
    mu_.pz = mu->pz();
    mu_.nChambers = mu->numberOfChambers();
    mu_.stationMask = mu->stationMask();
    mu_.nMatchedStations = mu->numberOfMatchedStations();

    if (mu->isGlobalMuon()) {
      reco::TrackRef glbTrack = mu->globalTrack();
      if (glbTrack.isNonnull()) {
        const reco::HitPattern& glbhit = glbTrack->hitPattern();
        mu_.normalizedChi2 = glbTrack->normalizedChi2();
        mu_.nValidHits = glbTrack->numberOfValidHits();
        mu_.nValidMuonHits = glbhit.numberOfValidMuonHits();
        mu_.qoverp = glbTrack->qoverp();
        mu_.theta = glbTrack->theta();
        mu_.lambda = glbTrack->lambda();
        mu_.dxy = glbTrack->dxy();
        mu_.d0 = glbTrack->d0();
        mu_.dsz = glbTrack->dsz();
        mu_.dz = glbTrack->dz();
        mu_.dxyBS = glbTrack->dxy(beamSpot.position());
        mu_.dszBS = glbTrack->dsz(beamSpot.position());
        mu_.dzBS = glbTrack->dz(beamSpot.position());
        mu_.vx = glbTrack->vx();
        mu_.vy = glbTrack->vy();
        mu_.vz = glbTrack->vz();
      }
      reco::TrackRef trackerTrack = mu->innerTrack();
      if (trackerTrack.isNonnull()) {
        const reco::HitPattern& inhit = trackerTrack->hitPattern();
        mu_.nValidTrackerHits = inhit.numberOfValidTrackerHits();
        mu_.nValidPixelHits = inhit.numberOfValidPixelHits();
        mu_.nTrackerLayers = inhit.trackerLayersWithMeasurement();
      }
    }
    else if (mu->isStandAloneMuon()) {
      reco::TrackRef muonTrack = mu->outerTrack();
      if (muonTrack.isNonnull()) {
        const reco::HitPattern& muonhit = muonTrack->hitPattern();
        mu_.nValidMuonHits = muonhit.numberOfValidMuonHits();
      }
    }
    else if (mu->isTrackerMuon()) {
      reco::TrackRef trackerTrack = mu->innerTrack();
      if (trackerTrack.isNonnull()) {
        const reco::HitPattern& inhit = trackerTrack->hitPattern();
        mu_.normalizedChi2 = trackerTrack->normalizedChi2();
        mu_.nValidHits = trackerTrack->numberOfValidHits();
        mu_.nValidTrackerHits = inhit.numberOfValidTrackerHits();
        mu_.nValidPixelHits = inhit.numberOfValidPixelHits();
        mu_.nTrackerLayers = inhit.trackerLayersWithMeasurement();
        mu_.qoverp = trackerTrack->qoverp();
        mu_.theta = trackerTrack->theta();
        mu_.lambda = trackerTrack->lambda();
        mu_.dxy = trackerTrack->dxy();
        mu_.d0 = trackerTrack->d0();
        mu_.dsz = trackerTrack->dsz();
        mu_.dz = trackerTrack->dz();
        mu_.dxyBS = trackerTrack->dxy(beamSpot.position());
        mu_.dszBS = trackerTrack->dsz(beamSpot.position());
        mu_.dzBS = trackerTrack->dz(beamSpot.position());
        mu_.vx = trackerTrack->vx();
        mu_.vy = trackerTrack->vy();
        mu_.vz = trackerTrack->vz();
      }
    }

    if (mu->muonBestTrack().isNonnull()) {
      mu_.muonBestTrack_nTrackerLayers = mu->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
      mu_.muonBestTrack_px = mu->muonBestTrack()->px();
      mu_.muonBestTrack_py = mu->muonBestTrack()->py();
      mu_.muonBestTrack_pz = mu->muonBestTrack()->pz();
      mu_.muonBestTrack_pt = mu->muonBestTrack()->pt();
      mu_.muonBestTrack_ptError = mu->muonBestTrack()->ptError();
      if (!PVHandle->empty() && !PVHandle->front().isFake()) {
        mu_.dxyVTX = mu->muonBestTrack()->dxy(vtx.position());
        mu_.dszVTX = mu->muonBestTrack()->dsz(vtx.position());
        mu_.dzVTX = mu->muonBestTrack()->dz(vtx.position());
      }
    }

    if (mu->innerTrack().isNonnull()) {
      mu_.innerTrack_nTrackerLayers = mu->innerTrack()->hitPattern().trackerLayersWithMeasurement();
      mu_.innerTrack_px = mu->innerTrack()->px();
      mu_.innerTrack_py = mu->innerTrack()->py();
      mu_.innerTrack_pz = mu->innerTrack()->pz();
      mu_.innerTrack_pt = mu->innerTrack()->pt();
      mu_.innerTrack_ptError = mu->innerTrack()->ptError();
      if (!PVHandle->empty() && !PVHandle->front().isFake()) {
        mu_.innerTrack_dxyVTX = mu->innerTrack()->dxy(vtx.position());
        mu_.innerTrack_dszVTX = mu->innerTrack()->dsz(vtx.position());
        mu_.innerTrack_dzVTX = mu->innerTrack()->dz(vtx.position());
      }
    }

    if (mu->tunePMuonBestTrack().isNonnull()) {
      mu_.tunePMuonBestTrack_nTrackerLayers = mu->tunePMuonBestTrack()->hitPattern().trackerLayersWithMeasurement();
      mu_.tunePMuonBestTrack_px = mu->tunePMuonBestTrack()->px();
      mu_.tunePMuonBestTrack_py = mu->tunePMuonBestTrack()->py();
      mu_.tunePMuonBestTrack_pz = mu->tunePMuonBestTrack()->pz();
      mu_.tunePMuonBestTrack_pt = mu->tunePMuonBestTrack()->pt();
      mu_.tunePMuonBestTrack_ptError = mu->tunePMuonBestTrack()->ptError();
      if (!PVHandle->empty() && !PVHandle->front().isFake()) {
        mu_.tunePMuonBestTrack_dxyVTX = mu->tunePMuonBestTrack()->dxy(vtx.position());
        mu_.tunePMuonBestTrack_dszVTX = mu->tunePMuonBestTrack()->dsz(vtx.position());
        mu_.tunePMuonBestTrack_dzVTX = mu->tunePMuonBestTrack()->dz(vtx.position());
      }
    }

    mu_.isolationR03_sumpt = mu->isolationR03().sumPt;
    mu_.isolationR03_hadEt = mu->isolationR03().hadEt;
    mu_.isolationR03_emEt = mu->isolationR03().emEt;
    mu_.isolationR05_sumpt = mu->isolationR05().sumPt;
    mu_.isolationR05_hadEt = mu->isolationR05().hadEt;
    mu_.isolationR05_emEt = mu->isolationR05().emEt;
    mu_.PfChargedHadronIsoR04 = mu->pfIsolationR04().sumChargedHadronPt;
    mu_.PfNeutralHadronIsoR04 = mu->pfIsolationR04().sumNeutralHadronEt;
    mu_.PfGammaIsoR04 = mu->pfIsolationR04().sumPhotonEt;
    mu_.PfPUPtR04 = mu->pfIsolationR04().sumPUPt;
    mu_.PfChargedHadronIsoR03 = mu->pfIsolationR03().sumChargedHadronPt;
    mu_.PfNeutralHadronIsoR03 = mu->pfIsolationR03().sumNeutralHadronEt;
    mu_.PfGammaIsoR03 = mu->pfIsolationR03().sumPhotonEt;
    mu_.PfPUPtR03 = mu->pfIsolationR03().sumPUPt;

    evt_.muons.push_back(mu_);
  }
}

void Ntuplizer::fillElectrons(const edm::Event& iEvent) {
  NtupleElectron el_;

  edm::Handle<double> rhoH;
  iEvent.getByToken(rhoToken,rhoH);
  evt_.rho = *rhoH;

  edm::Handle<std::vector<reco::Conversion>> conversions;
  iEvent.getByToken(conversionsToken, conversions);

  edm::Handle<std::vector<reco::GsfElectron>> electrons;
  iEvent.getByToken(electronsToken, electrons);

  // edm::Handle<edm::ValueMap<vid::CutFlowResult>> HEEPVIDResult;
  // iEvent.getByToken(HEEPVIDResultToken, HEEPVIDResult);

  edm::Handle<edm::ValueMap<int>> nrSatCrysMap;
  edm::Handle<edm::ValueMap<float>> trkIsoMap;
  iEvent.getByToken(trkIsoMapToken, trkIsoMap);
  iEvent.getByToken(nrSatCrysMapToken, nrSatCrysMap);

  edm::Handle<edm::ValueMap<int>> nrMatchedTrkMap;
  edm::Handle<edm::ValueMap<float>> rtMatchedTrkMap;
  edm::Handle<edm::ValueMap<bool>> addGsfTrkSelMap;
  edm::Handle<edm::ValueMap<double>> EcalRecHitIsoMap;
  // edm::Handle<edm::ValueMap<double>> HcalDepth1TowerSumEtMap;
  iEvent.getByToken(nrMatchedTrkToken, nrMatchedTrkMap);
  iEvent.getByToken(rtMatchedTrkToken, rtMatchedTrkMap);
  iEvent.getByToken(addGsfTrkSelToken, addGsfTrkSelMap);
  iEvent.getByToken(EcalRecHitIsoToken, EcalRecHitIsoMap);

  evt_.nElectrons = electrons->size();
  if(electrons->size() == 0) return;

  for (std::vector<reco::GsfElectron>::const_iterator el = electrons->begin(); el != electrons->end(); ++el) {
    edm::Ptr<reco::GsfElectron> elPtr(electrons,std::distance(electrons->begin(),el));

    double R = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y() +el->superCluster()->z()*el->superCluster()->z());
    double Rt = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y());

    el_.pt = el->pt();
    el_.eta = el->eta();
    el_.rap = el->rapidity();
    el_.phi = el->phi();
    el_.en = el->energy();
    el_.charge = el->charge();

    el_.enSC = el->superCluster()->energy();
    el_.preEnSC  = el->superCluster()->preshowerEnergy();
    el_.rawEnSC = el->superCluster()->rawEnergy();
    el_.etSC =  (el->superCluster()->energy())*(Rt/R);
    el_.etaSC = el->superCluster()->eta();
    el_.phiSC = el->superCluster()->phi();
    el_.etaSCWidth = el->superCluster()->etaWidth();
    el_.phiSCWidth = el->superCluster()->phiWidth();

    el_.full5x5_sigmaIetaIeta = el->full5x5_sigmaIetaIeta();
    el_.full5x5_sigmaIphiIphi = el->full5x5_sigmaIphiIphi();
    el_.full5x5_E1x5 = el->full5x5_e1x5();
    el_.full5x5_E2x5 = el->full5x5_e2x5Max();
    el_.full5x5_E5x5 = el->full5x5_e5x5();
    el_.full5x5_hOverE = el->full5x5_hcalOverEcal();
    el_.full5x5_hOverEBC = el->full5x5_hcalOverEcalBc();
    el_.full5x5_r9 = el->full5x5_r9();

    el_.dEtaIn = el->deltaEtaSuperClusterTrackAtVtx();
    el_.dPhiIn = el->deltaPhiSuperClusterTrackAtVtx();
    // el_.dEtaSeed = el->deltaEtaSeedClusterTrackAtCalo();
    el_.dPhiSeed = el->deltaPhiSeedClusterTrackAtCalo();
    el_.dEtaEle = el->deltaEtaEleClusterTrackAtCalo();
    el_.dPhiEle = el->deltaPhiEleClusterTrackAtCalo();
    el_.dEtaSeed = el->deltaEtaSeedClusterTrackAtVtx();
    el_.EseedOverP = el->eSeedClusterOverP();
    el_.EseedOverPout = el->eSeedClusterOverPout();
    el_.EOverP = el->eSuperClusterOverP();
    el_.EeleOverPout = el->eEleClusterOverPout();
    el_.pIn = el->trackMomentumAtVtx().R();
    el_.pOut = el->trackMomentumOut().R();
    el_.ptIn = el->trackMomentumAtVtx().Rho();
    el_.ptOut = el->trackMomentumOut().Rho();

    el_.ecalEn = el->ecalEnergy();
    el_.ecalErr = el->ecalEnergyError();
    el_.trkErr = el->trackMomentumError();
    el_.combErr = el->p4Error(reco::GsfElectron::P4_COMBINATION);
    el_.PFcombErr = el->p4Error(reco::GsfElectron::P4_PFLOW_COMBINATION);

    if(el->ecalEnergy() == 0) el_.ooEmooP = 1e30;
    else if(!isfinite(el->ecalEnergy())) el_.ooEmooP = 1e30;
    else el_.ooEmooP = fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy());

    // GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
    // el_.isoChargedHadrons = pfIso.sumChargedHadronPt;
    // el_.isoNeutralHadrons = pfIso.sumNeutralHadronEt;
    // el_.isoPhotons = pfIso.sumPhotonEt;
    // el_.isoChargedFromPU = pfIso.sumPUPt;
    el_.dr03EcalRecHitSumEt = el->dr03EcalRecHitSumEt();
    el_.dr03HcalDepth1TowerSumEt = el->dr03HcalDepth1TowerSumEt();
    // const vid::CutFlowResult& HEEPCutFlowResult = (*HEEPVIDResult)[elPtr];
    // el_.passHEEP = HEEPCutFlowResult.cutFlowPassed();
    // el_.passHEEPTrkIso = HEEPCutFlowResult.getCutResultByIndex(HEEPV70::TRKISO);
    // el_.passHEEPEMHad1Iso = HEEPCutFlowResult.getCutResultByIndex(HEEPV70::EMHADD1ISO);
    // el_.passHEEPEt = HEEPCutFlowResult.getCutResultByIndex(HEEPV70::ET);
    // el_.passHEEPdEtaInSeed = HEEPCutFlowResult.getCutResultByIndex(HEEPV70::DETAINSEED);
    // el_.passHEEPdPhiIn = HEEPCutFlowResult.getCutResultByIndex(HEEPV70::DPHIIN);
    // el_.passHEEPE2x5OverE5x5 = HEEPCutFlowResult.getCutResultByIndex(HEEPV70::E2X5OVER5X5);
    // el_.passHEEPHoverE = HEEPCutFlowResult.getCutResultByIndex(HEEPV70::HADEM);
    // el_.passHEEPdxy = HEEPCutFlowResult.getCutResultByIndex(HEEPV70::DXY);
    // el_.passHEEPMissingHits = HEEPCutFlowResult.getCutResultByIndex(HEEPV70::MISSHITS);
    // el_.passHEEPEcalDriven = HEEPCutFlowResult.getCutResultByIndex(HEEPV70::ECALDRIVEN);
    // el_.passN1HEEPTrkIso = HEEPCutFlowResult.getCutFlowResultMasking(HEEPV70::TRKISO).cutFlowPassed();
    // el_.passN1HEEPEMHad1Iso = HEEPCutFlowResult.getCutFlowResultMasking(HEEPV70::EMHADD1ISO).cutFlowPassed();
    // el_.passN1HEEPEt = HEEPCutFlowResult.getCutFlowResultMasking(HEEPV70::ET).cutFlowPassed();
    // el_.passN1HEEPdEtaInSeed = HEEPCutFlowResult.getCutFlowResultMasking(HEEPV70::DETAINSEED).cutFlowPassed();
    // el_.passN1HEEPdPhiIn = HEEPCutFlowResult.getCutFlowResultMasking(HEEPV70::DPHIIN).cutFlowPassed();
    // el_.passN1HEEPE2x5OverE5x5 = HEEPCutFlowResult.getCutFlowResultMasking(HEEPV70::E2X5OVER5X5).cutFlowPassed();
    // el_.passN1HEEPHoverE = HEEPCutFlowResult.getCutFlowResultMasking(HEEPV70::HADEM).cutFlowPassed();
    // el_.passN1HEEPdxy = HEEPCutFlowResult.getCutFlowResultMasking(HEEPV70::DXY).cutFlowPassed();
    // el_.passN1HEEPMissingHits = HEEPCutFlowResult.getCutFlowResultMasking(HEEPV70::MISSHITS).cutFlowPassed();
    // el_.passN1HEEPEcalDriven = HEEPCutFlowResult.getCutFlowResultMasking(HEEPV70::ECALDRIVEN).cutFlowPassed();
    // el_.HEEPTrkIso = HEEPCutFlowResult.getValueCutUpon(HEEPV70::TRKISO);
    // el_.HEEPEMHad1Iso = HEEPCutFlowResult.getValueCutUpon(HEEPV70::EMHADD1ISO);
    el_.HEEPnrSatCrysValue = (*nrSatCrysMap)[elPtr];
    el_.HEEPTrkIsoValue = (*trkIsoMap)[elPtr];
    el_.HEEPnrMatchedTrk = (*nrMatchedTrkMap)[elPtr];
    el_.HEEPrtMatchedTrk = (*rtMatchedTrkMap)[elPtr];
    el_.HEEPaddGsfTrkSel = (*addGsfTrkSelMap)[elPtr];
    el_.HEEPEcalRecHitIsoValue = (*EcalRecHitIsoMap)[elPtr];

    reco::GsfTrackRef TrackRef = el->gsfTrack();
    // el_.lostHits = TrackRef->lost();
    el_.lostHits = TrackRef->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    el_.validHits = TrackRef->found();
    el_.nValidHits = TrackRef->hitPattern().numberOfValidHits();
    el_.nValidPixelHits = TrackRef->hitPattern().numberOfValidPixelHits();
    el_.chi2 = TrackRef->normalizedChi2();
    el_.GsfHits = TrackRef->hitPattern().trackerLayersWithMeasurement();
    el_.GsfCharge = TrackRef->charge();

    el_.d0 = TrackRef->d0();
    el_.d0Err = TrackRef->d0Error();
    el_.dxyErr = TrackRef->dxyError();
    el_.vz = TrackRef->vz();
    el_.dzErr = TrackRef->dzError();
    el_.dszErr = TrackRef->dszError();

    if (!PVHandle->empty() && !PVHandle->front().isFake()) {
      el_.dxy = TrackRef->dxy(vtx.position());
      el_.dz = TrackRef->dz(vtx.position());
      el_.dsz = TrackRef->dsz(vtx.position());
    } else {
      el_.dxy = TrackRef->dxy();
      el_.dz = TrackRef->dz();
      el_.dsz = TrackRef->dsz();
    }

    el_.Gsfpt = TrackRef->pt();
    el_.Gsfeta = TrackRef->eta();
    el_.Gsfphi = TrackRef->phi();
    el_.GsfptErr = TrackRef->ptError();
    el_.Gsfpx = TrackRef->px();
    el_.Gsfpy = TrackRef->py();
    el_.Gsfpz = TrackRef->pz();
    if (el->closestCtfTrackRef().isNonnull()) {
      el_.KFchi2 = el->closestCtfTrackRef()->normalizedChi2();
      el_.KFvalidHits = el->closestCtfTrackRef()->hitPattern().trackerLayersWithMeasurement();
      el_.KFpixHits = el->closestCtfTrackRef()->hitPattern().numberOfValidPixelHits();
      el_.KFstripHits = el->closestCtfTrackRef()->hitPattern().numberOfValidStripHits();
      el_.KFpt = el->closestCtfTrackRef()->pt();
      el_.KFeta = el->closestCtfTrackRef()->eta();
      el_.KFphi = el->closestCtfTrackRef()->phi();
      el_.KFcharge = el->closestCtfTrackRef()->charge();
      el_.KFptErr = el->closestCtfTrackRef()->ptError();
      el_.KFpx = el->closestCtfTrackRef()->px();
      el_.KFpy = el->closestCtfTrackRef()->py();
      el_.KFpz = el->closestCtfTrackRef()->pz();

      el_.KFd0 = el->closestCtfTrackRef()->d0();

      if (!PVHandle->empty() && !PVHandle->front().isFake()) {
        el_.KFdxy = el->closestCtfTrackRef()->dxy(vtx.position());
        el_.KFdz = el->closestCtfTrackRef()->dz(vtx.position());
        el_.KFdsz = el->closestCtfTrackRef()->dsz(vtx.position());
      } else {
        el_.KFdxy = el->closestCtfTrackRef()->dxy();
        el_.KFdz = el->closestCtfTrackRef()->dz();
        el_.KFdsz = el->closestCtfTrackRef()->dsz();
      }
    }
    // el_.expectedMissingInnerHits = el->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS); // 94X
    el_.expectedMissingInnerHits = el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS); // 80X

    reco::ConversionRef ConvRef = ConversionTools::matchedConversion(*el, conversions, beamSpotHandle->position());
    if (ConvRef.isNonnull()) {
      reco::Vertex ConvVtx = ConvRef.get()->conversionVertex();
      if (ConvVtx.isValid()) {
        el_.convVtxFitProb = TMath::Prob(ConvVtx.chi2(), ConvVtx.ndof());
        el_.convVtxChi2 = ConvVtx.normalizedChi2();
      }
    }
    el_.passConversionVeto = !ConversionTools::hasMatchedConversion(*el, conversions, beamSpotHandle->position());
    el_.convDist = el->convDist();
    el_.convDcot = el->convDcot();
    el_.convRadius = el->convRadius();

    el_.fbrem = el->fbrem();
    el_.fbremSC = el->superClusterFbrem();
    el_.nbrem = el->numberOfBrems();

    el_.isEB = el->isEB();
    el_.isEE = el->isEE();
    el_.isEBEEGap = el->isEBEEGap();
    el_.isEBEtaGap = el->isEBEtaGap();
    el_.isEBPhiGap = el->isEBPhiGap();
    el_.isEEDeeGap = el->isEEDeeGap();
    el_.isEERingGap = el->isEERingGap();

    el_.ecalDrivenSeed = el->ecalDrivenSeed();
    el_.trackerDrivenSeed = el->trackerDrivenSeed();

    evt_.electrons.push_back(el_);
  }
}

void Ntuplizer::fillPhotons(const edm::Event& iEvent) {
  NtuplePhoton pho_;

  edm::Handle<std::vector<reco::Photon>> photons;
  iEvent.getByToken(photonToken, photons);

  evt_.nPhotons = photons->size();
  if( evt_.nPhotons == 0 ) return;

  for (std::vector<reco::Photon>::const_iterator pho = photons->begin(); pho != photons->end(); ++pho) {
    pho_.pt = pho->pt();
    pho_.eta = pho->eta();
    pho_.phi = pho->phi();
    pho_.etaSC = pho->superCluster()->eta();
    pho_.phiSC = pho->superCluster()->phi();
    pho_.HoverE = pho->hadTowOverEm();
    pho_.hasPixelSeed = pho->hasPixelSeed();
    pho_.Full5x5_SigmaIEtaIEta = pho->full5x5_sigmaIetaIeta();
    pho_.ChIso = pho->chargedHadronIso();
    pho_.NhIso = pho->neutralHadronIso();
    pho_.PhIso = pho->photonIso();

    evt_.photons.push_back(pho_);
  }
}

void Ntuplizer::fillJets(const edm::Event& iEvent) {
  NtupleJet jet_;

  edm::Handle<std::vector<reco::PFJet>> jets;
  iEvent.getByToken(jetToken,jets);

  evt_.nJets = jets->size();
  if( evt_.nJets == 0 ) return;

  for (vector<reco::PFJet>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {
    jet_.pt        = jet->pt();
    jet_.eta       = jet->eta();
    jet_.phi       = jet->phi();
    jet_.et = jet->et();
    jet_.energy = jet->energy();
    jet_.charge = jet->charge();
    jet_.chHadEn = jet->chargedHadronEnergy();
    jet_.chHadFrac = jet->chargedHadronEnergyFraction();
    jet_.neHadEn = jet->neutralHadronEnergy();
    jet_.neHadFrac = jet->neutralHadronEnergyFraction();
    jet_.phEn = jet->photonEnergy();
    jet_.phFrac = jet->photonEnergyFraction();
    jet_.elEn = jet->electronEnergy();
    jet_.elFrac = jet->electronEnergyFraction();
    jet_.muEn = jet->muonEnergy();
    jet_.muFrac = jet->muonEnergyFraction();
    jet_.chHadMulti = jet->chargedHadronMultiplicity();
    jet_.neHadMulti = jet->neutralHadronMultiplicity();
    jet_.phMulti = jet->photonMultiplicity();
    jet_.elMulti = jet->electronMultiplicity();
    jet_.muMulti = jet->muonMultiplicity();
    jet_.chEmEn = jet->chargedEmEnergy();
    jet_.chEmFrac = jet->chargedEmEnergyFraction();
    jet_.chMuEn = jet->chargedMuEnergy();
    jet_.chMuFrac = jet->chargedMuEnergyFraction();
    jet_.neEmEn = jet->neutralEmEnergy();
    jet_.neEmFrac = jet->neutralEmEnergyFraction();
    jet_.chMulti = jet->chargedMultiplicity();
    jet_.neMulti = jet->neutralMultiplicity();

    evt_.jets.push_back(jet_);
  }
}

void Ntuplizer::fillMETs(const edm::Event& iEvent) {
  NtupleMET MET_;

  edm::Handle<std::vector<reco::PFMET>> METs;
  iEvent.getByToken(METToken,METs);

  evt_.nMETs = METs->size();
  if( evt_.nMETs == 0 ) return;

  for (vector<reco::PFMET>::const_iterator MET = METs->begin(); MET != METs->end(); ++MET) {
    MET_.Type1_pt    = METs->front().pt();
    MET_.Type1_phi   = METs->front().phi();
    MET_.Type1_px    = METs->front().px();
    MET_.Type1_py    = METs->front().py();
    MET_.Type1_sumEt = METs->front().sumEt();

    evt_.METs.push_back(MET_);
  }
}

void Ntuplizer::fillGenParticles(const edm::Event& iEvent) {
  NtupleGenParticle particle_;

  edm::Handle<std::vector<reco::GenParticle>> particles;
  iEvent.getByToken(genparticlesToken, particles);

  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(generatorToken, genInfo);

  evt_.weight = genInfo->weight();

  evt_.nGenParticles = particles->size();
  if( particles->size() == 0 ) return;

  std::vector<int> ExcludePdgID {211,111,213,321,311,113,223,221,2112,310,130,2212};
  for (std::vector<reco::GenParticle>::const_iterator particle = particles->begin(); particle != particles->end(); ++particle) {
    auto pid = std::abs(particle->pdgId());
    if ( std::find(ExcludePdgID.begin(),ExcludePdgID.end(),pid)!=ExcludePdgID.end() ) continue;
    if ( pid==22 || pid==21 || pid==1 || pid==2 || pid==3 || pid==4 || pid==5 ) {
      if ( !( particle->isHardProcess() || particle->fromHardProcessBeforeFSR() ) ) continue;
    }

    particle_.id = particle->pdgId();
    particle_.pt = particle->pt();
    particle_.px = particle->px();
    particle_.py = particle->py();
    particle_.pz = particle->pz();
    particle_.eta = particle->eta();
    particle_.phi = particle->phi();
    particle_.energy = particle->energy();
    particle_.mass = particle->mass();
    particle_.charge = particle->charge();
    particle_.status = particle->status();
    particle_.mother = particle->mother()->pdgId();
    particle_.motherPt = particle->mother()->pt();

    particle_.isDirectHardProcessTauDecayProductFinalState = particle->isDirectHardProcessTauDecayProductFinalState();
    particle_.isDirectPromptTauDecayProductFinalState = particle->isDirectPromptTauDecayProductFinalState();
    particle_.isHardProcess = particle->isHardProcess();
    particle_.isLastCopy = particle->isLastCopy();
    particle_.isLastCopyBeforeFSR = particle->isLastCopyBeforeFSR();
    particle_.isPromptDecayed = particle->isPromptDecayed();
    particle_.isPromptFinalState = particle->isPromptFinalState();
    particle_.fromHardProcessFinalState = particle->fromHardProcessFinalState();
    particle_.fromHardProcessDecayed = particle->fromHardProcessDecayed();
    particle_.fromHardProcessBeforeFSR = particle->fromHardProcessBeforeFSR();

    evt_.genparticles.push_back(particle_);
  }
}

void Ntuplizer::fillGenJets(const edm::Event& iEvent) {
  NtupleGenJet genjet_;

  edm::Handle<std::vector<reco::GenJet>> genjets;
  iEvent.getByToken(genjetToken, genjets);

  evt_.nGenJets = genjets->size();
  if(evt_.nGenJets == 0) return;

  for (std::vector<reco::GenJet>::const_iterator genjet = genjets->begin(); genjet != genjets->end(); ++genjet) {
    genjet_.pt = genjet->pt();
    genjet_.eta = genjet->eta();
    genjet_.phi = genjet->phi();
    genjet_.et = genjet->et();
    genjet_.energy = genjet->energy();
    genjet_.charge = genjet->charge();
    genjet_.nConst = genjet->nConstituents();
    genjet_.emEnergy = genjet->emEnergy();
    genjet_.hadEnergy = genjet->hadEnergy();
    genjet_.invisibleEnergy = genjet->invisibleEnergy();
    genjet_.auxiliaryEnergy = genjet->auxiliaryEnergy();

    evt_.genjets.push_back(genjet_);
  }
}

void Ntuplizer::fillGsfTracks(const edm::Event& iEvent) {
  NtupleGsfTrack trk_;

  edm::Handle<std::vector<reco::GsfTrack>> trks;
  iEvent.getByToken(GsfTrackToken, trks);

  evt_.nGsfTrks = trks->size();
  if( evt_.nGsfTrks == 0 ) return;

  for (std::vector<reco::GsfTrack>::const_iterator trk = trks->begin(); trk != trks->end(); ++trk) {
    trk_.lostHits = trk->lost();
    trk_.validHits = trk->found();
    trk_.nValidHits = trk->hitPattern().numberOfValidHits();
    trk_.nValidPixelHits = trk->hitPattern().numberOfValidPixelHits();
    trk_.chi2 = trk->normalizedChi2();
    trk_.GsfHits = trk->hitPattern().trackerLayersWithMeasurement();
    trk_.GsfCharge = trk->charge();
    trk_.d0 = trk->d0();
    trk_.d0Err = trk->d0Error();
    trk_.dxy = trk->dxy();
    trk_.dxyErr = trk->dxyError();
    trk_.dz = trk->dz();
    trk_.vz = trk->vz();
    trk_.dzErr = trk->dzError();
    trk_.dsz = trk->dsz();
    trk_.dszErr = trk->dszError();
    trk_.Gsfpt = trk->pt();
    trk_.Gsfeta = trk->eta();
    trk_.Gsfphi = trk->phi();
    trk_.GsfptErr = trk->ptError();
    trk_.Gsfpx = trk->px();
    trk_.Gsfpy = trk->py();
    trk_.Gsfpz = trk->pz();

    evt_.GsfTracks.push_back(trk_);
  }
}

void Ntuplizer::fillGeneralTracks(const edm::Event& iEvent) {
  NtupleGeneralTrack trk_;

  edm::Handle<std::vector<reco::Track>> trks;
  iEvent.getByToken(GeneralTrackToken, trks);

  evt_.nTrks = trks->size();
  if(evt_.nTrks == 0) return;

  for (std::vector<reco::Track>::const_iterator trk = trks->begin(); trk != trks->end(); ++trk) {
    if(trk->pt() < 5.) continue;
    trk_.pt = trk->pt();
    trk_.eta = trk->eta();
    trk_.phi = trk->phi();
    trk_.chi2 = trk->chi2();
    trk_.ndof = trk->ndof();
    trk_.normChi2 = trk->normalizedChi2();
    trk_.charge = trk->charge();
    trk_.qoverp = trk->qoverp();
    trk_.qoverpErr = trk->qoverpError();
    trk_.dxy = trk->dxy();
    trk_.d0 = trk->d0();
    trk_.dz = trk->dz();
    trk_.dsz = trk->dsz();
    trk_.vz = trk->vz();
    trk_.ptErr = trk->ptError();
    trk_.nValidHits = trk->hitPattern().numberOfValidHits();
    trk_.nValidPixelHits = trk->hitPattern().numberOfValidPixelHits();
    trk_.isHighPurity = trk->quality(TrackBase::highPurity);
    trk_.isJetCoreRegionalAlgo = trk->isAlgoInMask(TrackBase::jetCoreRegionalStep);
    trk_.qualMask = static_cast<int>(trk->qualityMask());
    trk_.algoMask = static_cast<int>(trk->algoMaskUL());

    evt_.tracks.push_back(trk_);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
