#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TEfficiency.h>
#include <iostream>
#include "Analysis/Ntuplizer/interface/Ntuplizer.h"
#include "Physics.h"

using namespace std;

void analysis(TString DataSet)
{
  TChain* chain = new TChain("tree/PhysicsTree");
  if(DataSet=="QCD_Pt_50to80") chain->Add("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sako/QCD/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/QCD_Pt-50to80_Ntuplizer/200318_142748/0000/Ntuple_*.root");
  if(DataSet=="QCD_EMEnriched_Pt_30to50") chain->Add("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sako/QCD/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/QCD_EMEnriched_Pt-30to50_Ntuplizer/200409_043734/0000/Ntuple_*.root");
  if(DataSet=="QCD_EMEnriched_Pt_50to80") chain->Add("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sako/QCD/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/QCD_EMEnriched_Pt-50to80_Ntuplizer_resub_v3/200424_094932/0000/Ntuple_*.root");
  if(DataSet=="QCD_EMEnriched_Pt_80to120") chain->Add("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sako/QCD/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/QCD_EMEnriched_Pt-80to120_Ntuplizer/200406_060357/0000/Ntuple_*.root");
  if(DataSet=="QCD_EMEnriched_Pt_120to170") chain->Add("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sako/QCD/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/QCD_EMEnriched_Pt-120to170_Ntuplizer_resub_v3/200424_095123/0000/Ntuple_*.root");
  if(DataSet=="QCD_EMEnriched_Pt_170to300") chain->Add("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sako/QCD/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/QCD_EMEnriched_Pt-170to300_Ntuplizer/200406_060252/0000/Ntuple_*.root");
  if(DataSet=="QCD_EMEnriched_Pt_300toInf") chain->Add("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sako/QCD/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/QCD_EMEnriched_Pt-300toInf_Ntuplizer/200402_040649/0000/Ntuple_*.root");

  if(DataSet=="ZZTo4L") chain->Add("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sako/ZZTo4L/ZZTo4L_13TeV_powheg_pythia8/ZZTo4L_Ntuplizer_resub_v4/201027_110459/0000/Ntuple_*.root");
  if(DataSet=="ZGTo2NuG_PtG130") chain->Add("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sako/ZG/ZGTo2NuG_PtG-130_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZGTo2NuG_PtG-130_Ntuplizer_resub_v4/201027_110730/0000/Ntuple_*.root");
  if(DataSet=="ZGTo2NuG") chain->Add("dcap://cluster142.knu.ac.kr//pnfs/knu.ac.kr/data/cms/store/user/sako/ZG/ZGTo2NuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/ZGTo2NuG_Ntuplizer_resub_v4/201027_110905/0000/Ntuple_*.root");

  if(DataSet=="H200A1") chain->Add("/u/user/sako/ModHEEP/CMSSW_8_0_32/src/Analysis/Ntuplizer/configs/H200A1.root");
  if(DataSet=="H800A1") chain->Add("/u/user/sako/ModHEEP/CMSSW_8_0_32/src/Analysis/Ntuplizer/configs/H800A1.root");
  if(DataSet=="H2000A1") chain->Add("/u/user/sako/ModHEEP/CMSSW_8_0_32/src/Analysis/Ntuplizer/configs/H2000A1.root");

  if(DataSet=="cmsShow") chain->Add("/u/user/sako/ModHEEP/CMSSW_8_0_32/src/Analysis/Ntuplizer/configs/ZGPt130_cmsShow.root");

  if(DataSet=="test") chain->Add("Ntuple_1.root");

  NtupleEvent* evt = new NtupleEvent();
  chain->SetBranchAddress("Event",&evt);

  double ptHigh = 1000.;
  double invMass = 10.;
  if (DataSet=="ZZTo4L") invMass = 200;

  TH1D* tNevt = new TH1D("nEvt","Num of Evt;;",1,0.,2.);
  TH1D* tNevt_1 = new TH1D("nEvt_1","Num of Eles;;",1,0.,2.);
  TH1D* tNevt_2 = new TH1D("nEvt_2","Num of Eles with add GSF track;;",1,0.,2.);
  TH1D* tNevt_3 = new TH1D("nEvt_3","Num of Eles with HEEP modified Iso;;",1,0.,2.);

  TH1D* tNevt_4 = new TH1D("nEvt_4","Num of Eles;;",1,0.,2.);
  TH1D* tNevt_5 = new TH1D("nEvt_5","Num of Eles;;",1,0.,2.);
  TH1D* tNevt_6 = new TH1D("nEvt_6","Num of Eles;;",1,0.,2.);
  TH1D* tNevt_7 = new TH1D("nEvt_7","Num of Eles;;",1,0.,2.);
  TH1D* tNevt_8 = new TH1D("nEvt_8","Num of Eles;;",1,0.,2.);
  // TH1D* tGENPartonPt = new TH1D("GENPartonPt","GEN-lv parton Pt;p_{T};",200,0.0,ptHigh);
  // tGENPartonPt->Sumw2();
  // TH1D* tGENPartonEta = new TH1D("GENPartonEta","GEN-lv parton #eta;#eta;",120,-3.,3.);
  // tGENPartonEta->Sumw2();
  // TH1D* tGENPartonPhi = new TH1D("GENPartonPhi","GEN-lv parton #phi;#phi;",100,-4.,4.);
  // tGENPartonPhi->Sumw2();

  TH1D* tElePt = new TH1D("ElePt","RECO-lv electron Pt;p_{T};",200,0.0,ptHigh);
  tElePt->Sumw2();
  TH1D* tEleEta = new TH1D("EleEta","RECO-lv electron #eta;#eta;",120,-3.,3.);
  tEleEta->Sumw2();
  // TH1D* tHEEPnoIsoPt = new TH1D("HEEPnoIsoPt","HEEP electron Pt w/o Iso;p_{T};",200,0.,ptHigh);
  // tHEEPnoIsoPt->Sumw2();
  // TH1D* tHEEPnoIsoEta = new TH1D("HEEPnoIsoEta","HEEP electron #eta w/o Iso;#eta;",120,-3.,3.);
  // tHEEPnoIsoEta->Sumw2();
  // TH1D* tHEEPPt = new TH1D("HEEPPt","HEEP electron Pt;p_{T};",200,0.,ptHigh);
  // tHEEPPt->Sumw2();
  // TH1D* tHEEPEta = new TH1D("HEEPEta","HEEP electron #eta;#eta;",120,-3.,3.);
  // tHEEPEta->Sumw2();
  // TH1D* tHEEPmodIsoPt = new TH1D("HEEPmodIsoPt","HEEP electron Pt w/ modified Iso;p_{T};",200,0.,ptHigh);
  // tHEEPmodIsoPt->Sumw2();
  // TH1D* tHEEPmodIsoEta = new TH1D("HEEPmodIsoEta","HEEP electron #eta w/ modified Iso;#eta;",120,-3.,3.);
  // tHEEPmodIsoEta->Sumw2();

  TH1D* t_etaSCWidth = new TH1D("etaSCWidth","#eta SC Width",500,0.,0.1); t_etaSCWidth->Sumw2();
  TH1D* t_phiSCWidth = new TH1D("phiSCWidth","#phi SC Width",500,0.,0.1); t_phiSCWidth->Sumw2();
  TH1D* t_full5x5_sigmaIetaIeta = new TH1D("full5x5_sigmaIetaIeta","full5x5 #sigma_{I#etaI#eta}",400,0.,0.04); t_full5x5_sigmaIetaIeta->Sumw2();
  TH1D* t_full5x5_sigmaIphiIphi = new TH1D("full5x5_sigmaIphiIphi","full5x5 #sigma_{I#phi#phi}",1000,0.,0.5); t_full5x5_sigmaIphiIphi->Sumw2();
  TH1D* t_full5x5_E1x5 = new TH1D("full5x5_E1x5","full5x5_E1x5",200,0.,1.); t_full5x5_E1x5->Sumw2();
  TH1D* t_full5x5_E2x5 = new TH1D("full5x5_E2x5","full5x5_E2x5",200,0.,1.); t_full5x5_E2x5->Sumw2();
  TH1D* t_full5x5_E5x5 = new TH1D("full5x5_E5x5","full5x5_E5x5",1000,0.,1000.); t_full5x5_E5x5->Sumw2();
  TH1D* t_full5x5_r9 = new TH1D("full5x5_r9","full5x5_r9",200,0.,1.); t_full5x5_r9->Sumw2();
  TH1D* t_dPhiIn = new TH1D("dPhiIn","d#phi_{In}",200,-0.1,0.1); t_dPhiIn->Sumw2();
  TH1D* t_dEtaSeed = new TH1D("dEtaSeed","d#eta_{Seed}",200,-0.02,0.02); t_dEtaSeed->Sumw2();

  TH1D* t_dxy = new TH1D("dxy","dxy",500,-0.05,0.05); t_dxy->Sumw2();
  TH1D* t_dz = new TH1D("dz","dz",500,-0.1,0.1); t_dz->Sumw2();

  TH1D* t_KFdxy = new TH1D("KFdxy","KFdxy",500,-0.05,0.05); t_KFdxy->Sumw2();
  TH1D* t_KFdz = new TH1D("KFdz","KFdz",500,-0.1,0.1); t_KFdz->Sumw2();

  TH1D* t_HEEPnoSelectedGsfTrk = new TH1D("HEEPnoSelectedGsfTrk","# of additional Gsf trks candidates",10,0,10); t_HEEPnoSelectedGsfTrk->Sumw2();
  TH1D* t_HEEPaddGsfTrk_Gsfpt = new TH1D("HEEPaddGsfTrk_Gsfpt","p_{T} of 2nd Gsf trk",200,0.0,ptHigh); t_HEEPaddGsfTrk_Gsfpt->Sumw2();
  TH1D* t_HEEPaddGsfTrk_Gsfeta = new TH1D("HEEPaddGsfTrk_Gsfeta","#eta of 2nd Gsf trk",120,-3.,3.); t_HEEPaddGsfTrk_Gsfeta->Sumw2();
  TH1D* t_HEEPaddGsfTrk_Gsfphi = new TH1D("HEEPaddGsfTrk_Gsfphi","#phi of 2nd Gsf trk",100,-4.,4.); t_HEEPaddGsfTrk_Gsfphi->Sumw2();
  TH1D* t_HEEPaddGsfTrk_GsfptErr = new TH1D("HEEPaddGsfTrk_GsfptErr","Err(p_{T})/p_{T} of 2nd Gsf trk",200,0.,1.); t_HEEPaddGsfTrk_GsfptErr->Sumw2();
  TH1D* t_HEEPaddGsfTrk_lostHits = new TH1D("HEEPaddGsfTrk_lostHits","# of list hits of 2nd Gsf trk",20,0.,20.); t_HEEPaddGsfTrk_lostHits->Sumw2();
  TH1D* t_HEEPaddGsfTrk_nValidHits = new TH1D("HEEPaddGsfTrk_nValidHits","# of valid hits of 2nd Gsf trk",50,0.,50.); t_HEEPaddGsfTrk_nValidHits->Sumw2();
  TH1D* t_HEEPaddGsfTrk_nValidPixelHits = new TH1D("HEEPaddGsfTrk_nValidPixelHits","# of valid pixel hits of 2nd Gsf trk",10,0.,10.); t_HEEPaddGsfTrk_nValidPixelHits->Sumw2();
  TH1D* t_HEEPaddGsfTrk_chi2 = new TH1D("HEEPaddGsfTrk_chi2","Normalized #chi^2 of 2nd Gsf trk",200,0.,20.); t_HEEPaddGsfTrk_chi2->Sumw2();
  TH1D* t_HEEPaddGsfTrk_d0 = new TH1D("HEEPaddGsfTrk_d0","d0 of 2nd Gsf trk",500,-0.5,0.5); t_HEEPaddGsfTrk_d0->Sumw2();
  TH1D* t_HEEPaddGsfTrk_dxy = new TH1D("HEEPaddGsfTrk_dxy","dxy of 2nd Gsf trk",500,-0.1,0.1); t_HEEPaddGsfTrk_dxy->Sumw2();
  TH1D* t_HEEPaddGsfTrk_dz = new TH1D("HEEPaddGsfTrk_dz","dz of 2nd Gsf trk",500,-0.1,0.1); t_HEEPaddGsfTrk_dz->Sumw2();

  TH1D* t_convVtxFitProb = new TH1D("convVtxFitProb","convVtxFitProb",500,0.,1.); t_convVtxFitProb->Sumw2();
  TH1D* t_convVtxChi2 = new TH1D("convVtxChi2","convVtxChi2",200,0.,20.); t_convVtxChi2->Sumw2();
  TH1D* t_passConversionVeto = new TH1D("passConversionVeto","passConversionVeto",2,0.,2.); t_passConversionVeto->Sumw2();
  TH1D* t_convDist = new TH1D("convDist","convDist",500,-0.5,0.5); t_convDist->Sumw2();
  TH1D* t_convDcot = new TH1D("convDcot","convDcot",500,-0.5,0.5); t_convDcot->Sumw2();
  TH1D* t_convRadius = new TH1D("convRadius","convRadius",500,-100.,100.); t_convRadius->Sumw2();

  TH1D* t_HEEPaddVtx_isValid = new TH1D("HEEPaddVtx_isValid","HEEPaddVtx_isValid",2,0.,2.); t_HEEPaddVtx_isValid->Sumw2();
  TH1D* t_HEEPaddVtx_dx = new TH1D("HEEPaddVtx_dx","HEEPaddVtx_dx",500,-0.5,0.5); t_HEEPaddVtx_dx->Sumw2();
  TH1D* t_HEEPaddVtx_dy = new TH1D("HEEPaddVtx_dy","HEEPaddVtx_dy",500,-0.5,0.5); t_HEEPaddVtx_dy->Sumw2();
  TH1D* t_HEEPaddVtx_dz = new TH1D("HEEPaddVtx_dz","HEEPaddVtx_dz",500,-0.5,0.5); t_HEEPaddVtx_dz->Sumw2();
  TH1D* t_HEEPaddVtx_chi2 = new TH1D("HEEPaddVtx_chi2","HEEPaddVtx_chi2",200,0.,20.); t_HEEPaddVtx_chi2->Sumw2();
  TH1D* t_HEEPaddVtx_pt = new TH1D("HEEPaddVtx_pt","HEEPaddVtx_pt",200,0.0,ptHigh); t_HEEPaddVtx_pt->Sumw2();
  TH1D* t_HEEPaddVtx_rapidity = new TH1D("HEEPaddVtx_rapidity","HEEPaddVtx_rapidity",120,-3.,3.); t_HEEPaddVtx_rapidity->Sumw2();
  TH1D* t_HEEPaddVtx_phi = new TH1D("HEEPaddVtx_phi","HEEPaddVtx_phi",100,-4.,4.); t_HEEPaddVtx_phi->Sumw2();
  TH1D* t_HEEPaddVtx_M = new TH1D("HEEPaddVtx_M","HEEPaddVtx_M",500,0.,invMass); t_HEEPaddVtx_M->Sumw2();

  //
  // TEfficiency* tHEEPmodIsoPtEff = new TEfficiency("HEEPmodIsoPtEff","HEEP modified Iso eff vs Pt;p_{T};#epsilon",100,0.,ptHigh);
  // TEfficiency* tHEEPmodIsoEtaEff = new TEfficiency("HEEPmodIsoEtaEff","HEEP modified Iso eff vs #eta;#eta;#epsilon",60,-3.,3.);
  // TEfficiency* tHEEPIsoPtEff = new TEfficiency("HEEPIsoPtEff","HEEP Iso eff vs Pt;p_{T};#epsilon",100,0.,ptHigh);
  // TEfficiency* tHEEPIsoEtaEff = new TEfficiency("HEEPIsoEtaEff","HEEP Iso eff vs #eta;#eta;#epsilon",60,-3.,3.);

  for (unsigned int i = 0; i < chain->GetEntries(); i++) {
    chain->GetEntry(i);

    if (i % 100000 == 0) printf("Analyzing %dth event ...\n", i);

    std::vector<PhysicsGenParticle> genptcs;
    std::vector<PhysicsElectron> eles;

    tNevt->Fill(1);

    // for (auto ngenptc : evt->genparticles) {
    //   PhysicsGenParticle* genptc = (PhysicsGenParticle*)&ngenptc;
    //   int pid = std::abs(genptc->id);
    //
    //   if ( genptc->isHardProcess // fromHardProcessBeforeFSR
    //     && ( pid==1 || pid==2 || pid==3 || pid==4 || pid==5 || pid==21 ) ) {
    //
    //     genptcs.push_back(*genptc);
    //   }
    // }

    for (auto nele : evt->electrons) {
      PhysicsElectron* ele = (PhysicsElectron*)&nele;
      if ( ele->Accep(50., 1.4442) && ele->TransitionVeto() ) eles.push_back(*ele); // originally accep |eta| < 2.5
    }

    // for (auto genptc : genptcs) {
    //   tGENPartonPt->Fill(genptc.pt);
    //   tGENPartonEta->Fill(genptc.eta);
    //   tGENPartonPhi->Fill(genptc.phi);
    // }

    for (auto ele : eles) {

      tNevt_1->Fill(1);

      // if (ele.HEEPnoIso()) {
      //   tHEEPnoIsoPt->Fill(ele.pt);
      //   tHEEPnoIsoEta->Fill(ele.eta);
      //
      //   if (ele.HEEPaddGsfTrkSel) {
      //     if ( (ele.HEEPTrkIsoValue < 5.) && ele.EMHad1Iso(evt->rho) ) {
      //       tHEEPmodIsoPt->Fill(ele.pt);
      //       tHEEPmodIsoEta->Fill(ele.eta);
      //     }
      //     tHEEPmodIsoPtEff->Fill( (ele.HEEPTrkIsoValue < 5.) && ele.EMHad1Iso(evt->rho) , ele.pt );
      //     tHEEPmodIsoEtaEff->Fill( (ele.HEEPTrkIsoValue < 5.) && ele.EMHad1Iso(evt->rho) , ele.eta );
      //   } else {
      //     if ( (ele.HEEPTrkIsoValue < 5.) && ele.EMHad1IsoCustom(ele.dr03EcalRecHitSumEt,ele.dr03HcalDepth1TowerSumEt,evt->rho) ) {
      //       tHEEPPt->Fill(ele.pt);
      //       tHEEPEta->Fill(ele.eta);
      //     }
      //     tHEEPIsoPtEff->Fill( (ele.HEEPTrkIsoValue < 5.) && ele.EMHad1IsoCustom(ele.dr03EcalRecHitSumEt,ele.dr03HcalDepth1TowerSumEt,evt->rho) , ele.pt );
      //     tHEEPIsoEtaEff->Fill( (ele.HEEPTrkIsoValue < 5.) && ele.EMHad1IsoCustom(ele.dr03EcalRecHitSumEt,ele.dr03HcalDepth1TowerSumEt,evt->rho) , ele.eta );
      //   }
      // }

      if (ele.HEEPaddGsfTrkSel) {
        tNevt_2->Fill(1);
        if ((ele.HEEPTrkIsoValue < 5.) && ele.EMHad1Iso(evt->rho)) {
          tNevt_3->Fill(1);

          if (ele.ecalDrivenSeed) {
            tNevt_4->Fill(1);

            if (std::abs(ele.dPhiIn) < 0.06) {
              tNevt_5->Fill(1);

              if (ele.lostHits <= 1) {
                tNevt_6->Fill(1);

                if ( ele.full5x5_hOverE < (1./ele.en + 0.05) ) {
                  tNevt_7->Fill(1);

                  if (std::abs(ele.dxy) < 0.02) {
                    tNevt_8->Fill(1);
                  }
                }
              }
            }
          }

          if (ele.HEEPnoIn()) {
            tElePt->Fill(ele.pt);
            tEleEta->Fill(ele.eta);
            t_etaSCWidth->Fill(ele.etaSCWidth);
            t_phiSCWidth->Fill(ele.phiSCWidth);
            t_full5x5_sigmaIetaIeta->Fill(ele.full5x5_sigmaIetaIeta);
            t_full5x5_sigmaIphiIphi->Fill(ele.full5x5_sigmaIphiIphi);
            t_full5x5_E1x5->Fill(ele.full5x5_E1x5/ele.full5x5_E5x5);
            t_full5x5_E2x5->Fill(ele.full5x5_E2x5/ele.full5x5_E5x5);
            t_full5x5_E5x5->Fill(ele.full5x5_E5x5);
            t_full5x5_r9->Fill(ele.full5x5_r9);
            t_dPhiIn->Fill(ele.dPhiIn);
            t_dEtaSeed->Fill(ele.dEtaSeed);

            t_dxy->Fill(ele.dxy);
            t_dz->Fill(ele.dz);


            t_KFdxy->Fill(ele.KFdxy);
            t_KFdz->Fill(ele.KFdz);

            t_HEEPnoSelectedGsfTrk->Fill(ele.HEEPnoSelectedGsfTrk);
            t_HEEPaddGsfTrk_Gsfpt->Fill(ele.HEEPaddGsfTrk_Gsfpt);
            t_HEEPaddGsfTrk_Gsfeta->Fill(ele.HEEPaddGsfTrk_Gsfeta);
            t_HEEPaddGsfTrk_Gsfphi->Fill(ele.HEEPaddGsfTrk_Gsfphi);
            t_HEEPaddGsfTrk_GsfptErr->Fill(ele.HEEPaddGsfTrk_GsfptErr/ele.HEEPaddGsfTrk_Gsfpt);
            t_HEEPaddGsfTrk_lostHits->Fill(ele.HEEPaddGsfTrk_lostHits);
            t_HEEPaddGsfTrk_nValidHits->Fill(ele.HEEPaddGsfTrk_nValidHits);
            t_HEEPaddGsfTrk_nValidPixelHits->Fill(ele.HEEPaddGsfTrk_nValidPixelHits);
            t_HEEPaddGsfTrk_chi2->Fill(ele.HEEPaddGsfTrk_chi2);
            t_HEEPaddGsfTrk_d0->Fill(ele.HEEPaddGsfTrk_d0);
            t_HEEPaddGsfTrk_dxy->Fill(ele.HEEPaddGsfTrk_dxy);
            t_HEEPaddGsfTrk_dz->Fill(ele.HEEPaddGsfTrk_dz);

            t_convVtxFitProb->Fill(ele.convVtxFitProb);
            t_convVtxChi2->Fill(ele.convVtxChi2);
            t_passConversionVeto->Fill(ele.passConversionVeto);
            t_convDist->Fill(ele.convDist);
            t_convDcot->Fill(ele.convDcot);
            t_convRadius->Fill(ele.convRadius);

            t_HEEPaddVtx_isValid->Fill(ele.HEEPaddVtx_isValid);

            if (ele.HEEPaddVtx_isValid) {
              t_HEEPaddVtx_dx->Fill(ele.HEEPaddVtx_dx);
              t_HEEPaddVtx_dy->Fill(ele.HEEPaddVtx_dy);
              t_HEEPaddVtx_dz->Fill(ele.HEEPaddVtx_dz);
              t_HEEPaddVtx_chi2->Fill(ele.HEEPaddVtx_chi2);
              t_HEEPaddVtx_pt->Fill(ele.HEEPaddVtx_pt);
              t_HEEPaddVtx_rapidity->Fill(ele.HEEPaddVtx_rapidity);
              t_HEEPaddVtx_phi->Fill(ele.HEEPaddVtx_phi);
              t_HEEPaddVtx_M->Fill(ele.HEEPaddVtx_M);
            }

            if (DataSet=="cmsShow") std::cout << "run = " << evt->run << " lumi = " << evt->lumi << " event = " << evt->event << std::endl;
          }
        }
      }
    }
  }

  TFile* f = new TFile(DataSet+"_result.root","RECREATE");

  tNevt->Write();
  tNevt_1->Write();
  tNevt_2->Write();
  tNevt_3->Write();

  tNevt_4->Write();
  tNevt_5->Write();
  tNevt_6->Write();
  tNevt_7->Write();
  tNevt_8->Write();
  // tGENPartonPt->Write();
  // tGENPartonEta->Write();
  // tGENPartonPhi->Write();
  tElePt->Write();
  tEleEta->Write();

  t_etaSCWidth->Write();
  t_phiSCWidth->Write();
  t_full5x5_sigmaIetaIeta->Write();
  t_full5x5_sigmaIphiIphi->Write();
  t_full5x5_E1x5->Write();
  t_full5x5_E2x5->Write();
  t_full5x5_E5x5->Write();
  t_full5x5_r9->Write();
  t_dPhiIn->Write();
  t_dEtaSeed->Write();

  t_dxy->Write();
  t_dz->Write();
  t_KFdxy->Write();
  t_KFdz->Write();

  t_HEEPnoSelectedGsfTrk->Write();
  t_HEEPaddGsfTrk_Gsfpt->Write();
  t_HEEPaddGsfTrk_Gsfeta->Write();
  t_HEEPaddGsfTrk_Gsfphi->Write();
  t_HEEPaddGsfTrk_GsfptErr->Write();
  t_HEEPaddGsfTrk_lostHits->Write();
  t_HEEPaddGsfTrk_nValidHits->Write();
  t_HEEPaddGsfTrk_nValidPixelHits->Write();
  t_HEEPaddGsfTrk_chi2->Write();
  t_HEEPaddGsfTrk_d0->Write();
  t_HEEPaddGsfTrk_dxy->Write();
  t_HEEPaddGsfTrk_dz->Write();

  t_convVtxFitProb->Write();
  t_convVtxChi2->Write();
  t_passConversionVeto->Write();
  t_convDist->Write();
  t_convDcot->Write();
  t_convRadius->Write();

  t_HEEPaddVtx_isValid->Write();
  t_HEEPaddVtx_dx->Write();
  t_HEEPaddVtx_dy->Write();
  t_HEEPaddVtx_dz->Write();
  t_HEEPaddVtx_chi2->Write();
  t_HEEPaddVtx_pt->Write();
  t_HEEPaddVtx_rapidity->Write();
  t_HEEPaddVtx_phi->Write();
  t_HEEPaddVtx_M->Write();

  // tHEEPnoIsoPt->Write();
  // tHEEPnoIsoEta->Write();
  // tHEEPmodIsoPt->Write();
  // tHEEPmodIsoEta->Write();
  // tHEEPPt->Write();
  // tHEEPEta->Write();
  // tHEEPmodIsoPtEff->Write();
  // tHEEPmodIsoEtaEff->Write();
  // tHEEPIsoPtEff->Write();
  // tHEEPIsoEtaEff->Write();

  f->Close();

  std::cout << "Finished" << std::endl;

  return;
}
