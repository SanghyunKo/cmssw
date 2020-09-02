#ifndef Physics_H
#define Physics_H

#include <TLorentzVector.h>
#include "Analysis/Ntuplizer/interface/Ntuplizer.h"

const double muon_mass     = 0.1056583715;
const double electron_mass = 0.000510998;

class PhysicsElectron : public NtupleElectron
{
  public :
  bool Accep(Double_t Pt_, Double_t Eta_) {
    if( pt > Pt_ && fabs(etaSC) < Eta_) return true;
    else return false;
  }

  bool TransitionVeto() {
    if( fabs(etaSC) > 1.4442 && fabs(etaSC) < 1.566 ) return false;
    else return true;
  }

  double dEtOverEt(double EtGen1, double EtGen2) {
    double Et = fabs(etSC - EtGen1 - EtGen2)/(etSC);
    return Et;
  }

  bool TrackCut(double vz_, int nValidHits_, int nValidPixelHits_) {
    if(eta > 1.5) return fabs(vz-vz_)<0.5 && nValidHits_>7 && nValidPixelHits_>0;
    else return fabs(vz-vz_)<0.1 && nValidHits_>7 && nValidPixelHits_>0;
  }

  bool EMHad1Iso(double rho_) {
    double val = HEEPEcalRecHitIsoValue+dr03HcalDepth1TowerSumEt;
    if(eta < 1.5) return val < 2+0.03*etSC+0.28*rho_;
    else {
      if(etSC < 50) return val < 2.5+0.28*rho_;
      else return val < 2.5+0.03*(etSC-50.)+0.28*rho_;
    }
  }

  bool EMHad1IsoCustom(double ecal, double hcal, double rho_) {
    double val = ecal+hcal;
    if(eta < 1.5) return val < 2+0.03*etSC+0.28*rho_;
    else {
      if(etSC < 50) return val < 2.5+0.28*rho_;
      else return val < 2.5+0.03*(etSC-50.)+0.28*rho_;
    }
  }

  bool HEEPnoIso() {
    if ( etSC < 35. ) return false;
    if ( !ecalDrivenSeed ) return false;
    if ( !(std::abs(dPhiIn) < 0.06) ) return false;
    if ( !(lostHits <= 1) ) return false;
    if ( std::abs(etaSC) < 1.4442 ) {
      bool HoEcut = ( full5x5_hOverE < (1./en + 0.05) );
      bool E2x5cut = ( full5x5_E2x5/full5x5_E5x5 > 0.94 || full5x5_E1x5/full5x5_E5x5 > 0.83 );

      return ( (std::abs(dEtaSeed) < 0.004) && HoEcut && (std::abs(dxy) < 0.02) ); // E2x5 excluded
    } else if ( std::abs(etaSC) > 1.566 && std::abs(etaSC) < 2.5 ) {
      bool HoEcut = ( full5x5_hOverE < (5./en + 0.05) );
      bool E5x5cut = ( full5x5_E5x5 < 0.03 );

      return ( (std::abs(dEtaSeed) < 0.006) && HoEcut && E5x5cut && (std::abs(dxy) < 0.05) );
    }

    return false;
  }

  bool HEEPnoIn() {
    if ( etSC < 35. ) return false;
    if ( !ecalDrivenSeed ) return false;
    if ( !(std::abs(dPhiIn) < 0.06) ) return false;
    if ( !(lostHits <= 1) ) return false;
    if ( std::abs(etaSC) < 1.4442 ) {
      bool HoEcut = ( full5x5_hOverE < (1./en + 0.05) );
      // bool E2x5cut = ( full5x5_E2x5/full5x5_E5x5 > 0.94 || full5x5_E1x5/full5x5_E5x5 > 0.83 );

      return ( HoEcut && (std::abs(dxy) < 0.02) ); // E2x5 excluded
    } else if ( std::abs(etaSC) > 1.566 && std::abs(etaSC) < 2.5 ) {
      bool HoEcut = ( full5x5_hOverE < (5./en + 0.05) );
      // bool E5x5cut = ( full5x5_E5x5 < 0.03 );

      return ( HoEcut && (std::abs(dxy) < 0.05) ); // E5x5 excluded
    }

    return false;
  }
};

class PhysicsMuon : public NtupleMuon
{
  public :
  bool Accep(Double_t Pt_, Double_t Eta_) {
    if( pt > Pt_ && fabs(eta) < Eta_) return true;
    else return false;
  }
};

class PhysicsGenParticle : public NtupleGenParticle
{
  public :
  TLorentzVector p4() {
    TLorentzVector p4_;
    p4_.SetPtEtaPhiM(pt, eta, phi, mass);
    return p4_;
  }
};

#endif
