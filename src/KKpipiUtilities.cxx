// Martin Duy Tat 26th March 2021

// KKpipi
#include "KKpipi/KKpipiUtilities.h"
// Event information
#include "EvtRecEvent/EvtRecTrack.h"
#include "ExtEvent/RecExtTrack.h"
// ROOT
#include "TMath.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

CLHEP::HepLorentzVector KKpipiUtilities::GetPhoton4Vector(double Energy, double Theta, double Phi) {
  double Px = Energy*TMath::Sin(Theta)*TMath::Cos(Phi);
  double Py = Energy*TMath::Sin(Theta)*TMath::Sin(Phi);
  double Pz = Energy*TMath::Cos(Theta);
  return CLHEP::HepLorentzVector(Px, Py, Pz, Energy);
}

bool KKpipiUtilities::GetPhotonAngularSeparation(const CLHEP::Hep3Vector &EMCPosition, EvtRecTrackIterator Track_iter_begin, int TotalCharged, double &Angle, double &Theta, double &Phi) {
  Theta = 2*TMath::Pi();
  Phi = 2*TMath::Pi();
  Angle = 2*TMath::Pi();
  // Loop over all charged tracks
  for(int j = 0; j < TotalCharged; j++) {
    EvtRecTrackIterator Track_iter = Track_iter_begin + j;
    // Check if track is valid
    if(!(*Track_iter)->isExtTrackValid()) {
      continue;
    }
    // Get track extrapolated to the EMC
    RecExtTrack *ExternalTrack = (*Track_iter)->extTrack();
    // Check if external track is valid
    if(ExternalTrack->emcVolumeNumber() == -1) {
      continue;
    }
    // Get position of external track
    CLHEP::Hep3Vector ExternalPosition = ExternalTrack->emcPosition();
    // Find angle between track and shower
    double DeltaAngle = ExternalPosition.angle(EMCPosition);
    // Find polar angle between track and shower
    double DeltaTheta = ExternalPosition.theta() - EMCPosition.theta();
    // Find azimuthal angle between track and shower
    double DeltaPhi = ExternalPosition.deltaPhi(EMCPosition);
    if(DeltaAngle < Angle) {
      Theta = DeltaTheta;
      Phi = DeltaPhi;
      Angle = DeltaAngle;
    }
  }
  if(Angle == 2*TMath::Pi()) {
    return false;
  } else {
    return true;
  }
}
