// Martin Duy Tat 26th March 2021

// KKpipi
#include "KKpipi/KKpipiUtilities.h"
// ROOT
#include "TMath.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"

CLHEP::HepLorentzVector KKpipiUtilities::GetPhoton4Vector(double Energy, double Theta, double Phi) {
  double Px = Energy*TMath::Sin(Theta)*TMath::Cos(Phi);
  double Py = Energy*TMath::Sin(Theta)*TMath::Sin(Phi);
  double Pz = Energy*TMath::Cos(Theta);
  return CLHEP::HepLorentzVector(Px, Py, Pz, Energy);
}
