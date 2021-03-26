// Martin Duy Tat 26th March 2021
/**
 * KKpipiUtilities is a namespace with a few useful general purpose functions
 */

#ifndef KKPIPIUTILITIES
#define KKPIPIUTILITIES

// CLHEP
#include "CLHEP/Vector/LorentzVector.h"

namespace KKpipiUtilities {
  /**
   * Helper function to get four-momentum from shower information
   * @param Energy Energy of photon
   * @param Theta Polar angle of photon
   * @param Phi Azimuthal angle of photon
   */
  CLHEP::HepLorentzVector GetPhoton4Vector(double Energy, double Theta, double Phi);
}

#endif
