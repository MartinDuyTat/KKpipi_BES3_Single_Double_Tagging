// Martin Duy Tat 26th March 2021
/**
 * KKpipiUtilities is a namespace with a few useful general purpose functions
 */

#ifndef KKPIPIUTILITIES
#define KKPIPIUTILITIES

// Event information
#include "EvtRecEvent/EvtRecTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace KKpipiUtilities {
  /**
   * Helper function to get four-momentum from shower information
   * @param Energy Energy of photon
   * @param Theta Polar angle of photon
   * @param Phi Azimuthal angle of photon
   */
  CLHEP::HepLorentzVector GetPhoton4Vector(double Energy, double Theta, double Phi);
  /**
   * Helper function to calculate the angular separation between the EMC shower and the nearest charged track
   * @param EMCPosition Position of the EMC shower
   * @param TrackIter Iterator of the first charged track
   * @param TotalCharged Total number of charged tracks
   * @param Angle Output, the angular separation between the shower and the nearest charged track
   * @param Theta Output, the polar angle separation between the shower and the nearest charged track
   * @param Phi Output, the azimuthal angle separation between the shower and the nearest charged track
   * @return Returns true if the calculation was successful
   */
  bool GetPhotonAngularSeparation(const CLHEP::Hep3Vector &EMCPosition, EvtRecTrackIterator Track_iter_begin, int TotalCharged, double &Angle, double &Theta, double &Phi);
}

#endif
