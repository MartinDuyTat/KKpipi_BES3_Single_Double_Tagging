// Martin Duy Tat 12th February 2021, based on code by Yu Zhang
/**
 * FindPi0 is class for finding \f$\pi^0\f$ candidates and returning the unconstrained and kinematically constrained four-momenta of the photons
 */

#ifndef FINDPI0
#define FINDPI0
// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"

class FindPi0 {
  public: 
    /**
     * Default constructor, initializes everything to zero
     */
    FindPi0();
    /**
     * Trivial destructor
     */
    ~FindPi0();
    /**
     * Helper function to get four-momentum from shower information
     * @param Energy Energy of photon
     * @param Theta Polar angle of photon
     * @param Phi Azimuthal angle of photon
     */
    CLHEP::HepLorentzVector GetPhoton4Vector(double Energy, double Theta, double Phi) const;
    /**
     * Start looking for \f$\pi^0\f$ in the event
     * @param DTTool_iter DTagTool iterator pointing to the event with the tag
     * @param PiTrackIndex List of length 2 with track indices to the two pions in the event
     */
    StatusCode findPi0(DTagToolIterator &DTTool_iter, DTagTool &DTTool);
    /** 
     * Get the high energy photon unconstrained four-momentum
     * @param i Momentum component
     */
    double GetHighEPhotonP(int i) const;
    /** 
     * Get the high energy photon unconstrained four-momentum
     * @param i Momentum component
     */
    double GetLowEPhotonP(int i) const;
    /** 
     * Get the high energy photon kinematically constrained four-momentum
     * @param i Momentum component
     */
    double GetHighEPhotonPConstrained(int i) const;
    /** 
     * Get the low energy photon kinematically constrained four-momentum
     * @param i Momentum component
     */
    double GetLowEPhotonPConstrained(int i) const;
    /**
     * Get the \f$\chi^2\f$ of the kinematic fit
     */
    double GetChi2Fit() const;
  private:
    /**
     * The high energy photon unconstrained four-momentum
     */
    CLHEP::HepLorentzVector m_HighEPhotonP;
    /**
     * The low energy photon unconstrained four-momentum
     */
    CLHEP::HepLorentzVector m_LowEPhotonP;
    /**
     * The high energy photon kinematically constrained four-momentum
     */
    CLHEP::HepLorentzVector m_HighEPhotonPConstrained;
    /**
     * The high energy photon kinematically constrained four-momentum
     */
    CLHEP::HepLorentzVector m_LowEPhotonPConstrained;
    /**
     * Kinematic fit \f$\chi^2\f$
     */
    double m_Chi2Fit;
};

#endif
