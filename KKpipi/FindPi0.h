// Martin Duy Tat 12th February 2021, based on code by Yu Zhang
/**
 * FindPi0 is class for finding \f$\pi^0\f$ candidates and returning the unconstrained and kinematically constrained four-momenta of the photons
 * By default one \f$\pi^0\f$ is assumed, unless otherwise is specified in the constructor and the getters
 */

#ifndef FINDPI0
#define FINDPI0
// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include <vector>

class FindPi0 {
  public: 
    /**
     * Default constructor, initializes all properties to zero and the number of \f$\pi^0\f$ to one
     * @param npi0 Number of \f$\pi^0\f$ to look for
     */
    FindPi0(int npi0 = 1);
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
     * @param pi0_index 0 for first \f$\pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    double GetHighEPhotonP(int i, int pi0_index = 0) const;
    /** 
     * Get the high energy photon unconstrained four-momentum
     * @param i Momentum component
     * @param pi0_index 0 for first \f$\pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    double GetLowEPhotonP(int i, int pi0_index = 0) const;
    /** 
     * Get the high energy photon kinematically constrained four-momentum
     * @param i Momentum component
     * @param pi0_index 0 for first \f$\pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    double GetHighEPhotonPConstrained(int i, int pi0_index = 0) const;
    /** 
     * Get the low energy photon kinematically constrained four-momentum
     * @param i Momentum component
     * @param pi0_index 0 for first \f$\pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    double GetLowEPhotonPConstrained(int i, int pi0_index = 0) const;
    /**
     * Get the \f$\chi^2\f$ of the kinematic fit
     * @param pi0_index 0 for first \f$\pi^0\f$, 1 for second \f$\pi^0\f$...
     */
    double GetChi2Fit(int pi0_index = 0) const;
  private:
    /**
     * The high energy photon unconstrained four-momentum
     */
    std::vector<CLHEP::HepLorentzVector> m_HighEPhotonP;
    /**
     * The low energy photon unconstrained four-momentum
     */
    std::vector<CLHEP::HepLorentzVector> m_LowEPhotonP;
    /**
     * The high energy photon kinematically constrained four-momentum
     */
    std::vector<CLHEP::HepLorentzVector> m_HighEPhotonPConstrained;
    /**
     * The high energy photon kinematically constrained four-momentum
     */
    std::vector<CLHEP::HepLorentzVector> m_LowEPhotonPConstrained;
    /**
     * Kinematic fit \f$\chi^2\f$
     */
    std::vector<double> m_Chi2Fit;
    /**
     * Number of \f$\pi^0\f$ to look for
     */
    int m_npi0;
};

#endif
