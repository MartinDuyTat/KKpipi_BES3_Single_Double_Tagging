// Martin Duy Tat 28th January 2021, based on code by Yu Zhang
/**
 * FindMCInfo is a class for finding the MC truth info of an event
 */

#ifndef FINDMCINFO
#define FINDMCINFO

// Gaudi
#include "GaudiKernel/StatusCode.h"
// STL
#include<vector>

class FindMCInfo {
  public: 
    /**
     * Default constructor which initializes all variables to zero
     */
    FindMCInfo();
    /**
     * Trivial destructor
     */
    ~FindMCInfo();
    /**
     * Helper function to calculate the MC truth information
     */
    StatusCode CalculateMCInfo();
    /**
     * Get the number of particles
     */
    int GetNumberParticles();
    /**
     * Get vector of PDG IDs
     */
    std::vector<int> GetpdgID();
    /**
     * Get vector of mother indices
     */
    std::vector<int> GetMotherIndex();
    /**
     * Get MC mode label
     */
    int GetMCmode();
    /**
     * Get vector of true momenta in the \f$x\f$ direction
     */
    std::vector<int> GetTruePx();
    /**
     * Get vector of true momenta in the \f$y\f$ direction
     */
    std::vector<int> GetTruePy();
    /**
     * Get vector of true momenta in the \f$y\f$ direction
     */
    std::vector<int> GetTruePz();
    /**
     * Get vector of true energies
     */
    std::vector<int> GetTrueEnergy();
  private:
    /**
     * Number of particles in the decay chain
     * For example, \f$D^0\to K^-\pi^+\f$ counts as \f$3\f$ particles
     */
    int m_NumberParticles;
    /**
     * Array of particle IDs of all particles in the decay chain
     */
    std::vector<int> m_pdgID;
    /**
     * Vector of indices referring to the particle mother
     * For example, an index 0 means that ```m_pdgID[0]``` is the mother
     */
    std::vector<int> m_MotherIndex;
    /**
     * Generator label of the decay mode
     */
    int m_MCmode;
    /**
     * True x momentum of the \f$D\f$ meson
     */
    std::vector<double> m_TruePx;
    /**
     * True y momentum of the \f$D\f$ meson
     */
    std::vector<double> m_TruePy;
    /**
     * True z momentum of the \f$D\f$ meson
     */
    std::vector<double> m_TruePz;
    /**
     * True total energy of the \f$D\f$ meson
     */
    std::vector<double> m_TrueEnergy;
};

#endif
