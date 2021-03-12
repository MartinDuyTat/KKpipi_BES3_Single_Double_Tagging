// Martin Duy Tat 12th March 2021
/**
 * pipipi0SingleTag is a class for a BOSS algorithm
 * It runs over \f$D^0\bar{D^0}\f$ data and saves all events with a single \f$D\to \pi^+\pi^-\pi^0\f$ tag
 */

#ifndef PIPIPI0SINGLETAG
#define PIPIPI0SINGLETAG

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// STL
#include<string>

class pipipi0SingleTag: public Algorithm {
  public: 
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    pipipi0SingleTag(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~pipipi0SingleTag();
    /**
     * This function runs when algorithm is initialized
     */
    StatusCode initialize();
    /**
     * Execution of the algorithm
     */
    StatusCode execute();
    /**
     * This function runs when algorithm is finalized
     */
    StatusCode finalize();
    /**
     * Helper function to fill in information about the tag mode
     * @param DTTool_iter Iterator pointing to single tag candidate
     * @param DTTool DTagTool object with all the event information
     */
    StatusCode FillTuple(DTagToolIterator DTTool_iter, DTagTool &DTTool);
  private:
    /**
     * Dummy variable, placeholder for more important properties to be added later
     */
    int m_dummy;
    /**
     * The NTuple itself that is filled
     */
    NTuple::Tuple *m_tuple;
    /**
     * Run number (negative for MC event)
     */
    NTuple::Item<int> m_RunNumber;
    /**
     * Event number
     */
    NTuple::Item<int> m_EventNumber;
    /**
     * Number of particles in the decay chain
     * For example, \f$D^0\to K^-\pi^+\f$ counts as \f$3\f$ particles
     */
    NTuple::Item<int> m_NumberParticles;
    /**
     * Array of particle IDs of all particles in the decay chain
     */
    NTuple::Array<int> m_pdgID;
    /**
     * Array of indices referring to the particle mother
     * For example, an index 0 means that ```m_pdgID[0]``` is the mother
     */
    NTuple::Array<int> m_MotherIndex;
    /**
     * Generator label of the decay mode
     */
    NTuple::Item<int> m_MCmode;
    /**
     * Array of true x momenta
     */
    NTuple::Array<double> m_TruePx;
    /**
     * Array of true y momenta
     */
    NTuple::Array<double> m_TruePy;
    /**
     * Array of true z momenta
     */
    NTuple::Array<double> m_TruePz;
    /**
     * Array of true energies
     */
    NTuple::Array<double> m_TrueEnergy;
    /**
     * Invariant mass of \f$D\f$ meson
     */
    NTuple::Item<double> m_DMass;
    /**
     * Beam constrained mass
     */
    NTuple::Item<double> m_MBC;
    /**
     * \f$E_D - E_{\text{beam}}\f$
     */
    NTuple::Item<double> m_DeltaE;
    /**
     * Beam energy
     */
    NTuple::Item<double> m_BeamE;
    /**
     * \f$D\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_Dpx;
    /**
     * \f$D\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_Dpy;
    /**
     * \f$D\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_Dpz;
    /**
     * \f$D\f$ energy
     */
    NTuple::Item<double> m_Denergy;
    /**
     * \f$\pi^+\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_PiPluspx;
    /**
     * \f$\pi^+\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_PiPluspy;
    /**
     * \f$\pi^+\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_PiPluspz;
    /**
     * \f$\pi^+\f$ energy
     */
    NTuple::Item<double> m_PiPlusenergy;
    /**
     * \f$\pi^-\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_PiMinuspx;
    /**
     * \f$\pi^-\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_PiMinuspy;
    /**
     * \f$\pi^-\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_PiMinuspz;
    /**
     * \f$\pi^-\f$ energy
     */
    NTuple::Item<double> m_PiMinusenergy;
    /**
     * High energy photon from \f$\pi^0\f$ unconstrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_HighEPi0px;
    /**
     * High energy photon from \f$\pi^0\f$ unconstrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_HighEPi0py;
    /**
     * High energy photon from \f$\pi^0\f$ unconstrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_HighEPi0pz;
    /**
     * High energy photon from \f$\pi^0\f$ unconstrained energy
     */
    NTuple::Item<double> m_HighEPi0energy;
    /**
     * Low energy photon from \f$\pi^0\f$ unconstrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_LowEPi0px;
    /**
     * Low energy photon from \f$\pi^0\f$ unconstrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_LowEPi0py;
    /**
     * Low energy photon from \f$\pi^0\f$ unconstrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_LowEPi0pz;
    /**
     * Low energy photon from \f$\pi^0\f$ unconstrained energy
     */
    NTuple::Item<double> m_LowEPi0energy;
    /**
     * High energy photon from \f$\pi^0\f$ constrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_HighEPi0Constrainedpx;
    /**
     * High energy photon from \f$\pi^0\f$ constrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_HighEPi0Constrainedpy;
    /**
     * High energy photon from \f$\pi^0\f$ constrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_HighEPi0Constrainedpz;
    /**
     * High energy photon from \f$\pi^0\f$ constrained energy
     */
    NTuple::Item<double> m_HighEPi0Constrainedenergy;
    /**
     * Low energy photon from \f$\pi^0\f$ constrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_LowEPi0Constrainedpx;
    /**
     * Low energy photon from \f$\pi^0\f$ constrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_LowEPi0Constrainedpy;
    /**
     * Low energy photon from \f$\pi^0\f$ constrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_LowEPi0Constrainedpz;
    /**
     * Low energy photon from \f$\pi^0\f$ constrained energy
     */
    NTuple::Item<double> m_LowEPi0Constrainedenergy;
    /**
     * Tag \f$\pi^0\f$ kinematic fit \f$\chi^2\f$
     */
    NTuple::Item<double> m_Pi0Chi2Fit;
};

#endif