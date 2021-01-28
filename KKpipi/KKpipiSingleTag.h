// Martin Duy Tat 28th January 2021
/**
 * KKpipiSingleTag is a class for a BOSS algorithm
 * It runs over \f$D^0\bar{D^0}\f$ data and saves all events with a single \f$D\to K^+K^-\pi^+\pi^-\f$ tag
 */

#ifndef KKPIPISINGLETAG
#define KKPIPISINGLETAG

// Gaudi
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
// BOSS
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "DTagTool/DTagTool.h"
#include "SimplePIDSvc/ISimplePIDSvc.h"
// STL
#include<string>

class KKpipiSingleTag: public Algorithm {
  public: 
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    KKpipiSingleTag(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~KKpipiSingleTag();
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
  private:
    /**
     * The NTuple itself that is filled
     */
    NTuple::tuple *m_tuple;
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
     * Number of particles at generator level
     */
    NTuple::Item<int> m_GeneratorNumberParticles;
    /**
     * Array of particle IDs of all particles in the decay chain at a generator level
     */
    Ntuple::Array<int> m_GeneratorPDGID;
    /**
     * Array of particle IDs of the mother particles at a generator level
     */
    NTuple::Array<int> m_MotherID;
    /**
     * True total momentum of the \f$D\f$ meson
     */
    NTuple::Item<double> m_TrueMomentum;
    /**
     * True transverse momentum of the \f$D\f$ meson
     */
    NTuple::Item<double> m_TruePT;
    /**
     * True azimuthal angle \f$\phi\f$
     */
    NTuple::Item<double> m_TruePhi;
    /**
     * True polar angle \f$\theta\f$ 
     */
    NTuple::Item<double> m_TrueTheta;
    /**
     * Charm content of \f$D\f$ meson
     */
    NTuple::Item<int> m_Charm;
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
     * \f$K^+\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_KPluspx;
    /**
     * \f$K^+\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_KPluspy;
    /**
     * \f$K^+\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_KPluspz;
    /**
     * \f$K^+\f$ energy
     */
    NTuple::Item<double> m_KPlusenergy;
    /**
     * \f$K^-\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_KMinuspx;
    /**
     * \f$K^-\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_KMinuspy;
    /**
     * \f$K^-\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_KMinuspz;
    /**
     * \f$K^-\f$ energy
     */
    NTuple::Item<double> m_KMinusenergy;
};

#endif
