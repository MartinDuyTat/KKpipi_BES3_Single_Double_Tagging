// Martin Duy Tat 28th January 2021, based on code by Yu Zhang
/**
 * KKpipiSingleTag is a class for a BOSS algorithm
 * It runs over \f$D^0\bar{D^0}\f$ data and saves all events with a single \f$D\to K^+K^-\pi^+\pi^-\f$ tag
 */

#ifndef KKPIPISINGLETAG
#define KKPIPISINGLETAG

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
// BOSS
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "DTagTool/DTagTool.h"
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
    /**
     * Helper function to fill in information about the \f$D\f$ candidate
     * This function should be the same for all tag modes
     * @param DTTool_iter Iterator pointing to single tag candidate
     */
    StatusCode AssignTagInfo(DTagToolIterator DTTool_iter);
    /**
     * Helper function to fill in information about the \f$D\f$ daughters
     * @param DTTool_iter Iterator pointing to single tag candidate
     * @param DTTool DTagTool object with all the tag information
     */
    StatusCode AssignKKpipiDaughterInfo(DTagToolIterator DTTool_iter, const DTagTool &DTTool);
    /**
     * Enumeration to label daughter particles in the order K+ K- pi+ pi-
     */
    enum DaughterParticle {KPLUS, KMINUS, PIPLUS, PIMINUS};
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
     * Number of particles at generator level
     */
    NTuple::Item<int> m_GeneratorNumberParticles;
    /**
     * Array of particle IDs of all particles in the decay chain at a generator level
     */
    NTuple::Array<int> m_GeneratorPDGID;
    /**
     * Array of particle IDs of the mother particles at a generator level
     */
    NTuple::Array<int> m_MotherID;
    /**
     * True total momentum of the \f$D\f$ meson
     */
    NTuple::Array<double> m_TrueMomentum;
    /**
     * True transverse momentum of the \f$D\f$ meson
     */
    NTuple::Array<double> m_TruePT;
    /**
     * True azimuthal angle \f$\phi\f$
     */
    NTuple::Array<double> m_TruePhi;
    /**
     * True polar angle \f$\theta\f$ 
     */
    NTuple::Array<double> m_TrueTheta;
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
    /**
     * Flag equal to 1 for success and 0 for fail in the Kalman fit of tracks
     */
    NTuple::Item<int> m_KalmanFitSuccess;
    /**
     * \f$\chi^2\f$ of Kalman fit
     */
    NTuple::Item<double> m_KalmanFitChi2;
    /**
     * Kalman fitted \f$\pi^+\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_PiPluspxKalmanFit;
    /**
     * Kalman fitted \f$\pi^+\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_PiPluspyKalmanFit;
    /**
     * Kalman fitted \f$\pi^+\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_PiPluspzKalmanFit;
    /**
     * Kalman fitted \f$\pi^+\f$ energy
     */
    NTuple::Item<double> m_PiPlusenergyKalmanFit;
    /**
     * Kalman fitted \f$\pi^-\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_PiMinuspxKalmanFit;
    /**
     * Kalman fitted \f$\pi^-\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_PiMinuspyKalmanFit;
    /**
     * Kalman fitted \f$\pi^-\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_PiMinuspzKalmanFit;
    /**
     * Kalman fitted \f$\pi^-\f$ energy
     */
    NTuple::Item<double> m_PiMinusenergyKalmanFit;
    /**
     * Kalman fitted \f$K^+\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_KPluspxKalmanFit;
    /**
     * Kalman fitted \f$K^+\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_KPluspyKalmanFit;
    /**
     * Kalman fitted \f$K^+\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_KPluspzKalmanFit;
    /**
     * Kalman fitted \f$K^+\f$ energy
     */
    NTuple::Item<double> m_KPlusenergyKalmanFit;
    /**
     * Kalman fitted \f$K^-\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_KMinuspxKalmanFit;
    /**
     * Kalman fitted \f$K^-\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_KMinuspyKalmanFit;
    /**
     * Kalman fitted \f$K^-\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_KMinuspzKalmanFit;
    /**
     * Kalman fitted \f$K^-\f$ energy
     */
    NTuple::Item<double> m_KMinusenergyKalmanFit;
};

#endif
