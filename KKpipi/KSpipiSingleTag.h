// Martin Duy Tat 12th March 2021
/**
 * KSpipiSingleTag is a class for a BOSS algorithm
 * It runs over \f$D^0\bar{D^0}\f$ data and saves all events with a single \f$D\to K_S^0\pi^+\pi^-\f$ tag
 */

#ifndef KSPIPISINGLETAG
#define KSPIPISINGLETAG

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// STL
#include<string>

class KSpipiSingleTag: public Algorithm {
  public: 
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    KSpipiSingleTag(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~KSpipiSingleTag();
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
     * Number of particles in the decay chain without resonances
     */
    NTuple::Item<int> m_NumberParticlesStripped;
    /**
     * Array of particle IDs of all particles in the decay chain without resonances
     */
    NTuple::Array<int> m_pdgIDStripped;
    /**
     * Array of indices referring to the particle mother without resonances
     */
    NTuple::Array<int> m_MotherIndexStripped;
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
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tracks
     */
    NTuple::Item<int> m_KSFitSuccess;
    /**
     * The \f$K_S\f$ decay length, from VeeVertexAlg
     */
    NTuple::Item<double> m_DecayLengthVeeVertex;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_Chi2VeeVertex;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    NTuple::Item<double> m_KSMassVeeVertex;
    /**
     * The \f$K_S^0\f$ decay length, from fit
     */
    NTuple::Item<double> m_DecayLengthFit;
    /**
     * The \f$K_S^0\f$ decay length error, from fit
     */
    NTuple::Item<double> m_DecayLengthErrorFit;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    NTuple::Item<double> m_Chi2Fit;
    /**
     * \f$K_S^0\f$ momentum along \f$x\f$ after vertex fit
     */
    NTuple::Item<double> m_KSpx;
    /**
     * \f$K_S^0\f$ momentum along \f$y\f$ after vertex fit
     */
    NTuple::Item<double> m_KSpy;
    /**
     * \f$K_S^0\f$ momentum along \f$z\f$ after vertex fit
     */
    NTuple::Item<double> m_KSpz;
    /**
     * \f$K_S^0\f$ energy after vertex fit
     */
    NTuple::Item<double> m_KSenergy;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$x\f$ from the MDC track
     */
    NTuple::Item<double> m_KSPiPluspx;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$y\f$ from the MDC track
     */
    NTuple::Item<double> m_KSPiPluspy;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$z\f$ from the MDC track
     */
    NTuple::Item<double> m_KSPiPluspz;
    /**
     * The \f$\pi^+\f$ daughter energy from the MDC track
     */
    NTuple::Item<double> m_KSPiPlusenergy;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$x\f$ from the MDC track
     */
    NTuple::Item<double> m_KSPiMinuspx;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$y\f$ from the MDC track
     */
    NTuple::Item<double> m_KSPiMinuspy;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$z\f$ from the MDC track
     */
    NTuple::Item<double> m_KSPiMinuspz;
    /**
     * The \f$\pi^-\f$ daughter energy from the MDC track
     */
    NTuple::Item<double> m_KSPiMinusenergy;
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
     * Flag equal to 1 for success and 0 for fail in the Kalman fit of tag tracks
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
     * Kalman fitted \f$K_S^0\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_KSpxKalmanFit;
    /**
     * Kalman fitted \f$K_S^0\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_KSpyKalmanFit;
    /**
     * Kalman fitted \f$K_S^0\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_KSpzKalmanFit;
    /**
     * Kalman fitted \f$K_S^0\f$ energy
     */
    NTuple::Item<double> m_KSenergyKalmanFit;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tag tracks
     */
    NTuple::Item<int> m_pipiKSFitSuccess;
    /**
     * \f$K_S\f$ decay length, from VeeVertexAlg
     */
    NTuple::Item<double> m_pipiDecayLengthVeeVertex;
    /**
     * \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_pipiChi2VeeVertex;
    /**
     * \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    NTuple::Item<double> m_pipiKSMassVeeVertex;
    /**
     * \f$K_S^0\f$ decay length, from fit
     */
    NTuple::Item<double> m_pipiDecayLengthFit;
    /**
     * \f$K_S^0\f$ decay length error, from fit
     */
    NTuple::Item<double> m_pipiDecayLengthErrorFit;
    /**
     * \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    NTuple::Item<double> m_pipiChi2Fit;
    /**
     * Equal to 1 if the daughter tracks are from the same \f$D\f$ meson
     */
    NTuple::Item<int> m_IsSameDMother;
    /**
     * Equal to 1 if the daughter tracks are assigned a PID matching that of the MC truth
     */
    NTuple::Item<int> m_PIDTrue;
    /**
     * The \f$\pi^+\f$ from \f$K_S^0\f$ true PID
     */
    NTuple::Item<int> m_KSPiPlusTrueID;
    /**
     * The \f$\pi^-\f$ from \f$K_S^0\f$ true PID
     */
    NTuple::Item<int> m_KSPiMinusTrueID;
    /**
     * The \f$\pi^+\f$ true PID
     */
    NTuple::Item<int> m_PiPlusTrueID;
    /**
     * The \f$\pi^-\f$ true PID
     */
    NTuple::Item<int> m_PiMinusTrueID;
    /**
     * The \f$\pi^+\f$ from \f$K_S^0\f$ true mother PID
     */
    NTuple::Item<int> m_KSPiPlusMotherTrueID;
    /**
     * The \f$\pi^-\f$ from \f$K_S^0\f$ true mother PID
     */
    NTuple::Item<int> m_KSPiMinusMotherTrueID;
};

#endif
