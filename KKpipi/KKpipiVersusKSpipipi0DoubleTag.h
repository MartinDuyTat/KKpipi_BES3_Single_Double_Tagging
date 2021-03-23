// Martin Duy Tat 4th March 2021
/**
 * KKpipiVersusKSpipipi0DoubleTag is a class for a BOSS algorithm
 * It runs over \f$D^0\bar{D^0}\f$ data and saves all events with a double \f$D\to K^+K^-\pi^+\pi^-\f$ vs \f$D^0\to K_S^0\pi^+\pi^-\pi^0\f$ tag
 * It also runs a fit for the decay \f$K_S^0\to\pi^+\pi^-\f$ by refitting the primary and secondary vertex in the class FindKS, from this the flight significance is used to eliminate peaking background
 */

#ifndef KKPIPIVERSUSKSPIPIPI0DOUBLETAG
#define KKPIPIVERSUSKSPIPIPI0DOUBLETAG

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// STL
#include<string>

class KKpipiVersusKSpipipi0DoubleTag: public Algorithm {
  public: 
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    KKpipiVersusKSpipipi0DoubleTag(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~KKpipiVersusKSpipipi0DoubleTag();
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
     * @param DTTool_Signal_iter Iterator pointing to signal candidate
     * @param DTTool_Signal_iter Iterator pointing to tag candidate
     * @param DTTool DTagTool object with all the event information
     */
    StatusCode FillTuple(DTagToolIterator DTTool_Signal_iter, DTagToolIterator DTTool_Tag_iter, DTagTool &DTTool);
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
     * Invariant mass of signal \f$D\f$ meson
     */
    NTuple::Item<double> m_SignalDMass;
    /**
     * Beam constrained mass of signal
     */
    NTuple::Item<double> m_SignalMBC;
    /**
     * \f$E_D - E_{\text{beam}}\f$ of signal
     */
    NTuple::Item<double> m_SignalDeltaE;
    /**
     * Beam energy of signal
     */
    NTuple::Item<double> m_SignalBeamE;
    /**
     * Invariant mass of tag \f$D\f$ meson
     */
    NTuple::Item<double> m_TagDMass;
    /**
     * Beam constrained mass of tag
     */
    NTuple::Item<double> m_TagMBC;
    /**
     * \f$E_D - E_{\text{beam}}\f$ of tag
     */
    NTuple::Item<double> m_TagDeltaE;
    /**
     * Beam energy of tag
     */
    NTuple::Item<double> m_TagBeamE;
    /**
     * Signal \f$D\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalDpx;
    /**
     * Signal \f$D\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalDpy;
    /**
     * Signal \f$D\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalDpz;
    /**
     * Signal \f$D\f$ energy
     */
    NTuple::Item<double> m_SignalDenergy;
    /**
     * Tag \f$D\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagDpx;
    /**
     * Tag \f$D\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagDpy;
    /**
     * Tag \f$D\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagDpz;
    /**
     * Tag \f$D\f$ energy
     */
    NTuple::Item<double> m_TagDenergy;
    /**
     * Signal \f$\pi^+\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalPiPluspx;
    /**
     * Signal \f$\pi^+\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalPiPluspy;
    /**
     * Signal \f$\pi^+\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalPiPluspz;
    /**
     * Signal \f$\pi^+\f$ energy
     */
    NTuple::Item<double> m_SignalPiPlusenergy;
    /**
     * Signal \f$\pi^-\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalPiMinuspx;
    /**
     * Signal \f$\pi^-\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalPiMinuspy;
    /**
     * Signal \f$\pi^-\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalPiMinuspz;
    /**
     * Signal \f$\pi^-\f$ energy
     */
    NTuple::Item<double> m_SignalPiMinusenergy;
    /**
     * Signal \f$K^+\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalKPluspx;
    /**
     * Signal \f$K^+\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalKPluspy;
    /**
     * Signal \f$K^+\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalKPluspz;
    /**
     * Signal \f$K^+\f$ energy
     */
    NTuple::Item<double> m_SignalKPlusenergy;
    /**
     * Signal \f$K^-\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalKMinuspx;
    /**
     * Signal \f$K^-\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalKMinuspy;
    /**
     * Signal \f$K^-\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalKMinuspz;
    /**
     * Signal \f$K^-\f$ energy
     */
    NTuple::Item<double> m_SignalKMinusenergy;
    /**
     * Flag equal to 1 for success and 0 for fail in the Kalman fit of signal tracks
     */
    NTuple::Item<int> m_SignalKalmanFitSuccess;
    /**
     * Signal \f$\chi^2\f$ of Kalman fit
     */
    NTuple::Item<double> m_SignalKalmanFitChi2;
    /**
     * Signal Kalman fitted \f$\pi^+\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalPiPluspxKalmanFit;
    /**
     * Signal Kalman fitted \f$\pi^+\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalPiPluspyKalmanFit;
    /**
     * Signal Kalman fitted \f$\pi^+\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalPiPluspzKalmanFit;
    /**
     * Signal Kalman fitted \f$\pi^+\f$ energy
     */
    NTuple::Item<double> m_SignalPiPlusenergyKalmanFit;
    /**
     * Signal Kalman fitted \f$\pi^-\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalPiMinuspxKalmanFit;
    /**
     * Signal Kalman fitted \f$\pi^-\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalPiMinuspyKalmanFit;
    /**
     * Signal Kalman fitted \f$\pi^-\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalPiMinuspzKalmanFit;
    /**
     * Signal Kalman fitted \f$\pi^-\f$ energy
     */
    NTuple::Item<double> m_SignalPiMinusenergyKalmanFit;
    /**
     * Signal Kalman fitted \f$K^+\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalKPluspxKalmanFit;
    /**
     * Signal Kalman fitted \f$K^+\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalKPluspyKalmanFit;
    /**
     * Signal Kalman fitted \f$K^+\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalKPluspzKalmanFit;
    /**
     * Signal Kalman fitted \f$K^+\f$ energy
     */
    NTuple::Item<double> m_SignalKPlusenergyKalmanFit;
    /**
     * Signal Kalman fitted \f$K^-\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_SignalKMinuspxKalmanFit;
    /**
     * Signal Kalman fitted \f$K^-\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_SignalKMinuspyKalmanFit;
    /**
     * Signal Kalman fitted \f$K^-\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_SignalKMinuspzKalmanFit;
    /**
     * SignalKalman fitted \f$K^-\f$ energy
     */
    NTuple::Item<double> m_SignalKMinusenergyKalmanFit;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of signal tracks
     */
    NTuple::Item<int> m_SignalKSFitSuccess;
    /**
     * The signal \f$K_S\f$ decay length, from VeeVertexAlg
     */
    NTuple::Item<double> m_SignalDecayLengthVeeVertex;
    /**
     * The signal \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_SignalChi2VeeVertex;
    /**
     * The signal \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    NTuple::Item<double> m_SignalKSMassVeeVertex;
    /**
     * The signal \f$K_S^0\f$ decay length, from fit
     */
    NTuple::Item<double> m_SignalDecayLengthFit;
    /**
     * The signal \f$K_S^0\f$ decay length error, from fit
     */
    NTuple::Item<double> m_SignalDecayLengthErrorFit;
    /**
     * The signal \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    NTuple::Item<double> m_SignalChi2Fit;
    /**
     * The signal \f$K_S^0\f$ mass, from fit
     */
    NTuple::Item<double> m_SignalKSMassFit;
    /**
     * Equal to 1 if the signal daughter tracks are from the same \f$D\f$ meson
     */
    NTuple::Item<int> m_SignalIsSameDMother;
    /**
     * Equal to 1 if the signal daughter tracks are assigned a PID matching that of the MC truth
     */
    NTuple::Item<int> m_SignalPIDTrue;
    /**
     * The signal \f$K^+\f$ true PID
     */
    NTuple::Item<int> m_SignalKPlusTrueID;
    /**
     * The signal \f$K^-\f$ true PID
     */
    NTuple::Item<int> m_SignalKMinusTrueID;
    /**
     * The signal \f$\pi^+\f$ true PID
     */
    NTuple::Item<int> m_SignalPiPlusTrueID;
    /**
     * The signal \f$\pi^-\f$ true PID
     */
    NTuple::Item<int> m_SignalPiMinusTrueID;
    /**
     * The \f$K_S\f$ decay length, from VeeVertexAlg
     */
    NTuple::Item<double> m_TagDecayLengthVeeVertex;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_TagChi2VeeVertex;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    NTuple::Item<double> m_TagKSMassVeeVertex;
    /**
     * The \f$K_S^0\f$ decay length, from fit
     */
    NTuple::Item<double> m_TagDecayLengthFit;
    /**
     * The \f$K_S^0\f$ decay length error, from fit
     */
    NTuple::Item<double> m_TagDecayLengthErrorFit;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    NTuple::Item<double> m_TagChi2Fit;
    /**
     * The \f$K_S^0\f$ mass, from fit
     */
    NTuple::Item<double> m_TagKSMassFit;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$x\f$ from the MDC track
     */
    NTuple::Item<double> m_TagKSPiPluspx;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$y\f$ from the MDC track
     */
    NTuple::Item<double> m_TagKSPiPluspy;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$z\f$ from the MDC track
     */
    NTuple::Item<double> m_TagKSPiPluspz;
    /**
     * The \f$\pi^+\f$ daughter energy from the MDC track
     */
    NTuple::Item<double> m_TagKSPiPlusenergy;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$x\f$ from the MDC track
     */
    NTuple::Item<double> m_TagKSPiMinuspx;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$y\f$ from the MDC track
     */
    NTuple::Item<double> m_TagKSPiMinuspy;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$z\f$ from the MDC track
     */
    NTuple::Item<double> m_TagKSPiMinuspz;
    /**
     * The \f$\pi^-\f$ daughter energy from the MDC track
     */
    NTuple::Item<double> m_TagKSPiMinusenergy;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$x\f$ from the MDC track
     */
    NTuple::Item<double> m_TagKSPiPluspxFit;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$y\f$ after vertex fit
     */
    NTuple::Item<double> m_TagKSPiPluspyFit;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$z\f$ after vertex fit
     */
    NTuple::Item<double> m_TagKSPiPluspzFit;
    /**
     * The \f$\pi^+\f$ daughter energy after vertex fit
     */
    NTuple::Item<double> m_TagKSPiPlusenergyFit;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$x\f$ after vertex fit
     */
    NTuple::Item<double> m_TagKSPiMinuspxFit;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$y\f$ after vertex fit
     */
    NTuple::Item<double> m_TagKSPiMinuspyFit;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$z\f$ after vertex fit
     */
    NTuple::Item<double> m_TagKSPiMinuspzFit;
    /**
     * The \f$\pi^-\f$ daughter energy after vertex fit
     */
    NTuple::Item<double> m_TagKSPiMinusenergyFit;
    /**
     * Tag high energy photon from \f$\pi^0\f$ unconstrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagHighEPi0px;
    /**
     * Tag high energy photon from \f$\pi^0\f$ unconstrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagHighEPi0py;
    /**
     * Tag high energy photon from \f$\pi^0\f$ unconstrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagHighEPi0pz;
    /**
     * Tag high energy photon from \f$\pi^0\f$ unconstrained energy
     */
    NTuple::Item<double> m_TagHighEPi0energy;
    /**
     * Tag low energy photon from \f$\pi^0\f$ unconstrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagLowEPi0px;
    /**
     * Tag low energy photon from \f$\pi^0\f$ unconstrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagLowEPi0py;
    /**
     * Tag low energy photon from \f$\pi^0\f$ unconstrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagLowEPi0pz;
    /**
     * Tag low energy photon from \f$\pi^0\f$ unconstrained energy
     */
    NTuple::Item<double> m_TagLowEPi0energy;
    /**
     * Tag high energy photon from \f$\pi^0\f$ constrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagHighEPi0Constrainedpx;
    /**
     * Tag high energy photon from \f$\pi^0\f$ constrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagHighEPi0Constrainedpy;
    /**
     * Tag high energy photon from \f$\pi^0\f$ constrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagHighEPi0Constrainedpz;
    /**
     * Tag high energy photon from \f$\pi^0\f$ constrained energy
     */
    NTuple::Item<double> m_TagHighEPi0Constrainedenergy;
    /**
     * Tag low energy photon from \f$\pi^0\f$ constrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagLowEPi0Constrainedpx;
    /**
     * Tag low energy photon from \f$\pi^0\f$ constrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagLowEPi0Constrainedpy;
    /**
     * Tag low energy photon from \f$\pi^0\f$ constrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagLowEPi0Constrainedpz;
    /**
     * Tag low energy photon from \f$\pi^0\f$ constrained energy
     */
    NTuple::Item<double> m_TagLowEPi0Constrainedenergy;
    /**
     * Tag \f$\pi^0\f$ kinematic fit \f$\chi^2\f$
     */
    NTuple::Item<double> m_Pi0Chi2Fit;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of \f$\pi^+\pi^-\f$ tag tracks
     */
    NTuple::Item<int> m_TagpipiKSFitSuccess;
    /**
     * The tag \f$K_S\f$ decay length from \f$\pi^+\pi^-\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_TagpipiDecayLengthVeeVertex;
    /**
     * The tag \f$K_S^0\f$ \f$\chi^2\f$, from \f$\pi^+\pi^-\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_TagpipiChi2VeeVertex;
    /**
     * The tag \f$K_S^0\f$ mass, from \f$\pi^+\pi^-\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_TagpipiKSMassVeeVertex;
    /**
     * The tag \f$K_S^0\f$ decay length, from \f$\pi^+\pi^-\f$, from fit
     */
    NTuple::Item<double> m_TagpipiDecayLengthFit;
    /**
     * The tag \f$K_S^0\f$ decay length error, from \f$\pi^+\pi^-\f$, from fit
     */
    NTuple::Item<double> m_TagpipiDecayLengthErrorFit;
    /**
     * The tag \f$K_S^0\f$ \f$\chi^2\f$, from \f$\pi^+\pi^-\f$, from fit of primary vertex
     */
    NTuple::Item<double> m_TagpipiChi2Fit;
    /**
     * The tag \f$K_S^0\f$ mass, from \f$\pi^+\pi^-\f$, from fit
     */
    NTuple::Item<double> m_TagpipiKSMassFit;
    /**
     * Tag \f$\pi^+\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagPiPluspx;
    /**
     * Tag \f$\pi^+\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagPiPluspy;
    /**
     * Tag \f$\pi^+\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagPiPluspz;
    /**
     * Tag \f$\pi^=\f$ energy
     */
    NTuple::Item<double> m_TagPiPlusenergy;
    /**
     * Tag \f$\pi^-\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagPiMinuspx;
    /**
     * Tag \f$\pi^-\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagPiMinuspy;
    /**
     * Tag \f$\pi^-\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagPiMinuspz;
    /**
     * Tag \f$\pi^-\f$ energy
     */
    NTuple::Item<double> m_TagPiMinusenergy;
    /**
     * Equal to 1 if the tag daughter tracks are from the same \f$D\f$ meson
     */
    NTuple::Item<int> m_TagIsSameDMother;
    /**
     * Equal to 1 if the tag daughter tracks are assigned a PID matching that of the MC truth
     */
    NTuple::Item<int> m_TagPIDTrue;
    /**
     * The tag \f$\pi^+\f$ from \f$K_S^0\f$ true PID
     */
    NTuple::Item<int> m_TagKSPiPlusTrueID;
    /**
     * The tag \f$\pi^-\f$ from \f$K_S^0\f$ true PID
     */
    NTuple::Item<int> m_TagKSPiMinusTrueID;
    /**
     * The tag \f$\pi^+\f$ true PID
     */
    NTuple::Item<int> m_TagPiPlusTrueID;
    /**
     * The tag \f$\pi^-\f$ true PID
     */
    NTuple::Item<int> m_TagPiMinusTrueID;
};

#endif
