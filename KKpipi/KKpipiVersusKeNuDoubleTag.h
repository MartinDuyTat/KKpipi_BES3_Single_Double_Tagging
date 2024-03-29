// Martin Duy Tat 30th March 2021
/**
 * KKpipiVersusKeNuDoubleTag is a class for a BOSS algorithm
 * It runs over \f$D^0\bar{D^0}\f$ data and saves all events with a double \f$D\to K^+K^-\pi^+\pi^-\f$ vs \f$D^0\to Ke\nu\f$ tag
 * It also runs a fit for the decay \f$K_S^0\to\pi^+\pi^-\f$ by refitting the primary and secondary vertex in the class FindKS, from this the flight significance is used to eliminate peaking background
 */

#ifndef KKPIPIVERSUSKENUDOUBLETAG
#define KKPIPIVERSUSKENUDOUBLETAG

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// STL
#include<string>

class KKpipiVersusKeNuDoubleTag: public Algorithm {
  public: 
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    KKpipiVersusKeNuDoubleTag(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~KKpipiVersusKeNuDoubleTag();
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
    StatusCode FillTuple(DTagToolIterator DTTool_Signal_iter, DTagTool &DTTool);
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
     * Signal \f$\pi^+\pi^-\f$ invariant mass
     */
    NTuple::Item<double> m_SignalMpipi;
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
     * Tag \f$K\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagKpx;
    /**
     * Tag \f$K\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagKpy;
    /**
     * Tag \f$K\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagKpz;
    /**
     * Tag \f$K\f$ energy
     */
    NTuple::Item<double> m_TagKenergy;
    /**
     * Tag \f$K\f$ charge
     */
    NTuple::Item<int> m_TagKCharge;
    /**
     * Tag \f$e\f$ momentum, including FSR, along \f$x\f$
     */
    NTuple::Item<double> m_TagElectronpx;
    /**
     * Tag \f$e\f$ momentum, including FSR, along \f$y\f$
     */
    NTuple::Item<double> m_TagElectronpy;
    /**
     * Tag \f$e\f$ momentum, including FSR, along \f$z\f$
     */
    NTuple::Item<double> m_TagElectronpz;
    /**
     * Tag \f$e\f$, including FSR, energy
     */
    NTuple::Item<double> m_TagElectronenergy;
    /**
     * Tag \f$e\f$ charge
     */
    NTuple::Item<int> m_TagElectronCharge;
    /**
     * Tag FSR momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagFSRpx;
    /**
     * Tag FSR momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagFSRpy;
    /**
     * Tag FSR momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagFSRpz;
    /**
     * Tag FSR energy
     */
    NTuple::Item<double> m_TagFSRenergy;
    /**
     * Tag missing momentum along \f$x\f$
     */
    NTuple::Item<double> m_TagMissingpx;
    /**
     * Tag missing momentum along \f$y\f$
     */
    NTuple::Item<double> m_TagMissingpy;
    /**
     * Tag missing momentum along \f$z\f$
     */
    NTuple::Item<double> m_TagMissingpz;
    /**
     * Tag missing energy
     */
    NTuple::Item<double> m_TagMissingenergy;
    /**
     * \f$U_\text{miss} = E_\text{miss} - |p_\text{miss}|\f$
     */
    NTuple::Item<double> m_TagUMiss;
    /**
     * Number of good photons that are not FSR
     */
    NTuple::Item<int> m_TagNumberGamma;
    /**
     * Array of energies of good photon showers that are not FSR
     */
    NTuple::Array<double> m_TagExtraShowerEnergy;
    /**
     * Equal to 1 if the tag daughter tracks are from the same \f$D\f$ meson
     */
    NTuple::Item<int> m_TagIsSameDMother;
    /**
     * Equal to 1 if the tag daughter tracks are assigned a PID matching that of the MC truth
     */
    NTuple::Item<int> m_TagPIDTrue;
    /**
     * The tag \f$K\f$ true PID
     */
    NTuple::Item<int> m_TagKTrueID;
    /**
     * The tag \f$e\f$ true PID
     */
    NTuple::Item<int> m_TagElectronTrueID;
};

#endif
