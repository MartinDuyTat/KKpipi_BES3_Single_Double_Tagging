// Martin Duy Tat 26th march 2021
/**
 * KKpipiVersusKLDoubleTags is a class for a BOSS algorithm
 * It runs over \f$D^0\bar{D^0}\f$ data and saves all events with a double \f$D\to K^+K^-\pi^+\pi^-\f$ vs \f$K_L\pi^0\f$, \f$K_L\pi^0\pi^0\f$ and \f$K_L\omega\f$ tags
 * It also runs a fit for the decay \f$K_S^0\to\pi^+\pi^-\f$ by refitting the primary and secondary vertex in the class FindKS, from this the flight significance is used to eliminate peaking background
 */

#ifndef KKPIPIVERSUSKLDOUBLETAGS
#define KKPIPIVERSUSKLDOUBLETAGS

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// STL
#include<string>

class KKpipiVersusKLDoubleTags: public Algorithm {
  public: 
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    KKpipiVersusKLDoubleTags(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~KKpipiVersusKLDoubleTags();
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
    /**
     * Helper function that fills out the missing energy and and missing mass, after everything else has been filled
     */
    void FillMissingMassEnergy();
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
     * Tag flag equal to 1 if a \f$\pi^+\pi^-\f$ pair is found
     */
    NTuple::Item<int> m_TagFoundPionPair;
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
     * Tag \f$\pi^+\f$ energy
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
     * Number of \f$\pi^0\f$ found
     */
    NTuple::Item<int> m_TagNumberPi0;
    /**
     * Tag array of high energy photon from \f$\pi^0\f$ momentum along \f$x\f$
     */
    NTuple::Array<double> m_TagPi0HighEPhotonpx;
    /**
     * Tag array of high energy photon from \f$\pi^0\f$ momentum along \f$y\f$
     */
    NTuple::Array<double> m_TagPi0HighEPhotonpy;
    /**
     * Tag array of high energy photon from \f$\pi^0\f$ momentum along \f\z\f$
     */
    NTuple::Array<double> m_TagPi0HighEPhotonpz;
    /**
     * Tag array of high energy photon from \f$\pi^0\f$ energy
     */
    NTuple::Array<double> m_TagPi0HighEPhotonenergy;
    /**
     * Tag array of low energy photon from \f$\pi^0\f$ momentum along \f$x\f$
     */
    NTuple::Array<double> m_TagPi0LowEPhotonpx;
    /**
     * Tag array of low energy photon from \f$\pi^0\f$ momentum along \f$y\f$
     */
    NTuple::Array<double> m_TagPi0LowEPhotonpy;
    /**
     * Tag array of low energy photon from \f$\pi^0\f$ momentum along \f\z\f$
     */
    NTuple::Array<double> m_TagPi0LowEPhotonpz;
    /**
     * Tag array of low energy photon from \f$\pi^0\f$ energy
     */
    NTuple::Array<double> m_TagPi0LowEPhotonenergy;
    /**
     * Tag array of \f$\gamma\gamma\f$ invariant masses from \f$\pi^0\f$
     */
    NTuple::Array<double> m_TagPi0Mgammagamma;
    /**
     * Tag array of high energy photon from \f$\pi^0\f$ constrained momentum along \f$x\f$
     */
    NTuple::Array<double> m_TagPi0HighEPhotonpxConstrained;
    /**
     * Tag array of high energy photon from \f$\pi^0\f$ constrained momentum along \f$y\f$
     */
    NTuple::Array<double> m_TagPi0HighEPhotonpyConstrained;
    /**
     * Tag array of high energy photon from \f$\pi^0\f$ constrained momentum along \f\z\f$
     */
    NTuple::Array<double> m_TagPi0HighEPhotonpzConstrained;
    /**
     * Tag array of high energy photon from \f$\pi^0\f$ constrained energy
     */
    NTuple::Array<double> m_TagPi0HighEPhotonenergyConstrained;
    /**
     * Tag array of low energy photon from \f$\pi^0\f$ constrained momentum along \f$x\f$
     */
    NTuple::Array<double> m_TagPi0LowEPhotonpxConstrained;
    /**
     * Tag array of low energy photon from \f$\pi^0\f$ constrained momentum along \f$y\f$
     */
    NTuple::Array<double> m_TagPi0LowEPhotonpyConstrained;
    /**
     * Tag array of low energy photon from \f$\pi^0\f$ constrained momentum along \f\z\f$
     */
    NTuple::Array<double> m_TagPi0LowEPhotonpzConstrained;
    /**
     * Tag array of low energy photon from \f$\pi^0\f$ constrained energy
     */
    NTuple::Array<double> m_TagPi0LowEPhotonenergyConstrained;
    /**
     * Array of \f$\chi^2\f$ from kinematic fit of photons from \f$\pi^0\f$
     */
    NTuple::Array<double> m_TagPi0Chi2Fit;
    /**
     * Tag array of high energy photon from \f$\pi^0\f$ track ID
     */
    NTuple::Array<double> m_TagPi0HighEPhotonTrackID;
    /**
     * Tag array of low energy photon from \f$\pi^0\f$ track ID
     */
    NTuple::Array<double> m_TagPi0LowEPhotonTrackID;
    /**
     * Number of \f$\eta\f$ found
     */
    NTuple::Item<int> m_TagNumberEta;
    /**
     * Tag array of high energy photon from \f$\eta\f$ momentum along \f$x\f$
     */
    NTuple::Array<double> m_TagEtaHighEPhotonpx;
    /**
     * Tag array of high energy photon from \f$\eta\f$ momentum along \f$y\f$
     */
    NTuple::Array<double> m_TagEtaHighEPhotonpy;
    /**
     * Tag array of high energy photon from \f$\eta\f$ momentum along \f\z\f$
     */
    NTuple::Array<double> m_TagEtaHighEPhotonpz;
    /**
     * Tag array of high energy photon from \f$\eta\f$ energy
     */
    NTuple::Array<double> m_TagEtaHighEPhotonenergy;
    /**
     * Tag array of low energy photon from \f$\eta\f$ momentum along \f$x\f$
     */
    NTuple::Array<double> m_TagEtaLowEPhotonpx;
    /**
     * Tag array of low energy photon from \f$\eta\f$ momentum along \f$y\f$
     */
    NTuple::Array<double> m_TagEtaLowEPhotonpy;
    /**
     * Tag array of low energy photon from \f$\eta\f$ momentum along \f\z\f$
     */
    NTuple::Array<double> m_TagEtaLowEPhotonpz;
    /**
     * Tag array of low energy photon from \f$\eta\f$ energy
     */
    NTuple::Array<double> m_TagEtaLowEPhotonenergy;
    /**
     * Tag array of \f$\gamma\gamma\f$ invariant masses from \f$\eta\f$
     */
    NTuple::Array<double> m_TagEtaMgammagamma;
    /**
     * Tag array of high energy photon from \f$\eta\f$ constrained momentum along \f$x\f$
     */
    NTuple::Array<double> m_TagEtaHighEPhotonpxConstrained;
    /**
     * Tag array of high energy photon from \f$\eta\f$ constrained momentum along \f$y\f$
     */
    NTuple::Array<double> m_TagEtaHighEPhotonpyConstrained;
    /**
     * Tag array of high energy photon from \f$\eta\f$ constrained momentum along \f\z\f$
     */
    NTuple::Array<double> m_TagEtaHighEPhotonpzConstrained;
    /**
     * Tag array of high energy photon from \f$\eta\f$ constrained energy
     */
    NTuple::Array<double> m_TagEtaHighEPhotonenergyConstrained;
    /**
     * Tag array of low energy photon from \f$\eta\f$ constrained momentum along \f$x\f$
     */
    NTuple::Array<double> m_TagEtaLowEPhotonpxConstrained;
    /**
     * Tag array of low energy photon from \f$\eta\f$ constrained momentum along \f$y\f$
     */
    NTuple::Array<double> m_TagEtaLowEPhotonpyConstrained;
    /**
     * Tag array of low energy photon from \f$\eta\f$ constrained momentum along \f\z\f$
     */
    NTuple::Array<double> m_TagEtaLowEPhotonpzConstrained;
    /**
     * Tag array of low energy photon from \f$\eta\f$ constrained energy
     */
    NTuple::Array<double> m_TagEtaLowEPhotonenergyConstrained;
    /**
     * Array of \f$\chi^2\f$ from kinematic fit of photons from \f$\eta\f$
     */
    NTuple::Array<double> m_TagEtaChi2Fit;
    /**
     * Tag array of high energy photon from \f$\eta\f$ track ID
     */
    NTuple::Array<double> m_TagEtaHighEPhotonTrackID;
    /**
     * Tag array of low energy photon from \f$\eta\f$ track ID
     */
    NTuple::Array<double> m_TagEtaLowEPhotonTrackID;
    /**
     * Number of single photons found
     */
    NTuple::Item<int> m_TagNumberGamma;
    /**
     * Array of single photon energies
     */
    NTuple::Array<double> m_TagPhotonEnergy;
    /**
     * Array of single photon momenta along \f$x\f$
     */
    NTuple::Array<double> m_TagPhotonPx;
    /**
     * Array of single photon momenta along \f$y\f$
     */
    NTuple::Array<double> m_TagPhotonPy;
    /**
     * Array of single photon momenta along \f$z\f$
     */
    NTuple::Array<double> m_TagPhotonPz;
    /**
     * Array of single photon angular separation to nearest track
     */
    NTuple::Array<double> m_TagPhotonAngleSeparation;
    /**
     * Array of single photon polar angle separation to nearest track
     */
    NTuple::Array<double> m_TagPhotonThetaSeparation;
    /**
     * Array of single photon azimuthal angle separation to nearest track
     */
    NTuple::Array<double> m_TagPhotonPhiSeparation;
    /**
     * Tag array of single photon track ID
     */
    NTuple::Array<double> m_TagPhotonTrackID;
    /**
     * Missing energy on the tag side
     */
    NTuple::Item<double> m_TagMissingEnergy;
    /**
     * Missing mass squared on the tag side
     */
    NTuple::Item<double> m_TagMissingMass2;
    /**
     * Equal to 1 if the tag daughter tracks are from the same \f$D\f$ meson
     */
    NTuple::Item<int> m_TagIsSameDMother;
    /**
     * Equal to 1 if the tag daughter tracks are assigned a PID matching that of the MC truth
     */
    NTuple::Item<int> m_TagPIDTrue;
    /**
     * The tag \f$\pi^+\f$ true PID
     */
    NTuple::Item<int> m_TagPiPlusTrueID;
    /**
     * The tag \f$\pi^-\f$ true PID
     */
    NTuple::Item<int> m_TagPiMinusTrueID;
    /**
     * The tag \f$\pi^+\f$ \f$D\f$ origin
     */
    NTuple::Item<int> m_TagPiPlusDOrigin;
    /**
     * The tag \f$\pi^-\f$ \f$D\f$ origin
     */
    NTuple::Item<int> m_TagPiMinusDOrigin;
    /**
     * The tag high energy photon from \f$\pi^0\f$ true PID
     */
    NTuple::Array<int> m_TagPi0HighEPhotonTrueID;
    /**
     * The tag low energy photon from \f$\pi^0\f$ true PID
     */
    NTuple::Array<int> m_TagPi0LowEPhotonTrueID;
    /**
     * The tag high energy photon from \f$\eta\f$ true PID
     */
    NTuple::Array<int> m_TagEtaHighEPhotonTrueID;
    /**
     * The tag low energy photon from \f$\eta\f$ true PID
     */
    NTuple::Array<int> m_TagEtaLowEPhotonTrueID;
    /**
     * The tag single photon true PID
     */
    NTuple::Array<int> m_TagPhotonTrueID;
    /**
     * The tag high energy photon from \f$\pi^0\f$ true mother PID
     */
    NTuple::Array<int> m_TagPi0HighEPhotonMotherTrueID;
    /**
     * The tag low energy photon from \f$\pi^0\f$ true mother PID
     */
    NTuple::Array<int> m_TagPi0LowEPhotonMotherTrueID;
    /**
     * The tag high energy photon from \f$\eta\f$ true mother PID
     */
    NTuple::Array<int> m_TagEtaHighEPhotonMotherTrueID;
    /**
     * The tag low energy photon from \f$\eta\f$ true mother PID
     */
    NTuple::Array<int> m_TagEtaLowEPhotonMotherTrueID;
    /**
     * The tag high energy photon from \f$\pi^0\f$ \f$D\f$ meson origin
     */
    NTuple::Array<int> m_TagPi0HighEPhotonDOrigin;
    /**
     * The tag low energy photon from \f$\pi^0\f$ \f$D\f$ meson origin
     */
    NTuple::Array<int> m_TagPi0LowEPhotonDOrigin;
    /**
     * The tag high energy photon from \f$\eta\f$ \f$D\f$ meson origin
     */
    NTuple::Array<int> m_TagEtaHighEPhotonDOrigin;
    /**
     * The tag low energy photon from \f$\eta\f$ \f$D\f$ meson origin
     */
    NTuple::Array<int> m_TagEtaLowEPhotonDOrigin;
    /**
     * The tag single photon \f$D\f$ meson origin
     */
    NTuple::Array<int> m_TagPhotonDOrigin;
};

#endif
