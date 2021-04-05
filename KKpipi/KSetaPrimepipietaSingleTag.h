// Martin Duy Tat 12th March 2021
/**
 * KSetaPrimepipietaSingleTag is a class for a BOSS algorithm
 * It runs over \f$D^0\bar{D^0}\f$ data and saves all events with a single \f$D\to K_S\eta'(\pi\pi\eta)f$ tag
 */

#ifndef KSETAPRIMEPIPIETASINGLETAG
#define KSETAPRIMEPIPIETASINGLETAG

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// STL
#include<string>

class KSetaPrimepipietaSingleTag: public Algorithm {
  public: 
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    KSetaPrimepipietaSingleTag(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~KSetaPrimepipietaSingleTag();
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
     * The \f$K_S^0\f$ mass, from fit
     */
    NTuple::Item<double> m_KSMassFit;
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
     * The \f$\pi^+\f$ daughter momentum along \f$x\f$ from the MDC track
     */
    NTuple::Item<double> m_KSPiPluspxFit;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$y\f$ after vertex fit
     */
    NTuple::Item<double> m_KSPiPluspyFit;
    /**
     * The \f$\pi^+\f$ daughter momentum along \f$z\f$ after vertex fit
     */
    NTuple::Item<double> m_KSPiPluspzFit;
    /**
     * The \f$\pi^+\f$ daughter energy after vertex fit
     */
    NTuple::Item<double> m_KSPiPlusenergyFit;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$x\f$ after vertex fit
     */
    NTuple::Item<double> m_KSPiMinuspxFit;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$y\f$ after vertex fit
     */
    NTuple::Item<double> m_KSPiMinuspyFit;
    /**
     * The \f$\pi^-\f$ daughter momentum along \f$z\f$ after vertex fit
     */
    NTuple::Item<double> m_KSPiMinuspzFit;
    /**
     * The \f$\pi^-\f$ daughter energy after vertex fit
     */
    NTuple::Item<double> m_KSPiMinusenergyFit;
    /**
     * High energy photon from \f$\eta\f$ unconstrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_HighEEtapx;
    /**
     * High energy photon from \f$\eta\f$ unconstrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_HighEEtapy;
    /**
     * High energy photon from \f$\eta\f$ unconstrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_HighEEtapz;
    /**
     * High energy photon from \f$\eta\f$ unconstrained energy
     */
    NTuple::Item<double> m_HighEEtaenergy;
    /**
     * Low energy photon from \f$\eta\f$ unconstrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_LowEEtapx;
    /**
     * Low energy photon from \f$\eta\f$ unconstrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_LowEEtapy;
    /**
     * Low energy photon from \f$\eta\f$ unconstrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_LowEEtapz;
    /**
     * Low energy photon from \f$\eta\f$ unconstrained energy
     */
    NTuple::Item<double> m_LowEEtaenergy;
    /**
     * The \f$\gamma\gamma\f$ invariant mass
     */
    NTuple::Item<double> m_Mgammagamma;
    /**
     * High energy photon from \f$\eta\f$ constrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_HighEEtaConstrainedpx;
    /**
     * High energy photon from \f$\eta\f$ constrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_HighEEtaConstrainedpy;
    /**
     * High energy photon from \f$\eta\f$ constrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_HighEEtaConstrainedpz;
    /**
     * High energy photon from \f$\eta\f$ constrained energy
     */
    NTuple::Item<double> m_HighEEtaConstrainedenergy;
    /**
     * Low energy photon from \f$\eta\f$ constrained momentum along \f$x\f$
     */
    NTuple::Item<double> m_LowEEtaConstrainedpx;
    /**
     * Low energy photon from \f$\eta\f$ constrained momentum along \f$y\f$
     */
    NTuple::Item<double> m_LowEEtaConstrainedpy;
    /**
     * Low energy photon from \f$\eta\f$ constrained momentum along \f$z\f$
     */
    NTuple::Item<double> m_LowEEtaConstrainedpz;
    /**
     * Low energy photon from \f$\eta\f$ constrained energy
     */
    NTuple::Item<double> m_LowEEtaConstrainedenergy;
    /**
     *  \f$\eta\f$ kinematic fit \f$\chi^2\f$
     */
    NTuple::Item<double> m_EtaChi2Fit;
    /**
     * \f$\pi\pi\eta\f$ invariant mass
     */
    NTuple::Item<double> m_Mpipieta;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of \f$\pi^+\pi^-\f$ tag tracks
     */
    NTuple::Item<int> m_pipiKSFitSuccess;
    /**
     * The \f$K_S\f$ decay length from \f$\pi^+\pi^-\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_pipiDecayLengthVeeVertex;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from \f$\pi^+\pi^-\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_pipiChi2VeeVertex;
    /**
     * The \f$K_S^0\f$ mass, from \f$\pi^+\pi^-\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_pipiKSMassVeeVertex;
    /**
     * The \f$K_S^0\f$ decay length, from \f$\pi^+\pi^-\f$, from fit
     */
    NTuple::Item<double> m_pipiDecayLengthFit;
    /**
     * The \f$K_S^0\f$ decay length error, from \f$\pi^+\pi^-\f$, from fit
     */
    NTuple::Item<double> m_pipiDecayLengthErrorFit;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from \f$\pi^+\pi^-\f$, from fit of primary vertex
     */
    NTuple::Item<double> m_pipiChi2Fit;
    /**
     * The \f$K_S^0\f$ mass, from \f$\pi^+\pi^-\f$, from fit
     */
    NTuple::Item<double> m_pipiKSMassFit;
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
     * \f$\pi^=\f$ energy
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
     * The \f$\pi^+\f$ from \f$\eta'\f$ true PID
     */
    NTuple::Item<int> m_EtaPPiPlusTrueID;
    /**
     * The \f$\pi^-\f$ from \f$\eta'\f$ true PID
     */
    NTuple::Item<int> m_EtaPPiMinusTrueID;
};

#endif
