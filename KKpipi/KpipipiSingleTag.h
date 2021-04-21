// Martin Duy Tat 25th March
/**
 * KpipipiSingleTag is a class for a BOSS algorithm
 * It runs over \f$D^0\bar{D^0}\f$ data and saves all events with a single \f$D\to K\pi\pi\pi\f$ tag
 * It also runs a fit for the decay \f$K_S^0\to\pi^+\pi^-\f$ by refitting the primary and secondary vertex in the class FindKS, from this the flight significance is used to eliminate peaking background
 */

#ifndef KPIPIPISINGLETAG
#define KPIPIPISINGLETAG

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// STL
#include<string>

class KpipipiSingleTag: public Algorithm {
  public: 
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    KpipipiSingleTag(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~KpipipiSingleTag();
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
     * \f$K\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_Kpx;
    /**
     * \f$K\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_Kpy;
    /**
     * \f$K\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_Kpz;
    /**
     * \f$K\f$ energy
     */
    NTuple::Item<double> m_Kenergy;
    /**
     * \f$\pi\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_Pi1px;
    /**
     * \f$\pi\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_Pi1py;
    /**
     * \f$\pi\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_Pi1pz;
    /**
     * \f$\pi\f$ energy
     */
    NTuple::Item<double> m_Pi1energy;
    /**
     * \f$\pi\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_Pi2px;
    /**
     * \f$\pi\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_Pi2py;
    /**
     * \f$\pi\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_Pi2pz;
    /**
     * \f$\pi\f$ energy
     */
    NTuple::Item<double> m_Pi2energy;
    /**
     * \f$\pi\f$ momentum along \f$x\f$
     */
    NTuple::Item<double> m_Pi3px;
    /**
     * \f$\pi\f$ momentum along \f$y\f$
     */
    NTuple::Item<double> m_Pi3py;
    /**
     * \f$\pi\f$ momentum along \f$z\f$
     */
    NTuple::Item<double> m_Pi3pz;
    /**
     * \f$\pi\f$ energy
     */
    NTuple::Item<double> m_Pi3energy;
    /**
     * \f$K\f$ charge
     */
    NTuple::Item<int> m_KCharge;
    /**
     * \f$\pi\f$ charge
     */
    NTuple::Item<int> m_Pi1Charge;
    /**
     * \f$\pi\f$ charge
     */
    NTuple::Item<int> m_Pi2Charge;
    /**
     * \f$\pi\f$ charge
     */
    NTuple::Item<int> m_Pi3Charge;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tracks
     */
    NTuple::Item<int> m_12KSFitSuccess;
    /**
     * The \f$K_S\f$ decay length, from VeeVertexAlg
     */
    NTuple::Item<double> m_12DecayLengthVeeVertex;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_12Chi2VeeVertex;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    NTuple::Item<double> m_12KSMassVeeVertex;
    /**
     * The \f$K_S^0\f$ decay length, from fit
     */
    NTuple::Item<double> m_12DecayLengthFit;
    /**
     * The \f$K_S^0\f$ decay length error, from fit
     */
    NTuple::Item<double> m_12DecayLengthErrorFit;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    NTuple::Item<double> m_12Chi2Fit;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tracks
     */
    NTuple::Item<int> m_13KSFitSuccess;
    /**
     * The \f$K_S\f$ decay length, from VeeVertexAlg
     */
    NTuple::Item<double> m_13DecayLengthVeeVertex;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    NTuple::Item<double> m_13Chi2VeeVertex;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    NTuple::Item<double> m_13KSMassVeeVertex;
    /**
     * The \f$K_S^0\f$ decay length, from fit
     */
    NTuple::Item<double> m_13DecayLengthFit;
    /**
     * The \f$K_S^0\f$ decay length error, from fit
     */
    NTuple::Item<double> m_13DecayLengthErrorFit;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    NTuple::Item<double> m_13Chi2Fit;
    /**
     * Equal to 1 if the daughter tracks are from the same \f$D\f$ meson
     */
    NTuple::Item<int> m_IsSameDMother;
    /**
     * Equal to 1 if the daughter tracks are assigned a PID matching that of the MC truth
     */
    NTuple::Item<int> m_PIDTrue;
    /**
     * The \f$K\f$ true PID
     */
    NTuple::Item<int> m_KTrueID;
    /**
     * The \f$\pi\f$ true PID
     */
    NTuple::Item<int> m_Pi1TrueID;
    /**
     * The \f$\pi\f$ true PID
     */
    NTuple::Item<int> m_Pi2TrueID;
    /**
     * The \f$\pi\f$ true PID
     */
    NTuple::Item<int> m_Pi3TrueID;
};

#endif
