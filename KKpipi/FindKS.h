// Martin Duy Tat 5th February 2021, based on code by Yu Zhang
/**
 * FindKS is class for finding \f$K_S^0\f$ candidates and returning the flight distance, the flight distance error, the fit \f$\chi^2\f$ and the \f$K_S^0\f$ mass from VeeVertexAlg, plus the daughter four-momenta
 * In addition, the secondary vertex is fitted again using the \f$\pi^+\pi^-\f$ tracks and a primary vertex is fitted using the vertex parameters from the secondary vertex fit
 */

#ifndef FINDKS
#define FINDKS
// Gaudi
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/StatusCode.h"
// Event information
#include "EvtRecEvent/EvtRecTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>

class FindKS {
  public: 
    /**
     * Default constructor, initializes everything to zero
     */
    FindKS();
    /**
     * Trivial destructor
     */
    ~FindKS();
    /**
     * Start looking for \f$K_S\f$ in the event
     * @param DTTool_iter DTagTool iterator pointing to the event with the tag
     * @param PiTrackIndex List of length 2 with track indices to the two pions in the event
     */
    StatusCode findKS(const std::vector<SmartRefVector<EvtRecTrack>::iterator> &PiTrack_iter = std::vector<SmartRefVector<EvtRecTrack>::iterator>());
    /** 
     * Get decay length from VeeVertexAlg
     */
    double GetDecayLengthVeeVertex() const;
    /** 
     * Get \f$\chi^2\f$ from VeeVertexAlg
     */
    double GetChi2VeeVertex() const;
    /** 
     * Get \f$K_S^0\f$ mass from VeeVertexAlg
     */
    double GetKSMassVeeVertex() const;
    /** 
     * Get decay length from fit
     */
    double GetDecayLengthFit() const;
    /** 
     * Get decay length error from fit
     */
    double GetDecayLengthErrorFit() const;
    /** 
     * Get \f$\chi^2\f$ from fit of primary vertex
     */
    double GetChi2Fit() const;
    /** 
     * Get \f$K_S^0\f$ from fit
     */
    double GetKSMassFit() const;
    /**
     * Get the \f$\pi^+\f$ daughter four-momentum from the MDC track
     * @param i Momentum component
     */
    double GetKSPiPlusP(int i) const;
    /**
     * Get the \f$\pi^-\f$ daughter four-momentum from the MDC track
     * @param i Momentum component
     */
    double GetKSPiMinusP(int i) const;
    /**
     * Get the \f$\pi^+\f$ daughter four-momentum after vertex fit
     * @param i Momentum component
     */
    double GetKSPiPlusPFit(int i) const;
    /**
     * Get the \f$\pi^-\f$ daughter four-momentum after vertex fit
     * @param i Momentum component
     */
    double GetKSPiMinusPFit(int i) const;
  private:
    /**
     * The decay length, from VeeVertexAlg
     */
    double m_DecayLengthVeeVertex;
    /**
     * The \f$\chi^2\f$, from VeeVertexAlg
     */
    double m_Chi2VeeVertex;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double m_KSMassVeeVertex;
    /**
     * The decay length, from fit
     */
    double m_DecayLengthFit;
    /**
     * The decay length error, from fit
     */
    double m_DecayLengthErrorFit;
    /**
     * The \f$\chi^2\f$, from fit of primary vertex
     */
    double m_Chi2Fit;
    /**
     * The \f$K_S^0\f$ mass, from fit
     */
    double m_KSMassFit;
    /**
     * The \f$\pi^+\f$ daughter four-momentum from the MDC track
     */
    CLHEP::HepLorentzVector m_KSPiPlusP;
    /**
     * The \f$\pi^-\f$ daughter four-momentum from the MDC track
     */
    CLHEP::HepLorentzVector m_KSPiMinusP;
    /**
     * The \f$\pi^+\f$ daughter four-momentum after vertex fit
     */
    CLHEP::HepLorentzVector m_KSPiPlusPFit;
    /**
     * The \f$\pi^-\f$ daughter four-momentum after vertex fit
     */
    CLHEP::HepLorentzVector m_KSPiMinusPFit;
};

#endif
