// Martin Duy Tat 5th February 2021, based on code by Yu Zhang
/**
 * FindKS is class for finding \f$K_S^0\f$ candidates and returning the flight distance, the flight distance error, the fit \f$\chi^2\f$ and the \f$K_S^0\f$ mass from VeeVertexAlg
 * In addition, the secondary vertex is fitted again using the \f$\pi^+\pi^-\f$ tracks and a primary vertex is fitted using the vertex parameters from the secondary vertex fit
 */

#ifndef FINDKS
#define FINDKS
// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
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
    StatusCode findKS(DTagToolIterator &DTTool_iter, const std::vector<int> &PiTrackIndex);
    /** 
     * Get decay length from VeeVertexAlg
     */
    double getDecayLengthVeeVertex() const;
    /** 
     * Get \f$\chi^2\f$ from VeeVertexAlg
     */
    double getChi2VeeVertex() const;
    /** 
     * Get \f$K_S^0\f$ mass from VeeVertexAlg
     */
    double getKSMassVeeVertex() const;
    /** 
     * Get decay length from fit
     */
    double getDecayLengthFit() const;
    /** 
     * Get decay length error from fit
     */
    double getDecayLengthErrorFit() const;
    /** 
     * Get \f$\chi^2\f$ from fit of primary vertex
     */
    double getChi2PrimaryVertexFit() const;
    /** 
     * Get \f$\chi^2\f$ from fit of secondary vertex
     */
    double getChi2SecondaryVertexFit() const;
    /** 
     * Get \f$K_S^0\f$ from fit
     */
    double getKSMassFit() const;
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
    double m_Chi2PrimaryVertexFit;
    /**
     * The \f$\chi^2\f$, from fit of secondary vertex
     */
    double m_Chi2SecondaryVertexFit;
    /**
     * The \f$K_S^0\f$ mass, from fit
     */
    double m_KSMassFit;
};

#endif
