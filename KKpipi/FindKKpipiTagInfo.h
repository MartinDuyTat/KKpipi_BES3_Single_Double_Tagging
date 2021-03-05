// Martin Duy Tat 11th February 2021, based on code by Yu Zhang
/**
 * FindKKpipiTag is a class for extracting all the variables of a KKpipi tag
 */

#ifndef FINDKKPIPITAG
#define FINDKKPIPITAG

// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>
#include<string>

class FindKKpipiTagInfo {
  public: 
    /**
     * Default constructor that initalizes all variables to zero
     */
    FindKKpipiTagInfo();
    /**
     * Trivial destructor
     */
    ~FindKKpipiTagInfo();
    /**
     * Function that calculates all the tag information and saves them
     * @param DTTool_iter Iterator pointing to tag candidate
     * @param DTTool DTagTool object with all the tag information
     * @param Returns true if successful
     */
    StatusCode CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool);
    /**
     * Enumeration to label daughter particles in the order K+ K- pi+ pi-
     */
    enum DaughterParticle {KPLUS, KMINUS, PIPLUS, PIMINUS};
    /**
     * Get \f$K^+\f$ momentum component
     * @param i Component
     */
    double GetKPlusP(int i) const;
    /**
     * Get \f$K^-\f$ momentum component
     * @param i Component
     */
    double GetKMinusP(int i) const;
    /**
     * Get \f$\pi^+\f$ momentum component
     * @param i Component
     */
    double GetPiPlusP(int i) const;
    /**
     * Get \f$\pi^-\f$ momentum component
     * @param i Component
     */
    double GetPiMinusP(int i) const;
    /**
     * Get flag of Kalman fit success
     */
    int GetKalmanFitSuccess() const;
    /**
     * Get Kalman fit \f$chi^2\f$
     */
    double GetKalmanFitChi2() const;
    /**
     * Get \f$K^+\f$ momentum component from Kalman fit
     * @param i Component
     */
    double GetKPlusPKalmanFit(int i) const;
    /**
     * Get \f$K^-\f$ momentum component from Kalman fit
     * @param i Component
     */
    double GetKMinusPKalmanFit(int i) const;
    /**
     * Get \f$\pi^+\f$ momentum component from Kalman fit
     * @param i Component
     */
    double GetPiPlusPKalmanFit(int i) const;
    /**
     * Get \f$\pi^-\f$ momentum component from Kalman fit
     * @param i Component
     */
    double GetPiMinusPKalmanFit(int i) const;
    /**
     * Get flag of \f$K_S^0\f$ fit success of tracks
     */
    int GetKSFitSuccess() const;
    /**
     * Get the \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double GetDecayLengthVeeVertex() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double GetChi2VeeVertex() const;
    /**
     * Get the \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double GetKSMassVeeVertex() const;
    /**
     * Get the \f$K_S^0\f$ decay length, from fit
     */
    double GetDecayLengthFit() const;
    /**
     * Get the \f$K_S^0\f$ decay length error, from fit
     */
    double GetDecayLengthErrorFit() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double GetChi2Fit() const;
    /**
     * Get the \f$K_S^0\f$ mass, from fit
     */
    double GetKSMassFit() const;
  private:
    /**
     * \f$K^+\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KPlusP;
    /**
     * \f$K^-\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KMinusP;
    /**
     * \f$\pi^+\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_PiPlusP;
    /**
     * \f$\pi^-\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_PiMinusP;
    /**
     * Flag equal to 1 for success and 0 for fail in the Kalman fit of tracks
     */
    int m_KalmanFitSuccess;
    /**
     * \f$\chi^2\f$ of Kalman fit
     */
    double m_KalmanFitChi2;
    /**
     * Kalman fitted \f$K^+\f$ four- momentum
     */
    CLHEP::HepLorentzVector m_KPlusPKalmanFit;
    /**
     * Kalman fitted \f$K^-\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KMinusPKalmanFit;
    /**
     * Kalman fitted \f$\pi^+\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_PiPlusPKalmanFit;
    /**
     * Kalman fitted \f$\pi^-\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_PiMinusPKalmanFit;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tracks
     */
    int m_KSFitSuccess;
    /**
     * The \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double m_DecayLengthVeeVertex;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double m_Chi2VeeVertex;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double m_KSMassVeeVertex;
    /**
     * The \f$K_S^0\f$ decay length, from fit
     */
    double m_DecayLengthFit;
    /**
     * The \f$K_S^0\f$ decay length error, from fit
     */
    double m_DecayLengthErrorFit;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double m_Chi2Fit;
    /**
     * The \f$K_S^0\f$ mass, from fit
     */
    double m_KSMassFit;
};

#endif
