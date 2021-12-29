// Martin Duy Tat 5th March 2021
/**
 * FindKShhTagInfo is a class for extracting all the variables of a KSpipi or KSKK tag
 */

#ifndef FINDKSPIPITAGINFO
#define FINDKSPIPITAGINFO

//KKpipi
#include "KKpipi/FindKS.h"
// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>
#include<string>

class FindKShhTagInfo {
  public: 
    /**
     * Default constructor that initalizes all variables to zero
     */
    FindKShhTagInfo(const std::string &TagMode = "KSpipi");
    /**
     * Trivial destructor
     */
    ~FindKShhTagInfo();
    /**
     * Function that calculates all the tag information and saves them
     * @param DTTool_iter Iterator pointing to tag candidate
     * @param DTTool DTagTool object with all the tag information
     * @param Returns true if successful
     */
    StatusCode CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool);
    /**
     * Get the daughter track IDs, in the order (pi+ pi-)KS pi+ pi-
     */
    std::vector<int> GetDaughterTrackID() const;
    /**
     * Enumeration to label daughter particles in the order h+ h- KS0
     */
    enum DaughterParticle {HPLUS, HMINUS, KSHORT};
    /**
     * Get flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tracks
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
     * Get the \f$\pi^+\f$ daughter momentum component from the MDC track
     */
    double GetKSPiPlusP(int i) const;
    /**
     * Get the \f$\pi^-\f$ daughter momentum component from the MDC track
     */
    double GetKSPiMinusP(int i) const;
    /**
     * Get \f$K_S^0\f$ momentum component after vertex fit
     * @param i Component
     */
    double GetKShortP(int i) const;
    /**
     * Get \f$h^+\f$ momentum component
     * @param i Component
     */
    double GethPlusP(int i) const;
    /**
     * Get \f$h^-\f$ momentum component
     * @param i Component
     */
    double GethMinusP(int i) const;
    /**
     * Get flag of Kalman fit success
     */
    int GetKalmanFitSuccess() const;
    /**
     * Get Kalman fit \f$chi^2\f$
     */
    double GetKalmanFitChi2() const;
    /**
     * Get \f$K_S^0\f$ momentum component from Kalman fit
     * @param i Component
     */
    double GetKShortPKalmanFit(int i) const;
    /**
     * Get \f$h^+\f$ momentum component from Kalman fit
     * @param i Component
     */
    double GethPlusPKalmanFit(int i) const;
    /**
     * Get \f$h^-\f$ momentum component from Kalman fit
     * @param i Component
     */
    double GethMinusPKalmanFit(int i) const;
    /**
     * Get flag of \f$K_S^0\f$ from \f$\pi^+\pi^-\f$ fit of tracks success
     */
    int GetpipiKSFitSuccess() const;
    /**
     * Get the \f$K_S\f$ decay length from \f$\pi^+\pi^-\f$, from VeeVertexAlg
     */
    double GetpipiDecayLengthVeeVertex() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$ from \f$\pi^+\pi^-\f$, from VeeVertexAlg
     */
    double GetpipiChi2VeeVertex() const;
    /**
     * Get the \f$K_S^0\f$ mass from \f$\pi^+\pi^-\f$, from VeeVertexAlg
     */
    double GetpipiKSMassVeeVertex() const;
    /**
     * Get the \f$K_S^0\f$ decay length from \f$\pi^+\pi^-\f$, from fit
     */
    double GetpipiDecayLengthFit() const;
    /**
     * Get the \f$K_S^0\f$ decay length error from \f$\pi^+\pi^-\f$, from fit
     */
    double GetpipiDecayLengthErrorFit() const;
    /**
     * Get the \f$K_S^0\f$ from \f$\pi^+\pi^-\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double GetpipiChi2Fit() const;
  private:
    /**
     * "KSpipi" or "KSKK"
     */
    std::string m_TagMode;
    /**
     * Daughter track IDs, in the order (pi+ pi-)KS h+ h-
     */
    std::vector<int> m_DaughterTrackID;
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
     * The \f$\pi^+\f$ daughter momentum from the MDC track
     */
    CLHEP::HepLorentzVector m_KSPiPlusP;
    /**
     * The \f$\pi^-\f$ daughter momentum from the MDC track
     */
    CLHEP::HepLorentzVector m_KSPiMinusP;
    /**
     * \f$K_S^0\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KShortP;
    /**
     * \f$h^+\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_hPlusP;
    /**
     * \f$h^-\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_hMinusP;
    /**
     * Flag equal to 1 for success and 0 for fail in the Kalman fit of tracks
     */
    int m_KalmanFitSuccess;
    /**
     * \f$\chi^2\f$ of Kalman fit
     */
    double m_KalmanFitChi2;
    /**
     * Kalman fitted \f$K_S^0\f$ four- momentum
     */
    CLHEP::HepLorentzVector m_KShortPKalmanFit;
    /**
     * Kalman fitted \f$h^+\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_hPlusPKalmanFit;
    /**
     * Kalman fitted \f$h^-\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_hMinusPKalmanFit;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ from \f$\pi^+\pi^-\f$ fit of tracks
     */
    int m_pipiKSFitSuccess;
    /**
     * The \f$K_S\f$ from \f$\pi^+\pi^-\f$ decay length, from VeeVertexAlg
     */
    double m_pipiDecayLengthVeeVertex;
    /**
     * The \f$K_S^0\f$ from \f$\pi^+\pi^-\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double m_pipiChi2VeeVertex;
    /**
     * The \f$K_S^0\f$ from \f$\pi^+\pi^-\f$ mass, from VeeVertexAlg
     */
    double m_pipiKSMassVeeVertex;
    /**
     * The \f$K_S^0\f$ from \f$\pi^+\pi^-\f$ decay length, from fit
     */
    double m_pipiDecayLengthFit;
    /**
     * The \f$K_S^0\f$ from \f$\pi^+\pi^-\f$ decay length error, from fit
     */
    double m_pipiDecayLengthErrorFit;
    /**
     * The \f$K_S^0\f$ from \f$\pi^+\pi^-\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double m_pipiChi2Fit;
};

#endif
