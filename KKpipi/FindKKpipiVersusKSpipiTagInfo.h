// Martin Duy Tat 5th March 2021
/**
 * FindKKpipiVersusKSpipiTagInfo is a class for extracting all the variables of a KSpipi or KSKK tag
 */

#ifndef FINDKKPIPIVERSUSKSPIPITAGINFO
#define FINDKKPIPIVERSUSKSPIPITAGINFO

//KKpipi
#include "KKpipi/FindKS.h"
// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
#include "VertexFit/WTrackParameter.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>
#include<string>

struct WTrackParameters {
  WTrackParameter SignalKaon;
  WTrackParameter SignalPiPlus;
  WTrackParameter SignalPiMinus;
  WTrackParameter TagKSPiPlus;
  WTrackParameter TagKSPiMinus;
  WTrackParameter TagPiPlus;
  WTrackParameter TagPiMinus;
}

class FindKKpipiVersusKSpipiTagInfo {
  public: 
    /**
     * Default constructor that initalizes all variables to zero
     */
    FindKKpipiVersusKSpipiTagInfo(const std::string &TagMode = "KSpipi");
    /**
     * Trivial destructor
     */
    ~FindKKpipiVersusKSpipiTagInfo();
    /**
     * Function that calculates all the tag information and saves them
     * @param DTTool_iter Iterator pointing to tag candidate
     * @param DTTool DTagTool object with all the tag information
     * @param Returns true if successful
     */
    StatusCode CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool);
    /**
     * Function that performs the Kalman kinematic fit of the full decay tree
     * @param TrackParameters The track parameters of all particles for the Kalman fit
     * @param RecKCharge Charge of the reconstructed kaon
     */
    void DoKalmanFit(const WTrackParameters &TrackParameters, int RecKCharge);
    /**
     * Get the daughter track IDs, in the order K pi+ pi-
     */
    std::vector<int> GetDaughterTrackID_KKpipi() const;
    /**
     * Get the daughter track IDs, in the order (pi+ pi-)KS pi+ pi-
     */
    std::vector<int> GetDaughterTrackID_KSpipi() const;
    /**
     * Get flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tracks
     */
    int GetKSFitSuccess_KSpipi() const;
    /**
     * Get the \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double GetDecayLengthVeeVertex_KSpipi() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double GetChi2VeeVertex_KSpipi() const;
    /**
     * Get the \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double GetKSMassVeeVertex_KSpipi() const;
    /**
     * Get the \f$K_S^0\f$ decay length, from fit
     */
    double GetDecayLengthFit_KSpipi() const;
    /**
     * Get the \f$K_S^0\f$ decay length error, from fit
     */
    double GetDecayLengthErrorFit_KSpipi() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double GetChi2Fit_KSpipi() const;
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
    int GetKalmanFitSuccess_KSpipi() const;
    /**
     * Get Kalman fit \f$chi^2\f$
     */
    double GetKalmanFitChi2_KSpipi() const;
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
     * Get the \f$\pi^+\pi^-\f$ invariant mass
     */
    double GetMpipi_KKpipi() const;
    /**
     * Get flag of \f$K_S^0\f$ fit success of tracks
     */
    int GetKSFitSuccess_KKpipi() const;
    /**
     * Get the \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double GetDecayLengthVeeVertex_KKpipi() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double GetChi2VeeVertex_KKpipi() const;
    /**
     * Get the \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double GetKSMassVeeVertex_KKpipi() const;
    /**
     * Get the \f$K_S^0\f$ decay length, from fit
     */
    double GetDecayLengthFit_KKpipi() const;
    /**
     * Get the \f$K_S^0\f$ decay length error, from fit
     */
    double GetDecayLengthErrorFit_KKpipi() const;
    /**
     * Get the \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double GetChi2Fit_KKpipi() const;
    /**
     * Get the missing miss squared
     */
    double GetMissingMass2() const;
    /**
     * Get reconstructed kaon charge
     */
    int GetRecKCharge() const;
  private:
    /**
     * Daughter track IDs, in the order K+ K- pi+ pi-
     */
    std::vector<int> m_DaughterTrackID_KKpipi;
    /**
     * Daughter track IDs, in the order (pi+ pi-)KS h+ h-
     */
    std::vector<int> m_DaughterTrackID_KSpipi;
    /**
     * Flag equal to 1 for success and 0 for fail in the \f$K_S^0\f$ fit of tracks
     */
    int m_KSFitSuccess_KSpipi;
    /**
     * The \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double m_DecayLengthVeeVertex_KSpipi;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double m_Chi2VeeVertex_KSpipi;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double m_KSMassVeeVertex_KSpipi;
    /**
     * The \f$K_S^0\f$ decay length, from fit
     */
    double m_DecayLengthFit_KSpipi;
    /**
     * The \f$K_S^0\f$ decay length error, from fit
     */
    double m_DecayLengthErrorFit_KSpipi;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double m_Chi2Fit_KSpipi;
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
    int m_KSFitSuccess_KKpipi;
    /**
     * The \f$K_S\f$ decay length, from VeeVertexAlg
     */
    double m_DecayLengthVeeVertex_KKpipi;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from VeeVertexAlg
     */
    double m_Chi2VeeVertex_KKpipi;
    /**
     * The \f$K_S^0\f$ mass, from VeeVertexAlg
     */
    double m_KSMassVeeVertex_KKpipi;
    /**
     * The \f$K_S^0\f$ decay length, from fit
     */
    double m_DecayLengthFit_KKpipi;
    /**
     * The \f$K_S^0\f$ decay length error, from fit
     */
    double m_DecayLengthErrorFit_KKpipi;
    /**
     * The \f$K_S^0\f$ \f$\chi^2\f$, from fit of primary vertex
     */
    double m_Chi2Fit_KKpipi;
    /**
     * The missing mass squared
     */
    double m_MissingMass2;
    /**
     * The reconstructed kaon charge
     */
    int m_RecKCharge;
};

#endif
