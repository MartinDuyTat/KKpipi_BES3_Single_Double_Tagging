// Martin Duy Tat 25th March 2021, based on code by Yu Zhang
/**
 * FindKL is class for finding \f$K_L^0\f$ candidates in \f$K_L\pi^0\f$, \f$K_L\pi^0\pi^0\f$ and \f$K_L\omega\f$ tag
 * It reconstructs the pions on the other side of the signal and looks for missing energy and mass
 */

#ifndef FINDKL
#define FINDKL
// Gaudi
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/StatusCode.h"
// Boss
#include "DTagTool/DTagTool.h"
// Event information
#include "EvtRecEvent/EvtRecTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>

class FindKL {
  public: 
    /**
     * Default constructor, initializes everything to zero
     */
    FindKL();
    /**
     * Trivial destructor
     */
    ~FindKL();
    /**
     * Start looking for \f$K_L\f$ in the event
     * @param DTTool_iter DTagTool iterator pointing to the event with the tag
     * @param DTTool DTagTool object with all the tag information
     */
    StatusCode findKL(DTagToolIterator DTTool_iter, DTagTool DTTool);
    /**
     * Get flag that is 1 if a \f$\pi^\pi^-\f$ pair is found
     */
    int GetFoundPionPair() const;
    /**
     * Get the ith component of the \f$\pi^+\f$ four-momentum vector on the other side
     */
    double GetPiPlusP(int i) const;
    /**
     * Get the ith component of the \f$\pi^-\f$ four-momentum vector on the other side
     */
    double GetPiMinusP(int i) const;
    /**
     * Get the ith component of the jth high energy photon four-momentum from \f$\pi^0\f$
     */
    double GetPi0HighEPhotonP(int i, int j) const;
    /**
     * Get the ith component of the jth low energy photon four-momentum from \f$\pi^0\f$
     */
    double GetPi0LowEPhotonP(int i, int j) const;
    /**
     * Get the ith component of the jth high energy constrained photon four-momentum from \f$\pi^0\f$
     */
    double GetPi0HighEPhotonPConstrained(int i, int j) const;
    /**
     * Get the ith component of the jth low energy constrained photon four-momentum from \f$\pi^0\f$
     */
    double GetPi0LowEPhotonPConstrained(int i, int j) const;
    /**
     * Get the \f$\chi^2\f$ from the kinematic fit of the jth \f$\pi^0\f$
     */
    double GetPi0Chi2Fit(int j) const;
    /**
     * Get the jth high energy photon track ID from \f$\pi^0\f$
     */
    int GetPi0HighEPhotonTrackID(int j) const;
    /**
     * Get the jth low energy photon track ID from \f$\pi^0\f$
     */
    int GetPi0LowEPhotonTrackID(int j) const;
    /**
     * Get number of \f$\pi^0\f$ found
     */
    int GetNumberPi0() const;
    /**
     * Get the ith component of the jth high energy photon four-momentum from \f$\eta\f$
     */
    double GetEtaHighEPhotonP(int i, int j) const;
    /**
     * Get the ith component of the jth low energy photon four-momentum from \f$\eta\f$
     */
    double GetEtaLowEPhotonP(int i, int j) const;
    /**
     * Get the ith component of the jth high energy constrained photon four-momentum from \f$\eta\f$
     */
    double GetEtaHighEPhotonPConstrained(int i, int j) const;
    /**
     * Get the ith component of the jth low energy constrained photon four-momentum from \f$\eta\f$
     */
    double GetEtaLowEPhotonPConstrained(int i, int j) const;
    /**
     * Get the \f$\chi^2\f$ from the kinematic fit of the jth \f$\eta\f$
     */
    double GetEtaChi2Fit(int j) const;
    /**
     * Get the jth high energy photon track ID from \f$\eta\f$
     */
    int GetEtaHighEPhotonTrackID(int j) const;
    /**
     * Get the jth low energy photon track ID from \f$\eta\f$
     */
    int GetEtaLowEPhotonTrackID(int j) const;
    /**
     * Get number of \f$\eta\f$ found
     */
    int GetNumberEta() const;
    /**
     * Get ith component of jth single photon four momentum
     */
    double GetPhotonP(int i, int j) const;
    /**
     * Get jth single photon angular separation from the nearest track
     */
    double GetPhotonAngleSeparation(int j) const;
    /**
     * Get jth single photon polar angle separation from the nearest track
     */
    double GetPhotonThetaSeparation(int j) const;
    /**
     * Get jth single photon azimuthal angle separation from the nearest track
     */
    double GetPhotonPhiSeparation(int j) const;
    /**
     * Get jth single photon track ID
     */
    int GetPhotonTrackID(int j) const;
    /**
     * Get number of \f$\gamma\f$ found
     */
    int GetNumberGamma() const;
    /**
     * Get vector of track IDs of \f$\pi^+\pi^-\f$ pair
     */
    std::vector<int> GetDaughterTrackID() const;
  private:
    /**
     * Flag that is 1 if a \f$\pi^\pi^-\f$ pair is found
     */
    int m_FoundPionPair;
    /**
     * The \f$\pi^+\f$ four-momentum vector on the other side
     */
    CLHEP::HepLorentzVector m_PiPlusP;
    /**
     * The \f$\pi^-\f$ four-momentum vector on the other side
     */
    CLHEP::HepLorentzVector m_PiMinusP;
    /**
     * Vector of all high energy photon four-momentum from \f$\pi^0\f$
     */
    std::vector<CLHEP::HepLorentzVector> m_Pi0HighEPhotonP;
    /**
     * Vector of all low energy photon four-momentum from \f$\pi^0\f$
     */
    std::vector<CLHEP::HepLorentzVector> m_Pi0LowEPhotonP;
    /**
     * Vector of all high energy photon constrained four-momentum from \f$\pi^0\f$
     */
    std::vector<CLHEP::HepLorentzVector> m_Pi0HighEPhotonPConstrained;
    /**
     * Vector of all low energy photon constrained four-momentum from \f$\pi^0\f$
     */
    std::vector<CLHEP::HepLorentzVector> m_Pi0LowEPhotonPConstrained;
    /**
     * Vector of all \f$\chi^2\f$ from the kinematic fit of \f$\pi^0\f$
     */
    std::vector<double> m_Pi0Chi2Fit;
    /**
     * Vector of all high energy photon track ID from \f$\pi^0\f$
     */
    std::vector<int> m_Pi0HighEPhotonTrackID;
    /**
     * Vector of all low energy photon trackID from \f$\pi^0\f$
     */
    std::vector<int> m_Pi0LowEPhotonTrackID;
    /**
     * Number of \f$\pi^0\f$ found
     */
    int m_NumberPi0;
    /**
     * Vector of all high energy photon four-momentum from \f$\eta\f$
     */
    std::vector<CLHEP::HepLorentzVector> m_EtaHighEPhotonP;
    /**
     * Vector of all low energy photon four-momentum from \f$\eta\f$
     */
    std::vector<CLHEP::HepLorentzVector> m_EtaLowEPhotonP;
    /**
     * Vector of all high energy photon constrained four-momentum from \f$\eta\f$
     */
    std::vector<CLHEP::HepLorentzVector> m_EtaHighEPhotonPConstrained;
    /**
     * Vector of all low energy photon constrained four-momentum from \f$\eta\f$
     */
    std::vector<CLHEP::HepLorentzVector> m_EtaLowEPhotonPConstrained;
    /**
     * Vector of all \f$\chi^2\f$ from the kinematic fit of \f$\eta\f$
     */
    std::vector<double> m_EtaChi2Fit;
    /**
     * Vector of all high energy photon track ID from \f$\eta\f$
     */
    std::vector<int> m_EtaHighEPhotonTrackID;
    /**
     * Vector of all low energy photon trackID from \f$\eta\f$
     */
    std::vector<int> m_EtaLowEPhotonTrackID;
    /**
     * Number of \f$\eta\f$ found
     */
    int m_NumberEta;
    /**
     * Vector of all single photon four-momentum
     */
    std::vector<CLHEP::HepLorentzVector> m_PhotonP;
    /**
     * Vector of all single photon angular separation from the nearest track
     */
    std::vector<double> m_PhotonAngleSeparation;
    /**
     * Vector of all single photon polar angle separation from the nearest track
     */
    std::vector<double> m_PhotonThetaSeparation;
    /**
     * Vector of all single photon azimuthal angle separation from the nearest track
     */
    std::vector<double> m_PhotonPhiSeparation;
    /**
     * Vector of all single photon trackID
     */
    std::vector<int> m_PhotonTrackID;
    /**
     * Number of \f$\gamma\f$ found
     */
    int m_NumberGamma;
    /**
     * Vector of \f$\pi^+\pi^-\f$ daughter track IDs, in the order \f$\pi^+\f$ \f$\pi^-\f$
     */
    std::vector<int> m_DaughterTrackID;
};

#endif
