// Martin Duy Tat 12th February 2021
/**
 * FindhhTagInfo is a class for extracting all the variables of a \f$KK\f$ tag or a \f$\pi\pi\f$ tag
 */

#ifndef FINDHHTAGINFO
#define FINDHHTAGINFO

// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>
#include<string>

class FindhhTagInfo {
  public: 
    /**
     * Default constructor that initalizes all variables to zero, and specifies if it's a \f$KK\f$ or a \f$\pi\pi\f$ tag
     * @param TagMode String that either says "KK" or "pipi"
     * @param VetoTracks Vector of tracks that should be skipped (optional)
     * @param PIDrequirement Set this to false to turn off PID requirements for pions (optional)
     */
    FindhhTagInfo(std::string TagMode, const std::vector<int> &VetoTrackIDs = std::vector<int>(), bool PIDrequirement = true);
    /**
     * Trivial destructor
     */
    ~FindhhTagInfo();
    /**
     * Function that calculates all the tag information and saves them
     * @param DTTool_iter Iterator pointing to tag candidate
     * @param DTTool DTagTool object with all the tag information
     * @param Returns true if successful
     */
    StatusCode CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool);
    /**
     * Get the daughter track IDs
     */
    std::vector<int> GetDaughterTrackID() const;
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
     * Get the \f$hh\f$ invariant mass
     */
    double GetMhh() const;
    /**
     * Get track iterator of a \f$\pi^+\f$, if found
     */
    int GetPiPlusTrackID() const;
    /**
     * Get track iterator of a \f$\pi^-\f$, if found
     */
    int GetPiMinusTrackID() const;
  private:
    /**
     * Daughter track IDs
     */
    std::vector<int> m_DaughterTrackID;
    /**
     * Tag mode, either "KK" or "pipi"
     */
    std::string m_TagMode;
    /**
     * Vector of tracks that should be skipped
     * Use this when there are other \f$K_S\f$ or \f$\phi\f$ in the tag that are reconstructed separately
     */
    std::vector<int> m_VetoTrackIDs;
    /**
     * \f$h^+\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_hPlusP;
    /**
     * \f$h^-\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_hMinusP;
    /**
     * Track iterator of a \f$\pi^+\f$, if found
     */
    int m_PiPlusTrackID;
    /**
     * Track iterator of a \f$\pi^-\f$, if found
     */
    int m_PiMinusTrackID;
    /**
     * Flag that indicates whether or not pions are subject to PID requirements
     * Set to false for pions from eta' usually
     */
    bool m_PIDrequirement;
};

#endif
