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
     * @param VetoTracks Vector of tracks that should be skipped
     */
    FindhhTagInfo(std::string TagMode, const std::vector<int> &VetoTrackIDs = std::vector<int>());
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
     * Get \f$K^+\f$ momentum component
     * @param i Component
     */
    double GethPlusP(int i) const;
    /**
     * Get \f$K^-\f$ momentum component
     * @param i Component
     */
    double GethMinusP(int i) const;
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
    int m_PiPlusTrack_iter;
    /**
     * Track iterator of a \f$\pi^-\f$, if found
     */
    int m_PiMinusTrack_iter;
};

#endif
