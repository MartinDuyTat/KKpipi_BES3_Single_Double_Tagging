// Martin Duy Tat 12th February 2021
/**
 * FindKKTagInfo is a class for extracting all the variables of a KK tag
 */

#ifndef FINDKKTAGINFO
#define FINDKKTAGINFO

// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>
#include<string>

class FindKKTagInfo {
  public: 
    /**
     * Default constructor that initalizes all variables to zero
     */
    FindKKTagInfo();
    /**
     * Trivial destructor
     */
    ~FindKKTagInfo();
    /**
     * Function that calculates all the tag information and saves them
     * @param DTTool_iter Iterator pointing to tag candidate
     * @param DTTool DTagTool object with all the tag information
     * @param Returns true if successful
     */
    StatusCode CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool);
    /**
     * Enumeration to label daughter particles in the order K pi
     */
    enum DaughterParticle {KPLUS, KMINUS};
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
  private:
    /**
     * \f$K^+\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KPlusP;
    /**
     * \f$K^-\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KMinusP;
};

#endif
