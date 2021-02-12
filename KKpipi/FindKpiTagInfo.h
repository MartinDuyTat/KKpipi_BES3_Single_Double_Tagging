// Martin Duy Tat 12th February 2021, based on code by Yu Zhang
/**
 * FindKKpipiTag is a class for extracting all the variables of a Kpi tag
 */

#ifndef FINDKPITAG
#define FINDKPITAG

// Gaudi
#include "GaudiKernel/StatusCode.h"
// BOSS
#include "DTagTool/DTagTool.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// STL
#include<vector>
#include<string>

class FindKpiTagInfo {
  public: 
    /**
     * Default constructor that initalizes all variables to zero
     */
    FindKpiTagInfo();
    /**
     * Trivial destructor
     */
    ~FindKpiTagInfo();
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
    enum DaughterParticle {KAON, PION};
    /**
     * Get \f$K\f$ momentum component
     * @param i Component
     */
    double GetKP(int i);
    /**
     * Get \f$\pi\f$ momentum component
     * @param i Component
     */
    double GetPiP(int i);
    /**
     * Get \f$K\f$ charge
     */
    int GetKCharge();
    /**
     * Get \f$\pi\f$ charge
     */
    int GetPiCharge();
  private:
    /**
     * \f$K\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_KP;
    /**
     * \f$\pi\f$ four-momentum
     */
    CLHEP::HepLorentzVector m_PiP;
    /**
     * \f$K\f$ charge
     */
    int m_KCharge
    /**
     * \f$\pi\f$ charge
     */
    int m_PiCharge
};

#endif
