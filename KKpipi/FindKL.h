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
  private:
    /**
     * Flag that is 1 if a \f$\pi^\pi^-\f$ pair is found
     */
    int FoundPionPair;
    /**
     * The \f$\pi^+\f$ four-momentum vector on the other side
     */
    CLHEP::HepLorentzVector m_PiPlusP;
    /**
     * The \f$\pi^-\f$ four-momentum vector on the other side
     */
    CLHEP::HepLorentzVector m_PiMinusP;
};

#endif
