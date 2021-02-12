// Martin Duy Tat 12th February 2021

// KKpipi
#include "KKpipi/FindKKTagInfo.h"
// Gaudi
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/StatusCode.h"
// Event information
#include "EvtRecEvent/EvtRecTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
// STL
#include<vector>
// Particle masses
#include "KKpipi/ParticleMasses.h"

FindKKTagInfo::FindKKTagInfo() {
}

FindKKTagInfo::~FindKKTagInfo() {
}

StatusCode FindKKTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    if(DTTool.isKaon(*Track_iter) && MDCKalTrack->charge() > 0) {
      m_KPlusP = MDCKalTrack->p4(MASS::K_MASS);
    } else if(DTTool.isKaon(*Track_iter) && MDCKalTrack->charge() < 0) {
      m_KMinusP = MDCKalTrack->p4(MASS::K_MASS);
    }
  }
  return StatusCode::SUCCESS;
}

double FindKKTagInfo::GetKPlusP(int i) const {
  return m_KPlusP[i];
}

double FindKKTagInfo::GetKMinusP(int i) const {
  return m_KMinusP[i];
}
