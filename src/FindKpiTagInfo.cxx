// Martin Duy Tat 12th February 2021, based on code by Yu Zhang

// KKpipi
#include "KKpipi/FindKpiTagInfo.h"
// Gaudi
#include "GaudiKernel/SmartRefVector.h"
// Event information
#include "EvtRecEvent/EvtRecTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
// ROOT
#include "TMath.h"
// STL
#include<vector>
#include<string>
// Particle masses
#include "KKpipi/ParticleMasses.h"

FindKpiTagInfo::FindKpiTagInfo(): m_KCharge(0), m_PiCharge(0) {
}

FindKpiTagInfo::~FindKpiTagInfo() {
}

StatusCode FindKpiTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    if(DTTool.isKaon(*Track_iter)) {
      m_KP = MDCKalTrack->p4(MASS::K_MASS);
      m_KCharge = MDCKalTrack->charge();
    } else if(DTTool.isPion(*Track_iter)) {
      m_PiP = MDCKalTrack->p4(MASS::PI_MASS);
      m_PiCharge = MDCKalTrack->charge();
    }
  }
  return StatusCode::SUCCESS;
}

double FindKpiTagInfo::GetKP(int i) const {
  return m_KP[i];
}

double FindKpiTagInfo::GetPiP(int i) const {
  return m_PiP[i];
}

int FindKpiTagInfo::GetKCharge() const {
  return m_KCharge;
}

int FindKpiTagInfo::GetPiCharge() const {
  return m_PiCharge;
}
