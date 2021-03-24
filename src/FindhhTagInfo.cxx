// Martin Duy Tat 12th February 2021

// KKpipi
#include "KKpipi/FindhhTagInfo.h"
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
#include<algorithm>
// Particle masses
#include "KKpipi/ParticleMasses.h"

FindhhTagInfo::FindhhTagInfo(std::string TagMode, const std::vector<int> &VetoTrackIDs): m_DaughterTrackID(std::vector<int>(2)), m_TagMode(TagMode), m_VetoTrackIDs(VetoTrackIDs) {
}

FindhhTagInfo::~FindhhTagInfo() {
}

StatusCode FindhhTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    // If track is on the veto list, skip
    if(m_VetoTrackIDs.size() != 0 && std::find(m_VetoTrackIDs.begin(), m_VetoTrackIDs.end(), (*Track_iter)->trackId()) != m_VetoTrackIDs.end()) {
      continue;
    }
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    if(m_TagMode == "KK") {
      if(DTTool.isKaon(*Track_iter) && MDCKalTrack->charge() > 0) {
	m_hPlusP = MDCKalTrack->p4(MASS::K_MASS);
	m_DaughterTrackID[0] = (*Track_iter)->trackId();
      } else if(DTTool.isKaon(*Track_iter) && MDCKalTrack->charge() < 0) {
	m_hMinusP = MDCKalTrack->p4(MASS::K_MASS);
	m_DaughterTrackID[1] = (*Track_iter)->trackId();
      }
    } else if(m_TagMode == "pipi") {
      if(DTTool.isPion(*Track_iter) && MDCKalTrack->charge() > 0) {
	m_hPlusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_PiPlusTrackID = (*Track_iter)->trackId();
	m_DaughterTrackID[0] = (*Track_iter)->trackId();
      } else if(DTTool.isPion(*Track_iter) && MDCKalTrack->charge() < 0) {
	m_hMinusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_PiMinusTrackID = (*Track_iter)->trackId();
	m_DaughterTrackID[1] = (*Track_iter)->trackId();
      }
    }
  }
  return StatusCode::SUCCESS;
}

std::vector<int> FindhhTagInfo::GetDaughterTrackID() const {
  return m_DaughterTrackID;
}

double FindhhTagInfo::GethPlusP(int i) const {
  return m_hPlusP[i];
}

double FindhhTagInfo::GethMinusP(int i) const {
  return m_hMinusP[i];
}

int FindhhTagInfo::GetPiPlusTrackID() const {
  return m_PiPlusTrackID;
}

int FindhhTagInfo::GetPiMinusTrackID() const {
  return m_PiMinusTrackID;
}
