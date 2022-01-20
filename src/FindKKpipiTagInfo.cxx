// Martin Duy Tat 11th February 2021, based on code by Yu Zhang

// KKpipi
#include "KKpipi/FindKKpipiTagInfo.h"
#include "KKpipi/FindKS.h"
// Gaudi
#include "GaudiKernel/SmartRefVector.h"
// Event information
#include "EvtRecEvent/EvtRecTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
// ROOT
#include "TMath.h"
// STL
#include<vector>
#include<string>
// Particle masses
#include "KKpipi/ParticleMasses.h"
// Utilities
#include "KKpipi/KKpipiUtilities.h"

FindKKpipiTagInfo::FindKKpipiTagInfo(): m_DaughterTrackID(std::vector<int>(4)), m_KalmanFitSuccess(0), m_KalmanFitChi2(0.0), m_KSFitSuccess(0), m_DecayLengthVeeVertex(0.0), m_Chi2VeeVertex(0.0), m_KSMassVeeVertex(0.0), m_DecayLengthFit(0.0), m_DecayLengthErrorFit(0.0), m_Chi2Fit(0.0) {
}

FindKKpipiTagInfo::~FindKKpipiTagInfo() {
}

StatusCode FindKKpipiTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  std::vector<SmartRefVector<EvtRecTrack>::iterator> DaughterTrackIterators(4); // In the order K+ K- pi+ pi-
  std::vector<RecMdcKalTrack*> KalmanTracks(4); //In the order K+ K- pi+ pi-
  // Go through all tracks and find the daughter tracks
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    if(DTTool.isKaon(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	DaughterTrackIterators[KPLUS] = Track_iter;
	KalmanTracks[KPLUS] = MDCKalTrack;
	m_KPlusP = MDCKalTrack->p4(MASS::K_MASS);
	m_DaughterTrackID[KPLUS] = (*Track_iter)->trackId();
	KKpipiUtilities::GetIP(MDCKalTrack, m_KPlusIP_Vxy, m_KPlusIP_Vz);
      } else if (MDCKalTrack->charge() == -1) {
	DaughterTrackIterators[KMINUS] = Track_iter;
	KalmanTracks[KMINUS] = MDCKalTrack;
	m_KMinusP = MDCKalTrack->p4(MASS::K_MASS);
	m_DaughterTrackID[KMINUS] = (*Track_iter)->trackId();
	KKpipiUtilities::GetIP(MDCKalTrack, m_KMinusIP_Vxy, m_KMinusIP_Vz);
      }
    } else if(DTTool.isPion(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	DaughterTrackIterators[PIPLUS] = Track_iter;
	KalmanTracks[PIPLUS] = MDCKalTrack;
	m_PiPlusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID[PIPLUS] = (*Track_iter)->trackId();
      } else if(MDCKalTrack->charge() == -1) {
	DaughterTrackIterators[PIMINUS] = Track_iter;
	KalmanTracks[PIMINUS] = MDCKalTrack;
	m_PiMinusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID[PIMINUS] = (*Track_iter)->trackId();
      }
    }
  }
  // Prepare tracks for a Kalman kinematic fit
  WTrackParameter WTrackKplus(MASS::K_MASS, KalmanTracks[KPLUS]->getZHelix(), KalmanTracks[KPLUS]->getZError());
  WTrackParameter WTrackKminus(MASS::K_MASS, KalmanTracks[KMINUS]->getZHelix(), KalmanTracks[KMINUS]->getZError());
  WTrackParameter WTrackPIplus(MASS::PI_MASS, KalmanTracks[PIPLUS]->getZHelix(), KalmanTracks[PIPLUS]->getZError());
  WTrackParameter WTrackPIminus(MASS::PI_MASS, KalmanTracks[PIMINUS]->getZHelix(), KalmanTracks[PIMINUS]->getZError());
  KalmanKinematicFit *KalmanFit = KalmanKinematicFit::instance();
  KalmanFit->init();
  KalmanFit->AddTrack(KPLUS, WTrackKplus);
  KalmanFit->AddTrack(KMINUS, WTrackKminus);
  KalmanFit->AddTrack(PIPLUS, WTrackPIplus);
  KalmanFit->AddTrack(PIMINUS, WTrackPIminus);
  // Do a Kalman kinematic fit where the \f$D\f$ mass is constrained
  KalmanFit->AddResonance(0, MASS::D_MASS, KPLUS, KMINUS, PIPLUS, PIMINUS);
  m_KalmanFitSuccess = KalmanFit->Fit();
  if(m_KalmanFitSuccess) {
    m_KalmanFitChi2 = KalmanFit->chisq();
    m_KPlusPKalmanFit = KalmanFit->pfit(KPLUS);
    m_KMinusPKalmanFit = KalmanFit->pfit(KMINUS);
    m_PiPlusPKalmanFit = KalmanFit->pfit(PIPLUS);
    m_PiMinusPKalmanFit = KalmanFit->pfit(PIMINUS);
  }
  m_Mpipi = (m_PiPlusP + m_PiMinusP).m();
  m_KSFitSuccess = 0;
  // Check if the \f$\pi\pi\f$ pair is a \f$K_S\f$ in disguise
  if(m_Mpipi - MASS::KS_MASS < 0.050 && m_Mpipi - MASS::KS_MASS > -0.060) {
    FindKS findKS(false);
    std::vector<int> PionTrackIDs;
    PionTrackIDs.push_back((*DaughterTrackIterators[PIPLUS])->trackId());
    PionTrackIDs.push_back((*DaughterTrackIterators[PIMINUS])->trackId());
    StatusCode statuscode = findKS.findKS(DTTool_iter, DTTool, PionTrackIDs);
    m_DecayLengthVeeVertex = findKS.GetDecayLengthVeeVertex();
    m_Chi2VeeVertex = findKS.GetChi2VeeVertex();
    m_KSMassVeeVertex = findKS.GetKSMassVeeVertex();
    if(statuscode == StatusCode::SUCCESS) {
      m_KSFitSuccess = 1;
      m_DecayLengthFit = findKS.GetDecayLengthFit();
      m_DecayLengthErrorFit = findKS.GetDecayLengthErrorFit();
      m_Chi2Fit = findKS.GetChi2Fit();
    }
  }
  return StatusCode::SUCCESS;
}

std::vector<int> FindKKpipiTagInfo::GetDaughterTrackID() const {
  return m_DaughterTrackID;
}

double FindKKpipiTagInfo::GetKPlusP(int i) const {
  return m_KPlusP[i];
}

double FindKKpipiTagInfo::GetKMinusP(int i) const {
  return m_KMinusP[i];
}

double FindKKpipiTagInfo::GetPiPlusP(int i) const {
  return m_PiPlusP[i];
}

double FindKKpipiTagInfo::GetPiMinusP(int i) const {
  return m_PiMinusP[i];
}

int FindKKpipiTagInfo::GetKalmanFitSuccess() const {
  return m_KalmanFitSuccess;
}

double FindKKpipiTagInfo::GetKalmanFitChi2() const {
  return m_KalmanFitChi2;
}

double FindKKpipiTagInfo::GetKPlusPKalmanFit(int i) const {
  return m_KPlusPKalmanFit[i];
}

double FindKKpipiTagInfo::GetKMinusPKalmanFit(int i) const {
  return m_KMinusPKalmanFit[i];
}

double FindKKpipiTagInfo::GetPiPlusPKalmanFit(int i) const {
  return m_PiPlusPKalmanFit[i];
}

double FindKKpipiTagInfo::GetPiMinusPKalmanFit(int i) const {
  return m_PiMinusPKalmanFit[i];
}

double FindKKpipiTagInfo::GetMpipi() const {
  return m_Mpipi;
}

int FindKKpipiTagInfo::GetKSFitSuccess() const {
  return m_KSFitSuccess;
}

double FindKKpipiTagInfo::GetDecayLengthVeeVertex() const {
  return m_DecayLengthVeeVertex;
}

double FindKKpipiTagInfo::GetChi2VeeVertex() const {
  return m_Chi2VeeVertex;
}

double FindKKpipiTagInfo::GetKSMassVeeVertex() const {
  return m_KSMassVeeVertex;
}

double FindKKpipiTagInfo::GetDecayLengthFit() const {
  return m_DecayLengthFit;
}

double FindKKpipiTagInfo::GetDecayLengthErrorFit() const {
  return m_DecayLengthErrorFit;
}

double FindKKpipiTagInfo::GetChi2Fit() const {
  return m_Chi2Fit;
}

double FindKKpipiTagInfo::GetKPlusIP_Vxy() const {
  return m_KPlusIP_Vxy;
}

double FindKKpipiTagInfo::GetKPlusIP_Vz() const {
  return m_KPlusIP_Vz;
}

double FindKKpipiTagInfo::GetKMinusIP_Vxy() const {
  return m_KMinusIP_Vxy;
}

double FindKKpipiTagInfo::GetKMinusIP_Vz() const {
  return m_KMinusIP_Vz;
}
