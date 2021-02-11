// Martin Duy Tat 28th January 2021, based on code by Yu Zhang

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

FindKKpipiTagInfo::FindKKpipiTagInfo(): m_KalmanFitSuccess(0), m_KalmanFitChi2(0.0), m_KSFitSuccess(0), m_DecayLengthVeeVertex(0.0), m_Chi2VeeVertex(0.0), m_KSMassVeeVertex(0.0, m_DecayLengthFit(0.0), m_DecayLengthErrorFit(0.0), m_Chi2Fit(0.0), m_KSMassFit(0.0) {
}

FindKKpipiTagInfo::~FindKKpipiTagInfo() {
}

StatusCode FindKKpipiTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  std::vector<SmartRefVector<EvtRecTrack>::iterator> DaughterTrackIterators(4); // In the order K+ K- pi+ pi-
  std::vector<RecMdcKalTrack*> KalmanTracks(4); //In the order K+ K- pi+ pi-
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    if(DTTool.isKaon(*Track_iter)) {

      if(MDCKalTrack->charge() == +1) {
	DaughterTrackIterators[KPLUS] = Track_iter;
	KalmanTracks[KPLUS] = MDCKalTrack;
	m_KPlusP = MDCKalTrack->p4(MASS::K_MASS);
      } else if (MDCKalTrack->charge() == -1) {
	DaughterTrackIterators[KMINUS] = Track_iter;
	KalmanTracks[KMINUS] = MDCKalTrack;
	m_KMinusP = MDCKalTrack->p4(MASS::K_MASS);
      }
    } else if(DTTool.isPion(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	DaughterTrackIterators[PIPLUS] = Track_iter;
	KalmanTracks[PIPLUS] = MDCKalTrack;
	m_PiPlusP = MDCKalTrack->p4(MASS::PI_MASS);
      } else if(MDCKalTrack->charge() == -1) {
	DaughterTrackIterators[PIMINUS] = Track_iter;
	KalmanTracks[PIMINUS] = MDCKalTrack;
	m_PiMinusP = MDCKalTrack->p4(MASS::PI_MASS);
      }
    }
  }
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
  KalmanFit->AddResonance(0, MASS::D_MASS, KPLUS, KMINUS, PIPLUS, PIMINUS);
  m_KalmanFitSuccess = KalmanFit->Fit();
  if(m_KalmanFitSuccess) {
    m_KalmanFitChi2 = KalmanFit->chisq();
    m_KPlusPKalmanFit = KalmanFit->pfit(KPLUS);
    m_KMinusPKalmanFit = KalmanFit->pfit(KMINUS);
    m_PiPlusPKalmanFit = KalmanFit->pfit(PIPLUS);
    M_PiMinusPKalmanFit = KalmanFit->pfit(PIMINUS);
  }
  double Mpipi = (m_PiPlusP + m_PiMinusP).m();
  m_KSFitSuccess = 0;
  if(TMath::Abs(Mpipi - MASS::KS_MASS) < 0.020) {
    FindKS findKS;
    std::vector<SmartRefVector<EvtRecTrack>::iterator> PionTracks_iter;
    PionTracks_iter.push_back(DaughterTrackIterators[PIPLUS]);
    PionTracks_iter.push_back(DaughterTrackIterators[PIMINUS]);
    StatusCode statuscode = findKS.findKS(DTTool_iter, PionTracks_iter);
    if(statuscode == StatusCode::SUCCESS) {
      m_KSFitSuccess = 1;
      m_DecayLengthVeeVertex = findKS.getDecayLengthVeeVertex();
      m_Chi2VeeVertex = findKS.getChi2VeeVertex();
      m_KSMassVeeVertex = findKS.getKSMassVeeVertex();
      m_DecayLengthFit = findKS.getDecayLengthFit();
      m_DecayLengthErrorFit = findKS.getDecayLengthErrorFit();
      m_Chi2Fit = findKS.getChi2Fit();
      m_KSMassFit = findKS.getKSMassFit();
    }
  }
  return StatusCode::SUCCESS;
}

double FindKKpipiTagInfo::GetKPlusP(int i) {
  return m_KPlusP[i];
}

double FindKKpipiTagInfo::GetKMinusP(int i) {
  return m_KMinusP[i];
}

double FindKKpipiTagInfo::GetPiPlusP(int i) {
  return m_PiPlusP[i];
}

double FindKKpipiTagInfo::GetPiMinusP(int i) {
  return m_PiMinusP[i];
}

int FindKKpipiTagInfo::GetKalmanFitSuccess() {
  return m_KalmanFitSuccess;
}

double FindKKpipiTagInfo::GetKalmanFitChi2() {
  return m_KalmanFitChi2;
}

double FindKKpipiTagInfo::GetKPlusPKalmanFit(int i) {
  return m_KPlusPKalmanFit[i];
}

double FindKKpipiTagInfo::GetKMinusPKalmanFit(int i) {
  return m_KMinusPKalmanFit[i];
}

double FindKKpipiTagInfo::GetPiPlusPKalmanFit(int i) {
  return m_PiPlusPKalmanFit[i];
}

double FindKKpipiTagInfo::GetPiMinusPKalmanFit(int i) {
  return m_PiMinusPKalmanFit[i];
}

int FindKKpipiTagInfo::GetKSFitSuccess() {
  return m_KSFitSuccess;
}

double FindKKpipiTagInfo::GetDecayLengthVeeVertex() {
  return m_DecayLengthVeeVertex;
}

double FindKKpipiTagInfo::GetChi2VeeVertex() {
  return m_Chi2VeeVertex;
}

double FindKKpipiTagInfo::GetKSMassVeeVertex() {
  return m_KSMassVeeVertex;
}

double FindKKpipiTagInfo::GetDecayLengthFit() {
  return m_DecayLengthFit;
}

double FindKKpipiTagInfo::GetDecayLengthErrorFit() {
  return m_DecayLengthErrorFit;
}

double FindKKpipiTagInfo::GetChi2Fit() {
  return m_Chi2Fit;
}

double FindKKpipiTagInfo::GetKSMassFit() {
  return m_KSMassFit;
}
