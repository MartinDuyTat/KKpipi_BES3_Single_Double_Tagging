// Martin Duy Tat 5th March 2021

// KKpipi
#include "KKpipi/FindKSpipiTagInfo.h"
#include "KKpipi/FindKS.h"
#include "KKpipi/FindhhTagInfo.h"
#include "KKpipi/ParticleMasses.h"
// Gaudi
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/StatusCode.h"
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

FindKSpipiTagInfo::FindKSpipiTagInfo(): m_DaughterTrackID(std::vector<int>(4)), m_DecayLengthVeeVertex(0.0), m_Chi2VeeVertex(0.0), m_KSMassVeeVertex(0.0), m_DecayLengthFit(0.0), m_DecayLengthErrorFit(0.0), m_Chi2Fit(0.0), m_KSMassFit(0.0), m_KalmanFitSuccess(0), m_KalmanFitChi2(0.0), m_pipiKSFitSuccess(0), m_pipiDecayLengthVeeVertex(0.0), m_pipiChi2VeeVertex(0.0), m_pipiKSMassVeeVertex(0.0), m_pipiDecayLengthFit(0.0), m_pipiDecayLengthErrorFit(0.0), m_pipiChi2Fit(0.0), m_pipiKSMassFit(0.0) {
}

FindKSpipiTagInfo::~FindKSpipiTagInfo() {
}

StatusCode FindKSpipiTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  // Find the actual KS candidate
  FindKS findKS(true);
  StatusCode status = findKS.findKS(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  // Fill out the information about the KS candidate
  m_DecayLengthVeeVertex = findKS.GetDecayLengthVeeVertex();
  m_Chi2VeeVertex = findKS.GetChi2VeeVertex();
  m_KSMassVeeVertex = findKS.GetKSMassVeeVertex();
  m_DecayLengthFit = findKS.GetDecayLengthFit();
  m_DecayLengthErrorFit = findKS.GetDecayLengthErrorFit();
  m_Chi2Fit = findKS.GetChi2Fit();
  m_KSMassFit = findKS.GetKSMassFit();
  m_KSPiPlusP = CLHEP::HepLorentzVector(findKS.GetKSPiPlusP(0), findKS.GetKSPiPlusP(1), findKS.GetKSPiPlusP(2), findKS.GetKSPiPlusP(3));
  m_KSPiMinusP = CLHEP::HepLorentzVector(findKS.GetKSPiMinusP(0), findKS.GetKSPiMinusP(1), findKS.GetKSPiMinusP(2), findKS.GetKSPiMinusP(3));
  m_KSPiPlusPFit = CLHEP::HepLorentzVector(findKS.GetKSPiPlusPFit(0), findKS.GetKSPiPlusPFit(1), findKS.GetKSPiPlusPFit(2), findKS.GetKSPiPlusPFit(3));
  m_KSPiMinusPFit = CLHEP::HepLorentzVector(findKS.GetKSPiMinusPFit(0), findKS.GetKSPiMinusPFit(1), findKS.GetKSPiMinusPFit(2), findKS.GetKSPiMinusPFit(3));
  m_KShortP = findKS.GetKShortPFit();
  // Get the track ID of the KS daughters
  std::vector<int> KSDaughterTrackIDs = findKS.GetDaughterTrackIDs();
  m_DaughterTrackID[0] = KSDaughterTrackIDs[0];
  m_DaughterTrackID[1] = KSDaughterTrackIDs[1];
  // Get all tracks
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  std::vector<SmartRefVector<EvtRecTrack>::iterator> DaughterTrackIterators(2); // In the order pi+ pi-
  std::vector<RecMdcKalTrack*> KalmanTracks(2); //In the order pi+ pi-
  std::vector<RecMdcKalTrack*> KSDaughterKalmanTracks;
  // Loop over all tracks
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    // If track is from KS daughters, skip
    if(std::find(KSDaughterTrackIDs.begin(), KSDaughterTrackIDs.end(), (*Track_iter)->trackId()) != KSDaughterTrackIDs.end()) {
      KSDaughterKalmanTracks.push_back(MDCKalTrack);
      continue;
    }
    // Fill out track information
    if(DTTool.isPion(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	DaughterTrackIterators[PIPLUS] = Track_iter;
	KalmanTracks[PIPLUS] = MDCKalTrack;
	m_PiPlusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID[2] = (*Track_iter)->trackId();
      } else if(MDCKalTrack->charge() == -1) {
	DaughterTrackIterators[PIMINUS] = Track_iter;
	KalmanTracks[PIMINUS] = MDCKalTrack;
	m_PiMinusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID[3] = (*Track_iter)->trackId();
      }
    }
  }
  // Do a Kalman kinematic fit of the KS daughter tracks, and constrain the KS mass to its PDG value
  WTrackParameter WTrackKSPIplus(MASS::PI_MASS, KSDaughterKalmanTracks[PIPLUS]->getZHelix(), KSDaughterKalmanTracks[PIPLUS]->getZError());
  WTrackParameter WTrackKSPIminus(MASS::PI_MASS, KSDaughterKalmanTracks[PIMINUS]->getZHelix(), KSDaughterKalmanTracks[PIMINUS]->getZError());
  KalmanKinematicFit *KSKalmanFit = KalmanKinematicFit::instance();
  KSKalmanFit->init();
  KSKalmanFit->AddTrack(PIPLUS, WTrackKSPIplus);
  KSKalmanFit->AddTrack(PIMINUS, WTrackKSPIminus);
  KSKalmanFit->AddResonance(0, MASS::KS_MASS, PIPLUS, PIMINUS);
  bool KSKalmanFitSuccess = KSKalmanFit->Fit();
  // Do a Kalman kinematic fit of the D daughter tracks, and constrain the KS and D masses to their PDG values
  if(KSKalmanFitSuccess) {
    KSKalmanFit->BuildVirtualParticle(0);
    WTrackParameter WTrackKS = KSKalmanFit->wVirtualTrack(0);
    WTrackParameter WTrackPIplus(MASS::PI_MASS, KalmanTracks[PIPLUS]->getZHelix(), KalmanTracks[PIPLUS]->getZError());
    WTrackParameter WTrackPIminus(MASS::PI_MASS, KalmanTracks[PIMINUS]->getZHelix(), KalmanTracks[PIMINUS]->getZError());
    KalmanKinematicFit *KalmanFit = KalmanKinematicFit::instance();
    KalmanFit->init();
    KalmanFit->AddTrack(PIPLUS, WTrackPIplus);
    KalmanFit->AddTrack(PIMINUS, WTrackPIminus);
    KalmanFit->AddTrack(KSHORT, WTrackKS);
    KalmanFit->AddResonance(0, MASS::D_MASS, PIPLUS, PIMINUS, KSHORT);
    m_KalmanFitSuccess = KalmanFit->Fit();
    if(m_KalmanFitSuccess) {
      m_KalmanFitChi2 = KalmanFit->chisq();
      m_KShortPKalmanFit = KalmanFit->pfit(KSHORT);
      m_PiPlusPKalmanFit = KalmanFit->pfit(PIPLUS);
      m_PiMinusPKalmanFit = KalmanFit->pfit(PIMINUS);
    }
  }
  // Do a vertex fit of the pipi in the tag, to make sure they're not KS in disguise
  double Mpipi = (m_PiPlusP + m_PiMinusP).m();
  m_pipiKSFitSuccess = 0;
  if(TMath::Abs(Mpipi - MASS::KS_MASS) < 0.020) {
    FindKS pipifindKS(false);
    std::vector<int> PionTrackIDs;
    PionTrackIDs.push_back((*DaughterTrackIterators[PIPLUS])->trackId());
    PionTrackIDs.push_back((*DaughterTrackIterators[PIMINUS])->trackId());
    StatusCode statuscode = pipifindKS.findKS(DTTool_iter, DTTool, PionTrackIDs);
    if(statuscode == StatusCode::SUCCESS) {
      m_pipiKSFitSuccess = 1;
      m_pipiDecayLengthVeeVertex = pipifindKS.GetDecayLengthVeeVertex();
      m_pipiChi2VeeVertex = pipifindKS.GetChi2VeeVertex();
      m_pipiKSMassVeeVertex = pipifindKS.GetKSMassVeeVertex();
      m_pipiDecayLengthFit = pipifindKS.GetDecayLengthFit();
      m_pipiDecayLengthErrorFit = pipifindKS.GetDecayLengthErrorFit();
      m_pipiChi2Fit = pipifindKS.GetChi2Fit();
      m_pipiKSMassFit = pipifindKS.GetKSMassFit();
    }
  }
  return StatusCode::SUCCESS;
}

std::vector<int> FindKSpipiTagInfo::GetDaughterTrackID() const {
  return m_DaughterTrackID;
}

double FindKSpipiTagInfo::GetDecayLengthVeeVertex() const {
  return m_DecayLengthVeeVertex;
}

double FindKSpipiTagInfo::GetChi2VeeVertex() const {
  return m_Chi2VeeVertex;
}

double FindKSpipiTagInfo::GetKSMassVeeVertex() const {
  return m_KSMassVeeVertex;
}

double FindKSpipiTagInfo::GetDecayLengthFit() const {
  return m_DecayLengthFit;
}

double FindKSpipiTagInfo::GetDecayLengthErrorFit() const {
  return m_DecayLengthErrorFit;
}

double FindKSpipiTagInfo::GetChi2Fit() const {
  return m_Chi2Fit;
}

double FindKSpipiTagInfo::GetKSMassFit() const {
  return m_KSMassFit;
}

double FindKSpipiTagInfo::GetKSPiPlusP(int i) const {
  return m_KSPiPlusP[i];
}

double FindKSpipiTagInfo::GetKSPiMinusP(int i) const {
  return m_KSPiMinusP[i];
}

double FindKSpipiTagInfo::GetKSPiPlusPFit(int i) const {
  return m_KSPiPlusPFit[i];
}

double FindKSpipiTagInfo::GetKSPiMinusPFit(int i) const {
  return m_KSPiMinusPFit[i];
}

double FindKSpipiTagInfo::GetKShortP(int i) const {
  return m_KShortP[i];
}

double FindKSpipiTagInfo::GetPiPlusP(int i) const {
  return m_PiPlusP[i];
}

double FindKSpipiTagInfo::GetPiMinusP(int i) const {
  return m_PiMinusP[i];
}

int FindKSpipiTagInfo::GetKalmanFitSuccess() const {
  return m_KalmanFitSuccess;
}

double FindKSpipiTagInfo::GetKalmanFitChi2() const {
  return m_KalmanFitChi2;
}

double FindKSpipiTagInfo::GetKShortPKalmanFit(int i) const {
  return m_KShortPKalmanFit[i];
}

double FindKSpipiTagInfo::GetPiPlusPKalmanFit(int i) const {
  return m_PiPlusPKalmanFit[i];
}

double FindKSpipiTagInfo::GetPiMinusPKalmanFit(int i) const {
  return m_PiMinusPKalmanFit[i];
}

int FindKSpipiTagInfo::GetpipiKSFitSuccess() const {
  return m_pipiKSFitSuccess;
}

double FindKSpipiTagInfo::GetpipiDecayLengthVeeVertex() const {
  return m_pipiDecayLengthVeeVertex;
}

double FindKSpipiTagInfo::GetpipiChi2VeeVertex() const {
  return m_pipiChi2VeeVertex;
}

double FindKSpipiTagInfo::GetpipiKSMassVeeVertex() const {
  return m_pipiKSMassVeeVertex;
}

double FindKSpipiTagInfo::GetpipiDecayLengthFit() const {
  return m_pipiDecayLengthFit;
}

double FindKSpipiTagInfo::GetpipiDecayLengthErrorFit() const {
  return m_pipiDecayLengthErrorFit;
}

double FindKSpipiTagInfo::GetpipiChi2Fit() const {
  return m_pipiChi2Fit;
}

double FindKSpipiTagInfo::GetpipiKSMassFit() const {
  return m_pipiKSMassFit;
}
