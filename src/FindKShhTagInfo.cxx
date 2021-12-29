// Martin Duy Tat 5th March 2021

// KKpipi
#include "KKpipi/FindKShhTagInfo.h"
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

FindKShhTagInfo::FindKShhTagInfo(const std::string &TagMode): m_TagMode(TagMode), m_DaughterTrackID(std::vector<int>(4)), m_KSFitSuccess(0), m_DecayLengthVeeVertex(0.0), m_Chi2VeeVertex(0.0), m_KSMassVeeVertex(0.0), m_DecayLengthFit(0.0), m_DecayLengthErrorFit(0.0), m_Chi2Fit(0.0), m_KalmanFitSuccess(0), m_KalmanFitChi2(0.0), m_pipiKSFitSuccess(0), m_pipiDecayLengthVeeVertex(0.0), m_pipiChi2VeeVertex(0.0), m_pipiKSMassVeeVertex(0.0), m_pipiDecayLengthFit(0.0), m_pipiDecayLengthErrorFit(0.0), m_pipiChi2Fit(0.0) {
}

FindKShhTagInfo::~FindKShhTagInfo() {
}

StatusCode FindKShhTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  // Find the actual KS candidate
  FindKS findKS(true);
  StatusCode status = findKS.findKS(DTTool_iter, DTTool);
  if(status == StatusCode::SUCCESS) {
    m_KSFitSuccess = 1;
    m_DecayLengthFit = findKS.GetDecayLengthFit();
    m_DecayLengthErrorFit = findKS.GetDecayLengthErrorFit();
    m_Chi2Fit = findKS.GetChi2Fit();
  } else {
    m_KSFitSuccess = 0;
  }
  // Fill out the information about the KS candidate
  m_DecayLengthVeeVertex = findKS.GetDecayLengthVeeVertex();
  m_Chi2VeeVertex = findKS.GetChi2VeeVertex();
  m_KSMassVeeVertex = findKS.GetKSMassVeeVertex();
  m_KSPiPlusP = CLHEP::HepLorentzVector(findKS.GetKSPiPlusP(0), findKS.GetKSPiPlusP(1), findKS.GetKSPiPlusP(2), findKS.GetKSPiPlusP(3));
  m_KSPiMinusP = CLHEP::HepLorentzVector(findKS.GetKSPiMinusP(0), findKS.GetKSPiMinusP(1), findKS.GetKSPiMinusP(2), findKS.GetKSPiMinusP(3));
  m_KShortP = findKS.GetKShortPFit();
  // Get the track ID of the KS daughters
  std::vector<int> KSDaughterTrackIDs = findKS.GetDaughterTrackIDs();
  m_DaughterTrackID[0] = KSDaughterTrackIDs[0];
  m_DaughterTrackID[1] = KSDaughterTrackIDs[1];
  // Get all tracks
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  std::vector<SmartRefVector<EvtRecTrack>::iterator> DaughterTrackIterators(2); // In the order pi+ pi-
  std::vector<RecMdcKalTrack*> KalmanTracks(2); //In the order h+ h-
  std::vector<RecMdcKalTrack*> KSDaughterKalmanTracks;
  // Loop over all tracks
  double hMass = m_TagMode == "KSpipi" ? MASS::PI_MASS : MASS::K_MASS;
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    // If track is from KS daughters, skip
    if(std::find(KSDaughterTrackIDs.begin(), KSDaughterTrackIDs.end(), (*Track_iter)->trackId()) != KSDaughterTrackIDs.end()) {
      KSDaughterKalmanTracks.push_back(MDCKalTrack);
      continue;
    }
    // Fill out track information
    if((DTTool.isPion(*Track_iter) && m_TagMode == "KSpipi") ||
       (DTTool.isKaon(*Track_iter) && m_TagMode == "KSKK")) {
      if(MDCKalTrack->charge() == +1) {
	DaughterTrackIterators[HPLUS] = Track_iter;
	KalmanTracks[HPLUS] = MDCKalTrack;
	m_hPlusP = MDCKalTrack->p4(hMass);
	m_DaughterTrackID[2] = (*Track_iter)->trackId();
      } else if(MDCKalTrack->charge() == -1) {
	DaughterTrackIterators[HMINUS] = Track_iter;
	KalmanTracks[HMINUS] = MDCKalTrack;
	m_hMinusP = MDCKalTrack->p4(hMass);
	m_DaughterTrackID[3] = (*Track_iter)->trackId();
      }
    }
  }
  // Do a Kalman kinematic fit of the KS daughter tracks, and constrain the KS mass to its PDG value
  WTrackParameter WTrackKSPIplus(MASS::PI_MASS, KSDaughterKalmanTracks[0]->getZHelix(), KSDaughterKalmanTracks[0]->getZError());
  WTrackParameter WTrackKSPIminus(MASS::PI_MASS, KSDaughterKalmanTracks[1]->getZHelix(), KSDaughterKalmanTracks[1]->getZError());
  KalmanKinematicFit *KSKalmanFit = KalmanKinematicFit::instance();
  KSKalmanFit->init();
  KSKalmanFit->AddTrack(0, WTrackKSPIplus);
  KSKalmanFit->AddTrack(1, WTrackKSPIminus);
  KSKalmanFit->AddResonance(0, MASS::KS_MASS, 0, 1);
  bool KSKalmanFitSuccess = KSKalmanFit->Fit();
  // Do a Kalman kinematic fit of the D daughter tracks, and constrain the KS and D masses to their PDG values
  if(KSKalmanFitSuccess) {
    KSKalmanFit->BuildVirtualParticle(0);
    WTrackParameter WTrackKS = KSKalmanFit->wVirtualTrack(0);
    WTrackParameter WTrackhplus(hMass, KalmanTracks[HPLUS]->getZHelix(), KalmanTracks[HPLUS]->getZError());
    WTrackParameter WTrackhminus(hMass, KalmanTracks[HMINUS]->getZHelix(), KalmanTracks[HMINUS]->getZError());
    KalmanKinematicFit *KalmanFit = KalmanKinematicFit::instance();
    KalmanFit->init();
    KalmanFit->AddTrack(HPLUS, WTrackhplus);
    KalmanFit->AddTrack(HMINUS, WTrackhminus);
    KalmanFit->AddTrack(KSHORT, WTrackKS);
    KalmanFit->AddResonance(0, MASS::D_MASS, HPLUS, HMINUS, KSHORT);
    m_KalmanFitSuccess = KalmanFit->Fit();
    if(m_KalmanFitSuccess) {
      m_KalmanFitChi2 = KalmanFit->chisq();
      m_KShortPKalmanFit = KalmanFit->pfit(KSHORT);
      m_hPlusPKalmanFit = KalmanFit->pfit(HPLUS);
      m_hMinusPKalmanFit = KalmanFit->pfit(HMINUS);
    }
  }
  // For KSpipi, do a vertex fit of the pipi in the tag, to make sure they're not KS in disguise
  if(m_TagMode == "KSpipi") {
    double Mpipi = (m_hPlusP + m_hMinusP).m();
    m_pipiKSFitSuccess = 0;
    if(Mpipi - MASS::KS_MASS < 0.050 && Mpipi - MASS::KS_MASS > -0.060) {
      FindKS pipifindKS(false);
      std::vector<int> PionTrackIDs;
      PionTrackIDs.push_back((*DaughterTrackIterators[HPLUS])->trackId());
      PionTrackIDs.push_back((*DaughterTrackIterators[HMINUS])->trackId());
      StatusCode statuscode = pipifindKS.findKS(DTTool_iter, DTTool, PionTrackIDs);
      m_pipiDecayLengthVeeVertex = pipifindKS.GetDecayLengthVeeVertex();
      m_pipiChi2VeeVertex = pipifindKS.GetChi2VeeVertex();
      m_pipiKSMassVeeVertex = pipifindKS.GetKSMassVeeVertex();
      if(statuscode == StatusCode::SUCCESS) {
	m_pipiKSFitSuccess = 1;
	m_pipiDecayLengthFit = pipifindKS.GetDecayLengthFit();
	m_pipiDecayLengthErrorFit = pipifindKS.GetDecayLengthErrorFit();
	m_pipiChi2Fit = pipifindKS.GetChi2Fit();
      }
    }
  }
  return StatusCode::SUCCESS;
}

std::vector<int> FindKShhTagInfo::GetDaughterTrackID() const {
  return m_DaughterTrackID;
}

int FindKShhTagInfo::GetKSFitSuccess() const {
  return m_KSFitSuccess;
}

double FindKShhTagInfo::GetDecayLengthVeeVertex() const {
  return m_DecayLengthVeeVertex;
}

double FindKShhTagInfo::GetChi2VeeVertex() const {
  return m_Chi2VeeVertex;
}

double FindKShhTagInfo::GetKSMassVeeVertex() const {
  return m_KSMassVeeVertex;
}

double FindKShhTagInfo::GetDecayLengthFit() const {
  return m_DecayLengthFit;
}

double FindKShhTagInfo::GetDecayLengthErrorFit() const {
  return m_DecayLengthErrorFit;
}

double FindKShhTagInfo::GetChi2Fit() const {
  return m_Chi2Fit;
}

double FindKShhTagInfo::GetKSPiPlusP(int i) const {
  return m_KSPiPlusP[i];
}

double FindKShhTagInfo::GetKSPiMinusP(int i) const {
  return m_KSPiMinusP[i];
}

double FindKShhTagInfo::GetKShortP(int i) const {
  return m_KShortP[i];
}

double FindKShhTagInfo::GethPlusP(int i) const {
  return m_hPlusP[i];
}

double FindKShhTagInfo::GethMinusP(int i) const {
  return m_hMinusP[i];
}

int FindKShhTagInfo::GetKalmanFitSuccess() const {
  return m_KalmanFitSuccess;
}

double FindKShhTagInfo::GetKalmanFitChi2() const {
  return m_KalmanFitChi2;
}

double FindKShhTagInfo::GetKShortPKalmanFit(int i) const {
  return m_KShortPKalmanFit[i];
}

double FindKShhTagInfo::GethPlusPKalmanFit(int i) const {
  return m_hPlusPKalmanFit[i];
}

double FindKShhTagInfo::GethMinusPKalmanFit(int i) const {
  return m_hMinusPKalmanFit[i];
}

int FindKShhTagInfo::GetpipiKSFitSuccess() const {
  return m_pipiKSFitSuccess;
}

double FindKShhTagInfo::GetpipiDecayLengthVeeVertex() const {
  return m_pipiDecayLengthVeeVertex;
}

double FindKShhTagInfo::GetpipiChi2VeeVertex() const {
  return m_pipiChi2VeeVertex;
}

double FindKShhTagInfo::GetpipiKSMassVeeVertex() const {
  return m_pipiKSMassVeeVertex;
}

double FindKShhTagInfo::GetpipiDecayLengthFit() const {
  return m_pipiDecayLengthFit;
}

double FindKShhTagInfo::GetpipiDecayLengthErrorFit() const {
  return m_pipiDecayLengthErrorFit;
}

double FindKShhTagInfo::GetpipiChi2Fit() const {
  return m_pipiChi2Fit;
}
