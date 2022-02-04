// Martin Duy Tat 5th March 2021

// KKpipi
#include "KKpipi/FindKKpipiVersusKSpipiTagInfo.h"
#include "KKpipi/FindKS.h"
#include "KKpipi/FindhhTagInfo.h"
#include "KKpipi/ParticleMasses.h"
#include "KKpipi/KKpipiUtilities.h"
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
#include "VertexFit/WTrackParameter.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
// ROOT
#include "TMath.h"
// STL
#include<vector>
#include<string>
// Particle masses
#include "KKpipi/ParticleMasses.h"

FindKKpipiVersusKSpipiTagInfo::FindKKpipiVersusKSpipiTagInfo(): m_DaughterTrackID_KKpipi(std::vector<int>(3)),
				    m_DaughterTrackID_KSpipi(std::vector<int>(4)),
				    m_KSFitSuccess_KSpipi(0),
				    m_DecayLengthVeeVertex_KSpipi(0.0),
				    m_Chi2VeeVertex_KSpipi(0.0),
				    m_KSMassVeeVertex_KSpipi(0.0),
				    m_DecayLengthFit_KSpipi(0.0),
				    m_DecayLengthErrorFit_KSpipi(0.0),
				    m_Chi2Fit_KSpipi(0.0),
				    m_KalmanFitSuccess(0),
				    m_KalmanFitChi2(0.0),
				    m_KSFitSuccess_KKpipi(0),
				    m_DecayLengthVeeVertex_KKpipi(0.0),
				    m_Chi2VeeVertex_KKpipi(0.0),
				    m_KSMassVeeVertex_KKpipi(0.0),
				    m_DecayLengthFit_KKpipi(0.0),
				    m_DecayLengthErrorFit_KKpipi(0.0),
                                    m_Chi2Fit_KKpipi(0.0) {
}

FindKKpipiVersusKSpipiTagInfo::~FindKKpipiVersusKSpipiTagInfo() {
}

StatusCode FindKKpipiVersusKSpipiTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  // First start with tag KSpipi side
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
  WTrackParameters TrackParameters;
  // Loop over all tracks
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    // If track is from KS daughters, skip
    if(std::find(KSDaughterTrackIDs.begin(), KSDaughterTrackIDs.end(), (*Track_iter)->trackId()) != KSDaughterTrackIDs.end()) {
      if(MDCKalTrack->charge() == +1) {
	TrackParameters.TagKSPiPlus = WTrackParameter(MASS::PI_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
      } else {
	TrackParameters.TagKSPiMinus = WTrackParameter(MASS::PI_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
      }
      continue;
    }
    // Fill out track information
    if(DTTool.isPion(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	TrackParameters.TagPiPlus = WTrackParameter(MASS::PI_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
	m_hPlusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID[2] = (*Track_iter)->trackId();
      } else if(MDCKalTrack->charge() == -1) {
	TrackParameters.TagPiMinus = WTrackParameter(MASS::PI_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
	m_hMinusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID[3] = (*Track_iter)->trackId();
      }
    }
  }
  // Then do KKpipi side, with a missing kaon
  bool FoundKPlus = false, FoundKMinus = false, FoundPiPlus = false, FoundPiMinus = false;
  SmartRefVector<EvtRecTrack> OtherTracks = (*DTTool_iter)->otherTracks();
  // Go through all tracks and find the daughter tracks
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    if(DTTool.isKaon(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	FoundKPlus ? return StatusCode::FAILURE : FoundKPlus = true;
	TrackParameters.TagSignalKaon = WTrackParameter(MASS::K_MASS, MDCKalTrack->getZHelixK(), MDCKalTrack->getZErrorK());
	m_KPlusP = MDCKalTrack->p4(MASS::K_MASS);
	m_DaughterTrackID_KKpipi[0] = (*Track_iter)->trackId();
      } else if (MDCKalTrack->charge() == -1) {
	FoundKMinus ? return StatusCode::FAILURE : FoundKMinus = true;
	TrackParameters.TagSignalKaon = WTrackParameter(MASS::K_MASS, MDCKalTrack->getZHelixK(), MDCKalTrack->getZErrorK());
	m_KMinusP = MDCKalTrack->p4(MASS::K_MASS);
	m_DaughterTrackID_KKpipi[0] = (*Track_iter)->trackId();
      }
    } else if(DTTool.isPion(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	FoundPiPlus ? return StatusCode::FAILURE : FoundPiPlus = true;
	TrackParameters.TagSignalPiPlus = WTrackParameter(MASS::PI_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
	m_PiPlusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID_KKpipi[1] = (*Track_iter)->trackId();
      } else if(MDCKalTrack->charge() == -1) {
	FoundPiMinus ? return StatusCode::FAILURE : FoundPiMinus = true;
	TrackParameters.TagSignalPiMinus = WTrackParameter(MASS::PI_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
	m_PiMinusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID_KKpipi[2] = (*Track_iter)->trackId();
      }
    }
  }
  // Check if pions are found
  if(!FoundPiPlus || !FoundPiMinus) {
    return StatusCode::FAILURE;
  }
  // Check if a single kaon is found
  if((FoundKPlus && FoundKMinus) || (!FoundKPlus && !FoundKMinus)) {
    return StatusCode::FAILURE;
  }
  m_RecKCharge = FoundKPlus ? +1 : -1;
  if(mRecKCharge == +1) {
    CLHEP::HepLorentzVector P_X = m_KPlusP + m_PiPlusP + m_PiMinusP;
    m_KMinusP = KKpipiUtilities::GetMissingMomentum((*DTTool_iter)->p4(), P_X, (*DTTool_iter)->beamE());
    m_MissingMass2 = m_KMinusP.m2();
  } else {
    CLHEP::HepLorentzVector P_X = m_KMinusP + m_PiPlusP + m_PiMinusP;
    m_KPlusP = KKpipiUtilities::GetMissingMomentum((*DTTool_iter)->p4(), P_X, (*DTTool_iter)->beamE());
    m_MissingMass2 = m_KMinusP.m2();
  }
  // Finally do the Kalman kinematic fit
  DoKalmanFit(TrackParameters, m_RecKCharge);
  return StatusCode::SUCCESS;
}

void FindKKpipiVersusKSpipiTagInfo::DoKalmanFit(const WTrackParameters &TrackParameters, int RecKCharge) {    
  // Do a Kalman kinematic fit of the KS daughter tracks, and constrain the KS mass to its PDG value
  KalmanKinematicFit *KSKalmanFit = KalmanKinematicFit::instance();
  KSKalmanFit->init();
  KSKalmanFit->AddTrack(0, TrackParameters.TagKSPiPlus);
  KSKalmanFit->AddTrack(1, TrackParameters.TagKSPiMinus);
  KSKalmanFit->AddResonance(0, MASS::KS_MASS, 0, 1);
  bool KSKalmanFitSuccess = KSKalmanFit->Fit();
  // Do a Kalman kinematic fit of the D daughter tracks, and constrain the KS and D masses to their PDG values
  if(KSKalmanFitSuccess) {
    KSKalmanFit->BuildVirtualParticle(0);
    WTrackParameter WTrackKS = KSKalmanFit->wVirtualTrack(0);
    KalmanKinematicFit *KalmanFit = KalmanKinematicFit::instance();
    KalmanFit->init();
    KalmanFit->AddTrack(0, WTrackKS);
    KalmanFit->AddTrack(1, TrackParameters.TagPiPlus);
    KalmanFit->AddTrack(2, TrackParameters.TagPiMinus);
    KalmanFit->AddTrack(3, TrackParameters.SignalKaon);
    KalmanFit->AddMissTrack(4, MASS::K_MASS);
    KalmanFit->AddTrack(5, TrackParameters.SignalPiPlus);
    KalmanFit->AddTrack(6, TrackParameters.SignalPiMinus);
    KalmanFit->AddResonance(0, MASS::D_MASS, 0, 1, 2);
    KalmanFit->AddResonance(1, MASS::D_MASS, 3, 4, 5, 6);
    // Constrain total four-momentum to psi(3770) in lab frame
    CLHEP::HepLorentzVector Psipp_P(0.0, 0.0, 0.0, MASS::JPSI_MASS);
    // Crossing angle of 11 mrad
    Psipp_P.boost(0.011, 0.0, 0.0);
    KalmanFit->AddFourMomentum(2, Psipp_P);
    m_KalmanFitSuccess = KalmanFit->Fit();
    if(m_KalmanFitSuccess) {
      m_KalmanFitChi2 = KalmanFit->chisq();
      m_KShortPKalmanFit = KalmanFit->pfit(0);
      m_hPlusPKalmanFit = KalmanFit->pfit(1);
      m_hMinusPKalmanFit = KalmanFit->pfit(2);
      if(RecKCharge == +1) {
	m_KPlusPKalmanFit = KalmanFit->pfit(3);
	m_KMinusPKalmanFit = KalmanFit->pfit(4);
      } else {
	m_KPlusPKalmanFit = KalmanFit->pfit(4);
	m_KMinusPKalmanFit = KalmanFit->pfit(3);
      }
      m_PiPlusPKalmanFit = KalmanFit->pfit(5);
      m_PiMinusPKalmanFit = KalmanFit->pfit(6);
    }
  }
  
}

std::vector<int> FindKKpipiVersusKSpipiTagInfo::GetDaughterTrackID_KKpipi() const {
  return m_DaughterTrackID_KKpipi;
}

std::vector<int> FindKKpipiVersusKSpipiTagInfo::GetDaughterTrackID_KSpipi() const {
  return m_DaughterTrackID_KSpipi;
}

int FindKKpipiVersusKSpipiTagInfo::GetKSFitSuccess_KSpipi() const {
  return m_KSFitSuccesss_KSpipi;
}

double FindKKpipiVersusKSpipiTagInfo::GetDecayLengthVeeVertexs_KSpipi() const {
  return m_DecayLengthVeeVertexs_KSpipi;
}

double FindKKpipiVersusKSpipiTagInfo::GetChi2VeeVertexs_KSpipi() const {
  return m_Chi2VeeVertexs_KSpipi;
}

double FindKKpipiVersusKSpipiTagInfo::GetKSMassVeeVertexs_KSpipi() const {
  return m_KSMassVeeVertexs_KSpipi;
}

double FindKKpipiVersusKSpipiTagInfo::GetDecayLengthFits_KSpipi() const {
  return m_DecayLengthFits_KSpipi;
}

double FindKKpipiVersusKSpipiTagInfo::GetDecayLengthErrorFits_KSpipi() const {
  return m_DecayLengthErrorFits_KSpipi;
}

double FindKKpipiVersusKSpipiTagInfo::GetChi2Fits_KSpipi() const {
  return m_Chi2Fits_KSpipi;
}

double FindKKpipiVersusKSpipiTagInfo::GetKSPiPlusP(int i) const {
  return m_KSPiPlusP[i];
}

double FindKKpipiVersusKSpipiTagInfo::GetKSPiMinusP(int i) const {
  return m_KSPiMinusP[i];
}

double FindKKpipiVersusKSpipiTagInfo::GetKShortP(int i) const {
  return m_KShortP[i];
}

double FindKKpipiVersusKSpipiTagInfo::GethPlusP(int i) const {
  return m_hPlusP[i];
}

double FindKKpipiVersusKSpipiTagInfo::GethMinusP(int i) const {
  return m_hMinusP[i];
}

int FindKKpipiVersusKSpipiTagInfo::GetKalmanFitSuccess() const {
  return m_KalmanFitSuccess;
}

double FindKKpipiVersusKSpipiTagInfo::GetKalmanFitChi2() const {
  return m_KalmanFitChi2;
}

double FindKKpipiVersusKSpipiTagInfo::GetKShortPKalmanFit(int i) const {
  return m_KShortPKalmanFit[i];
}

double FindKKpipiVersusKSpipiTagInfo::GethPlusPKalmanFit(int i) const {
  return m_hPlusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiTagInfo::GethMinusPKalmanFit(int i) const {
  return m_hMinusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiTagInfo::GetKPlusP(int i) const {
  return m_KPlusP[i];
}

double FindKKpipiVersusKSpipiTagInfo::GetKMinusP(int i) const {
  return m_KMinusP[i];
}

double FindKKpipiVersusKSpipiTagInfo::GetPiPlusP(int i) const {
  return m_PiPlusP[i];
}

double FindKKpipiVersusKSpipiTagInfo::GetPiMinusP(int i) const {
  return m_PiMinusP[i];
}

double FindKKpipiVersusKSpipiTagInfo::GetKPlusPKalmanFit(int i) const {
  return m_KPlusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiTagInfo::GetKMinusPKalmanFit(int i) const {
  return m_KMinusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiTagInfo::GetPiPlusPKalmanFit(int i) const {
  return m_PiPlusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiTagInfo::GetPiMinusPKalmanFit(int i) const {
  return m_PiMinusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiTagInfo::GetMpipi() const {
  return (m_PiPlusP + m_PiMinusP).m();
}

int FindKKpipiVersusKKpipiTagInfo::GetKSFitSuccess_KKpipi() const {
  return m_KSFitSuccesss_KKpipi;
}

double FindKKpipiVersusKKpipiTagInfo::GetDecayLengthVeeVertexs_KKpipi() const {
  return m_DecayLengthVeeVertexs_KKpipi;
}

double FindKKpipiVersusKKpipiTagInfo::GetChi2VeeVertexs_KKpipi() const {
  return m_Chi2VeeVertexs_KKpipi;
}

double FindKKpipiVersusKKpipiTagInfo::GetKSMassVeeVertexs_KKpipi() const {
  return m_KSMassVeeVertexs_KKpipi;
}

double FindKKpipiVersusKKpipiTagInfo::GetDecayLengthFits_KKpipi() const {
  return m_DecayLengthFits_KKpipi;
}

double FindKKpipiVersusKKpipiTagInfo::GetDecayLengthErrorFits_KKpipi() const {
  return m_DecayLengthErrorFits_KKpipi;
}

double FindKKpipiVersusKKpipiTagInfo::GetChi2Fits_KKpipi() const {
  return m_Chi2Fits_KKpipi;
}

double FindKKpipiVersusKKpipiTagInfo::GetMissingMass2() const {
  return m_MissingMass2;
}

int FindKKpipiVersusKKpipiTagInfo::GetRecKCharge() const {
  return m_RecKCharge;
}
