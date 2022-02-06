// Martin Duy Tat 2nd February 2022

// KKpipi
#include "KKpipi/FindKKpipiVersusKSpipiPartRecoTagInfo.h"
#include "KKpipi/FindKS.h"
#include "KKpipi/FindhhTagInfo.h"
#include "KKpipi/ParticleMasses.h"
#include "KKpipi/KKpipiUtilities.h"
// Gaudi
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/SmartDataPtr.h"
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

FindKKpipiVersusKSpipiPartRecoTagInfo::FindKKpipiVersusKSpipiPartRecoTagInfo(): m_DaughterTrackID_KKpipi(std::vector<int>(3)),
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

FindKKpipiVersusKSpipiPartRecoTagInfo::~FindKKpipiVersusKSpipiPartRecoTagInfo() {
}

StatusCode FindKKpipiVersusKSpipiPartRecoTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  // First start with tag KSpipi side
  // Find the actual KS candidate
  FindKS findKS(true);
  StatusCode status = findKS.findKS(DTTool_iter, DTTool);
  if(status == StatusCode::SUCCESS) {
    m_KSFitSuccess_KSpipi = 1;
    m_DecayLengthFit_KSpipi = findKS.GetDecayLengthFit();
    m_DecayLengthErrorFit_KSpipi = findKS.GetDecayLengthErrorFit();
    m_Chi2Fit_KSpipi = findKS.GetChi2Fit();
  } else {
    m_KSFitSuccess_KSpipi = 0;
  }
  // Fill out the information about the KS candidate
  m_DecayLengthVeeVertex_KSpipi = findKS.GetDecayLengthVeeVertex();
  m_Chi2VeeVertex_KSpipi = findKS.GetChi2VeeVertex();
  m_KSMassVeeVertex_KSpipi = findKS.GetKSMassVeeVertex();
  m_KSPiPlusP = CLHEP::HepLorentzVector(findKS.GetKSPiPlusP(0), findKS.GetKSPiPlusP(1), findKS.GetKSPiPlusP(2), findKS.GetKSPiPlusP(3));
  m_KSPiMinusP = CLHEP::HepLorentzVector(findKS.GetKSPiMinusP(0), findKS.GetKSPiMinusP(1), findKS.GetKSPiMinusP(2), findKS.GetKSPiMinusP(3));
  m_KShortP = findKS.GetKShortPFit();
  // Get the track ID of the KS daughters
  std::vector<int> KSDaughterTrackIDs = findKS.GetDaughterTrackIDs();
  m_DaughterTrackID_KSpipi[0] = KSDaughterTrackIDs[0];
  m_DaughterTrackID_KSpipi[1] = KSDaughterTrackIDs[1];
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
	m_DaughterTrackID_KSpipi[2] = (*Track_iter)->trackId();
      } else if(MDCKalTrack->charge() == -1) {
	TrackParameters.TagPiMinus = WTrackParameter(MASS::PI_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
	m_hMinusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID_KSpipi[3] = (*Track_iter)->trackId();
      }
    }
  }
  // Then do KKpipi side, with a missing kaon
  bool FoundKPlus = false, FoundKMinus = false, FoundPiPlus = false, FoundPiMinus = false;
  std::vector<int> PionTrackIDs(2);
  SmartRefVector<EvtRecTrack> OtherTracks = (*DTTool_iter)->otherTracks();
  // Go through all tracks and find the daughter tracks
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = OtherTracks.begin(); Track_iter != OtherTracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    if(DTTool.isKaon(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	if(FoundKPlus) {
	  return StatusCode::FAILURE;
	} else {
	  FoundKPlus = true;
	}
	TrackParameters.SignalKaon = WTrackParameter(MASS::K_MASS, MDCKalTrack->getZHelixK(), MDCKalTrack->getZErrorK());
	m_KPlusP = MDCKalTrack->p4(MASS::K_MASS);
	m_DaughterTrackID_KKpipi[0] = (*Track_iter)->trackId();
      } else if (MDCKalTrack->charge() == -1) {
	if(FoundKMinus) {
	  return StatusCode::FAILURE;
	} else {
	  FoundKMinus = true;
	}
	TrackParameters.SignalKaon = WTrackParameter(MASS::K_MASS, MDCKalTrack->getZHelixK(), MDCKalTrack->getZErrorK());
	m_KMinusP = MDCKalTrack->p4(MASS::K_MASS);
	m_DaughterTrackID_KKpipi[0] = (*Track_iter)->trackId();
      }
    } else if(DTTool.isPion(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	if(FoundPiPlus) {
	  return StatusCode::FAILURE;
	} else {
	  FoundPiPlus = true;
	}
	TrackParameters.SignalPiPlus = WTrackParameter(MASS::PI_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
	m_PiPlusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID_KKpipi[1] = (*Track_iter)->trackId();
	PionTrackIDs[0] = (*Track_iter)->trackId();
      } else if(MDCKalTrack->charge() == -1) {
	if(FoundPiMinus) {
	  return StatusCode::FAILURE;
	} else {
	  FoundPiMinus = true;
	}
	TrackParameters.SignalPiMinus = WTrackParameter(MASS::PI_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
	m_PiMinusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID_KKpipi[2] = (*Track_iter)->trackId();
	PionTrackIDs[1] = (*Track_iter)->trackId();
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
  if(m_RecKCharge == +1) {
    CLHEP::HepLorentzVector P_X = m_KPlusP + m_PiPlusP + m_PiMinusP;
    m_KMinusP = KKpipiUtilities::GetMissingMomentum((*DTTool_iter)->p4(), P_X, (*DTTool_iter)->beamE());
    m_MissingMass2 = m_KMinusP.m2();
  } else {
    CLHEP::HepLorentzVector P_X = m_KMinusP + m_PiPlusP + m_PiMinusP;
    m_KPlusP = KKpipiUtilities::GetMissingMomentum((*DTTool_iter)->p4(), P_X, (*DTTool_iter)->beamE());
    m_MissingMass2 = m_KPlusP.m2();
  }
  double Mpipi = (m_PiPlusP + m_PiMinusP).m();
  m_KSFitSuccess_KKpipi = 0;
  // Check if the \f$\pi\pi\f$ pair is a \f$K_S\f$ in disguise
  if(Mpipi - MASS::KS_MASS < 0.050 && Mpipi - MASS::KS_MASS > -0.060) {
    FindKS findKS(false);
    StatusCode statuscode = findKS.findKS(DTTool_iter, DTTool, PionTrackIDs);
    m_DecayLengthVeeVertex_KKpipi = findKS.GetDecayLengthVeeVertex();
    m_Chi2VeeVertex_KKpipi = findKS.GetChi2VeeVertex();
    m_KSMassVeeVertex_KKpipi = findKS.GetKSMassVeeVertex();
    if(statuscode == StatusCode::SUCCESS) {
      m_KSFitSuccess_KKpipi = 1;
      m_DecayLengthFit_KKpipi = findKS.GetDecayLengthFit();
      m_DecayLengthErrorFit_KKpipi = findKS.GetDecayLengthErrorFit();
      m_Chi2Fit_KKpipi = findKS.GetChi2Fit();
    }
  }
  // Find the pi0
  m_NumberPi0 = FindPi0();
  // Finally do the Kalman kinematic fit
  DoKalmanFit(TrackParameters, m_RecKCharge);
  return StatusCode::SUCCESS;
}

void FindKKpipiVersusKSpipiPartRecoTagInfo::DoKalmanFit(const WTrackParameters &TrackParameters, int RecKCharge) {    
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

int FindKKpipiVersusKSpipiPartRecoTagInfo::FindPi0() {
  // Prepare event data service
  IDataProviderSvc *EventDataService = nullptr;
  Gaudi::svcLocator()->service("EventDataSvc", EventDataService);
  // Prepare pi0 service
  SmartDataPtr<EvtRecPi0Col> evtRecPi0Col(EventDataService, "/Event/EvtRec/EvtRecPi0Col");
  // Look for pi0
  int NumberPi0 = 0;
  for(EvtRecPi0Col::iterator Pi0_iter = evtRecPi0Col->begin(); Pi0_iter != evtRecPi0Col->end(); Pi0_iter++) {
    // Get photon tracks...?
    EvtRecTrack *HighEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Pi0_iter)->hiEnGamma());
    EvtRecTrack *LowEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Pi0_iter)->loEnGamma());
    // Get EM shower four-momenta of photons
    RecEmcShower *HighEPhotonShower = HighEnergyPhotonTrack->emcShower();
    RecEmcShower *LowEPhotonShower = LowEnergyPhotonTrack->emcShower();
    // Pi0 invariant mass
    CLHEP::HepLorentzVector HighEPhotonP = KKpipiUtilities::GetPhoton4Vector(HighEPhotonShower->energy(), HighEPhotonShower->theta(), HighEPhotonShower->phi());
    CLHEP::HepLorentzVector LowEPhotonP = KKpipiUtilities::GetPhoton4Vector(LowEPhotonShower->energy(), LowEPhotonShower->theta(), LowEPhotonShower->phi());
    double Mgammagamma = (HighEPhotonP + LowEPhotonP).m();
    if(Mgammagamma > 0.115 && Mgammagamma < 0.150) {
      NumberPi0++;
    }
  }
  return NumberPi0;
}

std::vector<int> FindKKpipiVersusKSpipiPartRecoTagInfo::GetDaughterTrackID_KKpipi() const {
  return m_DaughterTrackID_KKpipi;
}

std::vector<int> FindKKpipiVersusKSpipiPartRecoTagInfo::GetDaughterTrackID_KSpipi() const {
  return m_DaughterTrackID_KSpipi;
}

int FindKKpipiVersusKSpipiPartRecoTagInfo::GetKSFitSuccess_KSpipi() const {
  return m_KSFitSuccess_KSpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetDecayLengthVeeVertex_KSpipi() const {
  return m_DecayLengthVeeVertex_KSpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetChi2VeeVertex_KSpipi() const {
  return m_Chi2VeeVertex_KSpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetKSMassVeeVertex_KSpipi() const {
  return m_KSMassVeeVertex_KSpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetDecayLengthFit_KSpipi() const {
  return m_DecayLengthFit_KSpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetDecayLengthErrorFit_KSpipi() const {
  return m_DecayLengthErrorFit_KSpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetChi2Fit_KSpipi() const {
  return m_Chi2Fit_KSpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetKSPiPlusP(int i) const {
  return m_KSPiPlusP[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetKSPiMinusP(int i) const {
  return m_KSPiMinusP[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetKShortP(int i) const {
  return m_KShortP[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GethPlusP(int i) const {
  return m_hPlusP[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GethMinusP(int i) const {
  return m_hMinusP[i];
}

int FindKKpipiVersusKSpipiPartRecoTagInfo::GetKalmanFitSuccess() const {
  return m_KalmanFitSuccess;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetKalmanFitChi2() const {
  return m_KalmanFitChi2;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetKShortPKalmanFit(int i) const {
  return m_KShortPKalmanFit[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GethPlusPKalmanFit(int i) const {
  return m_hPlusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GethMinusPKalmanFit(int i) const {
  return m_hMinusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetKPlusP(int i) const {
  return m_KPlusP[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetKMinusP(int i) const {
  return m_KMinusP[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetPiPlusP(int i) const {
  return m_PiPlusP[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetPiMinusP(int i) const {
  return m_PiMinusP[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetKPlusPKalmanFit(int i) const {
  return m_KPlusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetKMinusPKalmanFit(int i) const {
  return m_KMinusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetPiPlusPKalmanFit(int i) const {
  return m_PiPlusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetPiMinusPKalmanFit(int i) const {
  return m_PiMinusPKalmanFit[i];
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetMpipi_KKpipi() const {
  return (m_PiPlusP + m_PiMinusP).m();
}

int FindKKpipiVersusKSpipiPartRecoTagInfo::GetKSFitSuccess_KKpipi() const {
  return m_KSFitSuccess_KKpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetDecayLengthVeeVertex_KKpipi() const {
  return m_DecayLengthVeeVertex_KKpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetChi2VeeVertex_KKpipi() const {
  return m_Chi2VeeVertex_KKpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetKSMassVeeVertex_KKpipi() const {
  return m_KSMassVeeVertex_KKpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetDecayLengthFit_KKpipi() const {
  return m_DecayLengthFit_KKpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetDecayLengthErrorFit_KKpipi() const {
  return m_DecayLengthErrorFit_KKpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetChi2Fit_KKpipi() const {
  return m_Chi2Fit_KKpipi;
}

double FindKKpipiVersusKSpipiPartRecoTagInfo::GetMissingMass2() const {
  return m_MissingMass2;
}

int FindKKpipiVersusKSpipiPartRecoTagInfo::GetRecKCharge() const {
  return m_RecKCharge;
}

int FindKKpipiVersusKSpipiPartRecoTagInfo::GetNumberPi0() const {
  return m_NumberPi0;
}
