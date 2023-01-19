// Martin Duy Tat 19th January 2023

// KKpipi
#include "KKpipi/FindKKpipiVersusKKPartRecoTagInfo.h"
#include "KKpipi/ParticleMasses.h"
#include "KKpipi/KKpipiUtilities.h"
#include "KKpipi/FindhhTagInfo.h"
// Gaudi
#include "GaudiKernel/Bootstrap.h"
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

FindKKpipiVersusKKPartRecoTagInfo::FindKKpipiVersusKKPartRecoTagInfo(): m_DaughterTrackID_KKpipi(std::vector<int>(3)),
                                    m_FindKKTagInfo("KK"),
				    m_KalmanFitSuccess(0),
				    m_KalmanFitChi2(0.0),
				    m_KSFitSuccess_KKpipi(0),
				    m_DecayLengthVeeVertex_KKpipi(0.0),
				    m_Chi2VeeVertex_KKpipi(0.0),
				    m_KSMassVeeVertex_KKpipi(0.0),
				    m_DecayLengthFit_KKpipi(0.0),
				    m_DecayLengthErrorFit_KKpipi(0.0),
                                    m_Chi2Fit_KKpipi(0.0),
                                    m_MissingMass2(0.0),
                                    m_RecKCharge(0),
                                    m_NumberPi0(0) {
}

FindKKpipiVersusKKPartRecoTagInfo::~FindKKpipiVersusKKPartRecoTagInfo() {
}

StatusCode FindKKpipiVersusKKPartRecoTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  // First start with tag KK side
  // Find the two kaons
  m_FindKKTagInfo = FindhhTagInfo("KK");
  StatusCode status = m_FindKKTagInfo.CalculateTagInfo(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  // Get all tracks
  WTrackParameters TrackParameters;
  // Loop over all tracks
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    // If track is from pi+ and pi-, save track information for Kalman fit
    if(DTTool.isKaon(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	TrackParameters.TagKPlus = WTrackParameter(MASS::K_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
      } else if(MDCKalTrack->charge() == -1) {
	TrackParameters.TagKMinus = WTrackParameter(MASS::K_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
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

void FindKKpipiVersusKKPartRecoTagInfo::DoKalmanFit(const WTrackParameters &TrackParameters, int RecKCharge) {    
  // Do a Kalman kinematic fit of the D daughter tracks, and constrain the KS and D masses to their PDG values
  KalmanKinematicFit *KalmanFit = KalmanKinematicFit::instance();
  KalmanFit->init();
  KalmanFit->AddTrack(0, TrackParameters.TagKPlus);
  KalmanFit->AddTrack(1, TrackParameters.TagKMinus);
  KalmanFit->AddTrack(2, TrackParameters.SignalKaon);
  KalmanFit->AddMissTrack(3, MASS::K_MASS);
  KalmanFit->AddTrack(4, TrackParameters.SignalPiPlus);
  KalmanFit->AddTrack(5, TrackParameters.SignalPiMinus);
  KalmanFit->AddResonance(0, MASS::D_MASS, 0, 1);
  KalmanFit->AddResonance(1, MASS::D_MASS, 2, 3, 4, 5);
  // Constrain total four-momentum to psi(3770) in lab frame
  CLHEP::HepLorentzVector Psipp_P(0.0, 0.0, 0.0, MASS::JPSI_MASS);
  // Crossing angle of 11 mrad
  Psipp_P.boost(0.011, 0.0, 0.0);
  KalmanFit->AddFourMomentum(2, Psipp_P);
  m_KalmanFitSuccess = KalmanFit->Fit();
  if(m_KalmanFitSuccess) {
    m_KalmanFitChi2 = KalmanFit->chisq();
    if(RecKCharge == +1) {
      m_KPlusPKalmanFit = KalmanFit->pfit(2);
      m_KMinusPKalmanFit = KalmanFit->pfit(3);
    } else {
      m_KPlusPKalmanFit = KalmanFit->pfit(3);
      m_KMinusPKalmanFit = KalmanFit->pfit(2);
    }
    m_PiPlusPKalmanFit = KalmanFit->pfit(4);
    m_PiMinusPKalmanFit = KalmanFit->pfit(5);
  }
}

int FindKKpipiVersusKKPartRecoTagInfo::FindPi0() {
  // Prepare event data service
  IDataProviderSvc *EventDataService = nullptr;
  Gaudi::svcLocator()->service("EventDataSvc", EventDataService);
  // Prepare pi0 service
  SmartDataPtr<EvtRecPi0Col> evtRecPi0Col(EventDataService, "/Event/EvtRec/EvtRecPi0Col");
  // Look for pi0
  int NumberPi0 = 0;
  for(EvtRecPi0Col::iterator Pi0_iter = evtRecPi0Col->begin(); Pi0_iter != evtRecPi0Col->end(); Pi0_iter++) {
    // Get photon tracks...? Also check if these are the same as those from the KK
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

std::vector<int> FindKKpipiVersusKKPartRecoTagInfo::GetDaughterTrackID_KKpipi() const {
  return m_DaughterTrackID_KKpipi;
}

int FindKKpipiVersusKKPartRecoTagInfo::GetKalmanFitSuccess() const {
  return m_KalmanFitSuccess;
}

double FindKKpipiVersusKKPartRecoTagInfo::GetKalmanFitChi2() const {
  return m_KalmanFitChi2;
}

double FindKKpipiVersusKKPartRecoTagInfo::GetKPlusP(int i) const {
  return m_KPlusP[i];
}

double FindKKpipiVersusKKPartRecoTagInfo::GetKMinusP(int i) const {
  return m_KMinusP[i];
}

double FindKKpipiVersusKKPartRecoTagInfo::GetPiPlusP(int i) const {
  return m_PiPlusP[i];
}

double FindKKpipiVersusKKPartRecoTagInfo::GetPiMinusP(int i) const {
  return m_PiMinusP[i];
}

double FindKKpipiVersusKKPartRecoTagInfo::GetKPlusPKalmanFit(int i) const {
  return m_KPlusPKalmanFit[i];
}

double FindKKpipiVersusKKPartRecoTagInfo::GetKMinusPKalmanFit(int i) const {
  return m_KMinusPKalmanFit[i];
}

double FindKKpipiVersusKKPartRecoTagInfo::GetPiPlusPKalmanFit(int i) const {
  return m_PiPlusPKalmanFit[i];
}

double FindKKpipiVersusKKPartRecoTagInfo::GetPiMinusPKalmanFit(int i) const {
  return m_PiMinusPKalmanFit[i];
}

double FindKKpipiVersusKKPartRecoTagInfo::GetMpipi_KKpipi() const {
  return (m_PiPlusP + m_PiMinusP).m();
}

int FindKKpipiVersusKKPartRecoTagInfo::GetKSFitSuccess_KKpipi() const {
  return m_KSFitSuccess_KKpipi;
}

double FindKKpipiVersusKKPartRecoTagInfo::GetDecayLengthVeeVertex_KKpipi() const {
  return m_DecayLengthVeeVertex_KKpipi;
}

double FindKKpipiVersusKKPartRecoTagInfo::GetChi2VeeVertex_KKpipi() const {
  return m_Chi2VeeVertex_KKpipi;
}

double FindKKpipiVersusKKPartRecoTagInfo::GetKSMassVeeVertex_KKpipi() const {
  return m_KSMassVeeVertex_KKpipi;
}

double FindKKpipiVersusKKPartRecoTagInfo::GetDecayLengthFit_KKpipi() const {
  return m_DecayLengthFit_KKpipi;
}

double FindKKpipiVersusKKPartRecoTagInfo::GetDecayLengthErrorFit_KKpipi() const {
  return m_DecayLengthErrorFit_KKpipi;
}

double FindKKpipiVersusKKPartRecoTagInfo::GetChi2Fit_KKpipi() const {
  return m_Chi2Fit_KKpipi;
}

double FindKKpipiVersusKKPartRecoTagInfo::GetMissingMass2() const {
  return m_MissingMass2;
}

int FindKKpipiVersusKKPartRecoTagInfo::GetRecKCharge() const {
  return m_RecKCharge;
}

int FindKKpipiVersusKKPartRecoTagInfo::GetNumberPi0() const {
  return m_NumberPi0;
}

const FindhhTagInfo& FindKKpipiVersusKKPartRecoTagInfo::GetKKTagInfo() const {
  return m_FindKKTagInfo;
}
