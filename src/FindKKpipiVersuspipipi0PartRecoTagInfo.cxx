// Martin Duy Tat 19th January 2023

// KKpipi
#include "KKpipi/FindKKpipiVersuspipipi0PartRecoTagInfo.h"
#include "KKpipi/FindKS.h"
#include "KKpipi/FindPi0Eta.h"
#include "KKpipi/ParticleMasses.h"
#include "KKpipi/KKpipiUtilities.h"
#include "KKpipi/FindPi0Eta.h"
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

FindKKpipiVersuspipipi0PartRecoTagInfo::FindKKpipiVersuspipipi0PartRecoTagInfo(): m_DaughterTrackID_KKpipi(std::vector<int>(3)),
				    m_KSFitSuccess_pipipi0(0),
				    m_DecayLengthVeeVertex_pipipi0(0.0),
				    m_Chi2VeeVertex_pipipi0(0.0),
				    m_KSMassVeeVertex_pipipi0(0.0),
				    m_DecayLengthFit_pipipi0(0.0),
				    m_DecayLengthErrorFit_pipipi0(0.0),
				    m_Chi2Fit_pipipi0(0.0),
                                    m_FindpipiTagInfo("pipi"),
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

FindKKpipiVersuspipipi0PartRecoTagInfo::~FindKKpipiVersuspipipi0PartRecoTagInfo() {
}

StatusCode FindKKpipiVersuspipipi0PartRecoTagInfo::CalculateTagInfo(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  // First start with tag pipipi0 side
  // Find the two pions
  m_FindpipiTagInfo = FindhhTagInfo("pipi");
  StatusCode status = m_FindpipiTagInfo.CalculateTagInfo(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  // Find KS from potential KSpi0 decay
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  double Mpipi = m_FindpipiTagInfo.GetMhh();
  m_KSFitSuccess_pipipi0 = 0;
  if(Mpipi - MASS::KS_MASS < 0.050 && Mpipi - MASS::KS_MASS > -0.060) {
    FindKS findKS(false);
    std::vector<int> PionTrackIDs;
    for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
      if(DTTool.isPion(*Track_iter)) {
	PionTrackIDs.push_back((*Track_iter)->trackId());
      }
    }
    status = findKS.findKS(DTTool_iter, DTTool, PionTrackIDs);
    m_DecayLengthVeeVertex_pipipi0 = findKS.GetDecayLengthVeeVertex();
    m_Chi2VeeVertex_pipipi0 = findKS.GetChi2VeeVertex();
    m_KSMassVeeVertex_pipipi0 = findKS.GetKSMassVeeVertex();
    if(status == StatusCode::SUCCESS) {
      m_KSFitSuccess_pipipi0 = 1;
      m_DecayLengthFit_pipipi0 = findKS.GetDecayLengthFit();
      m_DecayLengthErrorFit_pipipi0 = findKS.GetDecayLengthErrorFit();
      m_Chi2Fit_pipipi0 = findKS.GetChi2Fit();
    }
  }
  // Get all tracks
  WTrackParameters TrackParameters;
  // Loop over all tracks
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
    RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
    // If track is from pi+ and pi-, save track information for Kalman fit
    if(DTTool.isPion(*Track_iter)) {
      if(MDCKalTrack->charge() == +1) {
	TrackParameters.TagPiPlus = WTrackParameter(MASS::PI_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
      } else if(MDCKalTrack->charge() == -1) {
	TrackParameters.TagPiMinus = WTrackParameter(MASS::PI_MASS, MDCKalTrack->getZHelix(), MDCKalTrack->getZError());
      }
    }
  }
  // Find the pi0
  m_FindPi0 = FindPi0Eta();
  status = m_FindPi0.findPi0Eta(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return StatusCode::FAILURE;
  }
  TrackParameters.TagHighEShower = m_FindPi0.GetHighEShower();
  TrackParameters.TagLowEShower = m_FindPi0.GetLowEShower();
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
  Mpipi = (m_PiPlusP + m_PiMinusP).m();
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

void FindKKpipiVersuspipipi0PartRecoTagInfo::DoKalmanFit(const WTrackParameters &TrackParameters, int RecKCharge) {    
  // Do a Kalman kinematic fit of the pi0 showers, and constrain the pi0 mass to its PDG value
  KalmanKinematicFit *Pi0KalmanFit = KalmanKinematicFit::instance();
  Pi0KalmanFit->init();
  Pi0KalmanFit->AddTrack(0, 0.0, TrackParameters.TagHighEShower);
  Pi0KalmanFit->AddTrack(1, 0.0, TrackParameters.TagLowEShower);
  Pi0KalmanFit->AddResonance(0, MASS::PI0_MASS, 0, 1);
  bool Pi0KalmanFitSuccess = Pi0KalmanFit->Fit();
  Pi0KalmanFit->BuildVirtualParticle(0);
  WTrackParameter WTrackPi0 = Pi0KalmanFit->wVirtualTrack(0);
  // Do a Kalman kinematic fit of the D daughter tracks, and constrain the KS and D masses to their PDG values
  if(Pi0KalmanFitSuccess) {
    KalmanKinematicFit *KalmanFit = KalmanKinematicFit::instance();
    KalmanFit->init();
    KalmanFit->AddTrack(0, TrackParameters.TagPiPlus);
    KalmanFit->AddTrack(1, TrackParameters.TagPiMinus);
    KalmanFit->AddTrack(2, WTrackPi0);
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

int FindKKpipiVersuspipipi0PartRecoTagInfo::FindPi0() {
  // Prepare event data service
  IDataProviderSvc *EventDataService = nullptr;
  Gaudi::svcLocator()->service("EventDataSvc", EventDataService);
  // Prepare pi0 service
  SmartDataPtr<EvtRecPi0Col> evtRecPi0Col(EventDataService, "/Event/EvtRec/EvtRecPi0Col");
  // Look for pi0
  int NumberPi0 = 0;
  for(EvtRecPi0Col::iterator Pi0_iter = evtRecPi0Col->begin(); Pi0_iter != evtRecPi0Col->end(); Pi0_iter++) {
    // Get photon tracks...? Also check if these are the same as those from the pipipi0
    EvtRecTrack *HighEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Pi0_iter)->hiEnGamma());
    int HighEPhotonTrackID = HighEnergyPhotonTrack->trackId();
    EvtRecTrack *LowEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Pi0_iter)->loEnGamma());
    int LowEPhotonTrackID = LowEnergyPhotonTrack->trackId();
    if(HighEPhotonTrackID == m_FindPi0.GetHighEPhotonTrackID() && LowEPhotonTrackID == m_FindPi0.GetLowEPhotonTrackID()) {
      continue;
    }
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

std::vector<int> FindKKpipiVersuspipipi0PartRecoTagInfo::GetDaughterTrackID_KKpipi() const {
  return m_DaughterTrackID_KKpipi;
}

int FindKKpipiVersuspipipi0PartRecoTagInfo::GetKSFitSuccess_pipipi0() const {
  return m_KSFitSuccess_pipipi0;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetDecayLengthVeeVertex_pipipi0() const {
  return m_DecayLengthVeeVertex_pipipi0;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetChi2VeeVertex_pipipi0() const {
  return m_Chi2VeeVertex_pipipi0;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetKSMassVeeVertex_pipipi0() const {
  return m_KSMassVeeVertex_pipipi0;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetDecayLengthFit_pipipi0() const {
  return m_DecayLengthFit_pipipi0;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetDecayLengthErrorFit_pipipi0() const {
  return m_DecayLengthErrorFit_pipipi0;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetChi2Fit_pipipi0() const {
  return m_Chi2Fit_pipipi0;
}

int FindKKpipiVersuspipipi0PartRecoTagInfo::GetKalmanFitSuccess() const {
  return m_KalmanFitSuccess;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetKalmanFitChi2() const {
  return m_KalmanFitChi2;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetKPlusP(int i) const {
  return m_KPlusP[i];
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetKMinusP(int i) const {
  return m_KMinusP[i];
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetPiPlusP(int i) const {
  return m_PiPlusP[i];
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetPiMinusP(int i) const {
  return m_PiMinusP[i];
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetKPlusPKalmanFit(int i) const {
  return m_KPlusPKalmanFit[i];
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetKMinusPKalmanFit(int i) const {
  return m_KMinusPKalmanFit[i];
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetPiPlusPKalmanFit(int i) const {
  return m_PiPlusPKalmanFit[i];
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetPiMinusPKalmanFit(int i) const {
  return m_PiMinusPKalmanFit[i];
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetMpipi_KKpipi() const {
  return (m_PiPlusP + m_PiMinusP).m();
}

int FindKKpipiVersuspipipi0PartRecoTagInfo::GetKSFitSuccess_KKpipi() const {
  return m_KSFitSuccess_KKpipi;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetDecayLengthVeeVertex_KKpipi() const {
  return m_DecayLengthVeeVertex_KKpipi;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetChi2VeeVertex_KKpipi() const {
  return m_Chi2VeeVertex_KKpipi;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetKSMassVeeVertex_KKpipi() const {
  return m_KSMassVeeVertex_KKpipi;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetDecayLengthFit_KKpipi() const {
  return m_DecayLengthFit_KKpipi;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetDecayLengthErrorFit_KKpipi() const {
  return m_DecayLengthErrorFit_KKpipi;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetChi2Fit_KKpipi() const {
  return m_Chi2Fit_KKpipi;
}

double FindKKpipiVersuspipipi0PartRecoTagInfo::GetMissingMass2() const {
  return m_MissingMass2;
}

int FindKKpipiVersuspipipi0PartRecoTagInfo::GetRecKCharge() const {
  return m_RecKCharge;
}

int FindKKpipiVersuspipipi0PartRecoTagInfo::GetNumberPi0() const {
  return m_NumberPi0;
}

const FindPi0Eta& FindKKpipiVersuspipipi0PartRecoTagInfo::GetPi0Info() const {
  return m_FindPi0;
}

const FindhhTagInfo& FindKKpipiVersuspipipi0PartRecoTagInfo::GetpipiTagInfo() const {
  return m_FindpipiTagInfo;
}
