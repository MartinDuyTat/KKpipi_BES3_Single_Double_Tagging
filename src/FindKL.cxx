// Martin Duy Tat 25th March 2021, based on code by Yu Zhang

// KKpipi file
#include "KKpipi/FindKL.h"
#include "KKpipi/ParticleMasses.h"
#include "KKpipi/KKpipiUtilities.h"
// Gaudi
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/StatusCode.h"
// Event information
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "ExtEvent/RecExtTrack.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
// STL
#include <cmath>
#include <vector>
// ROOT
#include "TMath.h"

FindKL::FindKL(): m_FoundPionPair(0), m_FoundKaonPair(0), m_KalmanFitSuccess(false), m_NumberPi0(0), m_NumberEta(0), m_NumberGamma(0), m_DaughterTrackID(std::vector<int>(2)), m_KalmanTracks(std::vector<RecMdcKalTrack*>(2)) {
}

FindKL::~FindKL() {
}

StatusCode FindKL::findKL(DTagToolIterator DTTool_iter, DTagTool DTTool) {
  // Prepare message service
  IMessageSvc *msgSvc;
  Gaudi::svcLocator()->service("MessageSvc", msgSvc);
  MsgStream log(msgSvc, "FindKL");
  // Prepare event data service
  IDataProviderSvc *EventDataService = nullptr;
  Gaudi::svcLocator()->service("EventDataSvc", EventDataService);
  // Prepare pi0 service
  SmartDataPtr<EvtRecPi0Col> evtRecPi0Col(EventDataService, "/Event/EvtRec/EvtRecPi0Col");
  if(!evtRecPi0Col) {
    log << MSG::ERROR << "EvtRecPi0Col not found" << endreq;
  }
  // Prepare eta service
  SmartDataPtr<EvtRecEtaToGGCol> evtRecEtaToGGCol(EventDataService, "/Event/EvtRec/EvtRecEtaToGGCol");
  if(!evtRecEtaToGGCol) {
    log << MSG::ERROR << "EvtRecEtaToGGCol not found" << endreq;
  }
  // Get tracks on the other side of the reconstructed D meson
  SmartRefVector<EvtRecTrack> OtherTracks = (*DTTool_iter)->otherTracks();
  // Loop over all tracks on the other side to find pi+ pi- or K+ K-
  int NumberPiPlusTracks = 0, NumberPiMinusTracks = 0, NumberKPlusTracks = 0, NumberKMinusTracks = 0;
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = OtherTracks.begin(); Track_iter != OtherTracks.end(); Track_iter++) {
    // First check if track is valid
    if(!(*Track_iter)->isMdcTrackValid() || !(*Track_iter)->isMdcKalTrackValid()) {
      continue;
    }
    // Look for good pion tracks
    if(DTTool.isPion(*Track_iter)) {
      RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
      MDCKalTrack->setPidType(RecMdcKalTrack::pion);
      if(MDCKalTrack->charge() == +1) {
	NumberPiPlusTracks++;
	m_hPlusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID[0] = (*Track_iter)->trackId();
	KalmanTracks[0] = MDCKalTrack;
      } else if(MDCKalTrack->charge() == -1) {
	NumberPiMinusTracks++;
	m_hMinusP = MDCKalTrack->p4(MASS::PI_MASS);
	m_DaughterTrackID[1] = (*Track_iter)->trackId();
	KalmanTracks[1] = MDCKalTrack;
      } else {
	return StatusCode::FAILURE;
      }
    // Or look for good Kaon tracks
    } else if(DTTool.isKaon(*Track_iter)) {
      RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
      MDCKalTrack->setPidType(RecMdcKalTrack::kaon);
      if(MDCKalTrack->charge() == +1) {
	NumberKPlusTracks++;
	m_hPlusP = MDCKalTrack->p4(MASS::K_MASS);
	m_DaughterTrackID[0] = (*Track_iter)->trackId();
	KalmanTracks[0] = MDCKalTrack;
      } else if(MDCKalTrack->charge() == -1) {
	NumberKMinusTracks++;
	m_hMinusP = MDCKalTrack->p4(MASS::K_MASS);
	m_DaughterTrackID[1] = (*Track_iter)->trackId();
	KalmanTracks[1] = MDCKalTrack;
      } else {
	return StatusCode::FAILURE;
      }
    // If tracks are bad, we don't want the event, could be badly reconstructed KS instead of KL
    } else {
      return StatusCode::FAILURE;
    }
  }
  // If total number of good tracks is not 0 or 2, reject event
  int TotalGoodTracks = NumberPiPlusTracks + NumberPiMinusTracks + NumberKPlusTracks + NumberKMinusTracks;
  if(TotalGoodTracks != 0 && TotalGoodTracks != 2) {
    return StatusCode::FAILURE;
  }
  if(TotalGoodTracks == 0) {
    // If no good tracks are found
    m_FoundPionPair = false;
    m_FoundKaonPair = false;
  } else {
    // If two good tracks are found, it's either a pion pair or a kaon pair
    if(NumberPiPlusTracks == 1 && NumberPiMinusTracks == 1) {
      // Pion pair found
      m_FoundPionPair = true;
      m_FoundKaonPair = false;
    } else if(NumberKPlusTracks == 1 && Number KMinusTracks == 1) {
      // Kaon pair found
      m_FoundPionPair = false
      m_FoundPionPair = true;
    } else {
      // Else we don't want this event
      return StatusCode::FAILURE;
    }
  }
  // Look for pi0
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
    if(Mgammagamma < 0.115 || Mgammagamma > 0.150) {
      continue;
    }
    // Get photon track ID
    int HighEnergyPhotonTrackID = HighEnergyPhotonTrack->trackId();
    int LowEnergyPhotonTrackID = LowEnergyPhotonTrack->trackId();
    // Get photon momenta
    m_Pi0HighEPhotonP.push_back(KKpipiUtilities::GetPhoton4Vector(HighEPhotonShower->energy(), HighEPhotonShower->theta(), HighEPhotonShower->phi()));
    m_Pi0LowEPhotonP.push_back(KKpipiUtilities::GetPhoton4Vector(LowEPhotonShower->energy(), LowEPhotonShower->theta(), LowEPhotonShower->phi()));
    // Get kinematically constrained four-momenta of photons
    m_Pi0HighEPhotonPConstrained.push_back((*Pi0_iter)->hiPfit());
    m_Pi0LowEPhotonPConstrained.push_back((*Pi0_iter)->loPfit());
    m_Pi0Chi2Fit.push_back((*Pi0_iter)->chisq());
    m_Pi0HighEPhotonTrackID.push_back(HighEnergyPhotonTrack->trackId());
    m_Pi0LowEPhotonTrackID.push_back(LowEnergyPhotonTrack->trackId());
    m_NumberPi0++;
  }
  // Look for eta
  for(EvtRecEtaToGGCol::iterator Eta_iter = evtRecEtaToGGCol->begin(); Eta_iter != evtRecEtaToGGCol->end(); Eta_iter++) {
    // Get photon tracks...?
    EvtRecTrack *HighEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Eta_iter)->hiEnGamma());
    EvtRecTrack *LowEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Eta_iter)->loEnGamma());
    // Get EM shower four-momenta of photons
    RecEmcShower *HighEPhotonShower = HighEnergyPhotonTrack->emcShower();
    RecEmcShower *LowEPhotonShower = LowEnergyPhotonTrack->emcShower();
    // Eta invariant mass
    CLHEP::HepLorentzVector HighEPhotonP = KKpipiUtilities::GetPhoton4Vector(HighEPhotonShower->energy(), HighEPhotonShower->theta(), HighEPhotonShower->phi());
    CLHEP::HepLorentzVector LowEPhotonP = KKpipiUtilities::GetPhoton4Vector(LowEPhotonShower->energy(), LowEPhotonShower->theta(), LowEPhotonShower->phi());
    double Mgammagamma = (HighEPhotonP + LowEPhotonP).m();
    if(Mgammagamma < 0.480 || Mgammagamma > 0.580) {
      continue;
    }
    // Get photon track ID
    int HighEnergyPhotonTrackID = HighEnergyPhotonTrack->trackId();
    int LowEnergyPhotonTrackID = LowEnergyPhotonTrack->trackId();
    // Get photon momenta
    m_EtaHighEPhotonP.push_back(KKpipiUtilities::GetPhoton4Vector(HighEPhotonShower->energy(), HighEPhotonShower->theta(), HighEPhotonShower->phi()));
    m_EtaLowEPhotonP.push_back(KKpipiUtilities::GetPhoton4Vector(LowEPhotonShower->energy(), LowEPhotonShower->theta(), LowEPhotonShower->phi()));
    // Get kinematically constrained four-momenta of photons
    m_EtaHighEPhotonPConstrained.push_back((*Eta_iter)->hiPfit());
    m_EtaLowEPhotonPConstrained.push_back((*Eta_iter)->loPfit());
    m_EtaChi2Fit.push_back((*Eta_iter)->chisq());
    m_EtaHighEPhotonTrackID.push_back(HighEnergyPhotonTrack->trackId());
    m_EtaLowEPhotonTrackID.push_back(LowEnergyPhotonTrack->trackId());
    m_NumberEta++;
  }
  // Get showers on the other side of the reconstructed D meson
  SmartRefVector<EvtRecTrack> OtherShowers = (*DTTool_iter)->otherShowers();
  // Loop over all showers to find photons
  for(SmartRefVector<EvtRecTrack>::iterator Shower_iter = OtherShowers.begin(); Shower_iter != OtherShowers.end(); Shower_iter++) {
    // Check if shower is valid
    if(!(*Shower_iter)->isEmcShowerValid()) {
      continue;
    }
    // Get reconstructed EMC shower
    RecEmcShower *EMCShower = (*Shower_iter)->emcShower();
    if(EMCShower->module() == 1 && EMCShower->energy() < 0.025) {
      // Shower in the barrel must have energy larger than 25 MeV
      continue;
    } else if(EMCShower->module() != 1 && EMCShower->energy() < 0.050) {
      // Shower in endcap must have energy larger than 50 MeV
      continue;
    }
    if(EMCShower->time() < 0 || EMCShower->time() > 14) {
      // EMC shower time requirement 0 <= T <= 14 (in units of 50 ns)
      continue;
    }
    // Get EMC position of shower
    CLHEP::Hep3Vector EMCPosition(EMCShower->x(), EMCShower->y(), EMCShower->z());
    // Initialize angles to their maximum
    double Theta;
    double Phi;
    double Angle;
    if(!KKpipiUtilities::GetPhotonAngularSeparation(EMCPosition, Angle, Theta, Phi)) {
      continue;
    }
    m_PhotonP.push_back(KKpipiUtilities::GetPhoton4Vector(EMCShower->energy(), EMCShower->theta(), EMCShower->phi()));
    m_PhotonAngleSeparation.push_back(Angle);
    m_PhotonThetaSeparation.push_back(Theta);
    m_PhotonPhiSeparation.push_back(Phi);
    m_PhotonTrackID.push_back((*Shower_iter)->trackId());
    m_NumberGamma++;
  }
  GetMissingFourMomentum();
  if(m_FoundPionPair && m_NumberPi0 == 0 && m_NumberEta == 0) {
    DoKalmanKinematicFit();
  }
  return StatusCode::SUCCESS;
}

void FindKL::GetMissingFourMomentum() {
  CLHEP::HepLorentzVector P_X;
  if(m_FoundPionPair || m_FoundKaonPair) {
    P_X += m_hPlusP + m_hMinusP;
  }
  for(int i = 0; i < m_NumberPi0; i++) {
    P_X += m_Pi0HighEPhotonPConstrained[i] + m_Pi0LowEPhotonPConstrained[i];
  }
  m_KLongP = KKpipiUtilities::GetMissingMomentum((*DTTool_iter)->p4(), P_X, (*DTTool_iter)->beamE());
}

void FindKL::GetMissingMass2() {
  return m_KLongP.m2();
}

void FindKL::DoKalmanKinematicFit() {
  double hMass;
  if(m_FoundPionPair && !m_FoundKaonPair) {
    hMass = MASS::PI_MASS;
  } else if(!m_FoundPionPair && m_FoundKaonPair) {
    hMass = MASS::K_MASS;
  }
  WTrackParameter WTrackKplus(hMass, m_KalmanTracks[0]->getZHelix(), m_KalmanTracks[0]->getZError());
  WTrackParameter WTrackKminus(hMass, m_KalmanTracks[1]->getZHelix(), m_KalmanTracks[1]->getZError());
  WTrackParameter WTrackKLong(m_KLongP, m_KLongP.phi(), m_KLongP.theta(), m_KLongP.t());
  WTrackKLong.setCharge(0);
  WTrackKLong.setMass(MASS::KS_MASS);
  KalmanKinematicFit *KalmanFit = KalmanKinematicFit::instance();
  KalmanFit->init();
  KalmanFit->AddTrack(0, WTrackKplus);
  KalmanFit->AddTrack(1, WTrackKminus);
  KalmanFit->AddTrack(2, WTrackKLong);
  KalmanFit->AddResonance(0, MASS::D_MASS, 0, 1, 2);
  m_KalmanFitSuccess = KalmanFit->Fit();
  if(m_KalmanFitSuccess) {
    m_hPlusPKalmanFit = KalmanFit->pfit(0);
    m_hMinusPKalmanFit = KalmanFit->pfit(1);
    m_KLongPKalmanFit = KalmanFit->pfit(2);
  }
}

bool FindKL::FoundKLpipiTag() const {
  if(!m_FoundPionPair || m_FoundKaonPair || m_NumberPi0 != 0 || m_NumberEta != 0) {
    return false;
  } else {
    return true;
  }
}

bool FindKL::FoundKLKKTag() const {
  if(m_FoundPionPair || !m_FoundKaonPair || m_NumberPi0 != 0 || m_NumberEta != 0) {
    return false;
  } else {
    return true;
  }
}

bool FindKL::FoundKLpi0Tag() {
  if(m_FoundPionPair || m_FoundKaonPair || m_NumberEta != 0 || m_NumberPi0 != 1) {
    return false;
  } else {
    return true;
  }
}

bool FindKL::FoundKLpi0pi0Tag() {
  if(m_FoundPionPair || m_FoundKaonPair || m_NumberEta != 0 || m_NumberPi0 != 2) {
    return false;
  } else {
    // Make sure the two photons don't share common showers
    if(m_Pi0HighEPhotonTrackID[0] == m_Pi0HighEPhotonTrackID[1] ||
       m_Pi0HighEPhotonTrackID[0] == m_Pi0LowEPhotonTrackID[1] ||
       m_Pi0LowEPhotonTrackID[0] == m_Pi0HighEPhotonTrackID[1] ||
       m_Pi0LowEPhotonTrackID[0] == m_Pi0LowEPhotonTrackID[1] ||) {
      return false;
    } else {
      return true;
    }
  }
}

bool FindKL::FoundKLpipipi0Tag() {
  if(!m_FoundPionPair || m_FoundKaonPair || m_NumberEta != 0 || m_NumberPi0 != 1) {
    return false;
  } else {
    return true;
  }
}
  

bool FindKL::GetFoundPionPair() const {
  return m_FoundPionPair;
}

bool FindKL::GetFoundKaonPair() const {
  return m_FoundKaonPair;
}

double FindKL::GetPiPlusP(int i) const {
  return m_PiPlusP[i];
}

double FindKL::GetPiMinusP(int i) const {
  return m_PiMinusP[i];
}

double FindKL::GetPi0HighEPhotonP(int i, int j) const{
  return m_Pi0HighEPhotonP[j][i];
}

double FindKL::GetPi0LowEPhotonP(int i, int j) const{
  return m_Pi0LowEPhotonP[j][i];
}

double FindKL::GetPi0Mgammagamma(int j) const {
  return (m_Pi0HighEPhotonP[j] + m_Pi0LowEPhotonP[j]).m();
}

double FindKL::GetPi0HighEPhotonPConstrained(int i, int j) const{
  return m_Pi0HighEPhotonPConstrained[j][i];
}

double FindKL::GetPi0LowEPhotonPConstrained(int i, int j) const{
  return m_Pi0LowEPhotonPConstrained[j][i];
}

double FindKL::GetPi0Chi2Fit(int j) const {
  return m_Pi0Chi2Fit[j];
}

int FindKL::GetPi0HighEPhotonTrackID(int j) const{
  return m_Pi0HighEPhotonTrackID[j];
}

int FindKL::GetPi0LowEPhotonTrackID(int j) const{
  return m_Pi0LowEPhotonTrackID[j];
}

int FindKL::GetNumberPi0() const {
  return m_NumberPi0;
}

double FindKL::GetEtaHighEPhotonP(int i, int j) const{
  return m_EtaHighEPhotonP[j][i];
}

double FindKL::GetEtaLowEPhotonP(int i, int j) const{
  return m_EtaLowEPhotonP[j][i];
}

double FindKL::GetEtaMgammagamma(int j) const {
  return (m_EtaHighEPhotonP[j] + m_EtaLowEPhotonP[j]).m();
}

double FindKL::GetEtaHighEPhotonPConstrained(int i, int j) const{
  return m_EtaHighEPhotonPConstrained[j][i];
}

double FindKL::GetEtaLowEPhotonPConstrained(int i, int j) const{
  return m_EtaLowEPhotonPConstrained[j][i];
}

double FindKL::GetEtaChi2Fit(int j) const {
  return m_EtaChi2Fit[j];
}

int FindKL::GetEtaHighEPhotonTrackID(int j) const{
  return m_EtaHighEPhotonTrackID[j];
}

int FindKL::GetEtaLowEPhotonTrackID(int j) const{
  return m_EtaLowEPhotonTrackID[j];
}

int FindKL::GetNumberEta() const {
  return m_NumberEta;
}

double FindKL::GetPhotonP(int i, int j) const {
  return m_PhotonP[j][i];
}

double FindKL::GetPhotonAngleSeparation(int j) const {
  return m_PhotonAngleSeparation[j];
}

double FindKL::GetPhotonThetaSeparation(int j) const {
  return m_PhotonThetaSeparation[j];
}

double FindKL::GetPhotonPhiSeparation(int j) const {
  return m_PhotonPhiSeparation[j];
}

int FindKL::GetPhotonTrackID(int j) const {
  return m_PhotonTrackID[j];
}

int FindKL::GetNumberGamma() const {
  return m_NumberGamma;
}

std::vector<int> FindKL::GetDaughterTrackID() const {
  return m_DaughterTrackID;
}
