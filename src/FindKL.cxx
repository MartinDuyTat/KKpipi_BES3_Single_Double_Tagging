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

FindKL::FindKL(): m_FoundPionPair(0), m_NumberPi0(0), m_NumberEta(0), m_NumberGamma(0), m_DaughterTrackID(std::vector<int>(2)) {
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
  // Prepare reconstructed event service
  SmartDataPtr<EvtRecEvent> evtRecEvent(EventDataService, EventModel::EvtRec::EvtRecEvent);
  if(!evtRecEvent) {
    log << MSG::ERROR << "EvtRecEvent not found" << endreq;
  }
  // Prepare event tracks service
  SmartDataPtr<EvtRecTrackCol> evtRecTrackCol(EventDataService, "/Event/EvtRec/EvtRecTrackCol");
  if(!evtRecTrackCol) {
    log << MSG::ERROR << "EvtRecTrackCol not found" << endreq;
  }
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
  // Loop over all tracks on the other side to find pi+ pi-
  int NumberPiPlusTracks = 0, NumberPiMinusTracks = 0;
  for(SmartRefVector<EvtRecTrack>::iterator Track_iter = OtherTracks.begin(); Track_iter != OtherTracks.end(); Track_iter++) {
    // First check if track is valid
    if(!(*Track_iter)->isMdcTrackValid() || !(*Track_iter)->isMdcKalTrackValid()) {
      continue;
    }
    if(DTTool.isGoodTrack(*Track_iter)) {
      if(DTTool.isPion(*Track_iter)) {
	RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
	MDCKalTrack->setPidType(RecMdcKalTrack::pion);
	if(MDCKalTrack->charge() == +1) {
	  NumberPiPlusTracks++;
	  m_PiPlusP = MDCKalTrack->p4(MASS::PI_MASS);
	  m_DaughterTrackID[0] = (*Track_iter)->trackId();
	} else if(MDCKalTrack->charge() == -1) {
	  NumberPiMinusTracks++;
	  m_PiMinusP = MDCKalTrack->p4(MASS::PI_MASS);
	  m_DaughterTrackID[1] = (*Track_iter)->trackId();
	} else {
	  return StatusCode::FAILURE;
	}
      } else {
	return StatusCode::FAILURE;
      }
    }
  }
  // If no pions are found, or a pi+ pi- pair is found, keep going, otherwise reject event
  if(NumberPiPlusTracks == NumberPiMinusTracks) {
    if(NumberPiPlusTracks == 1) {
      m_FoundPionPair = 1;
    } else if (NumberPiPlusTracks != 0) {
      return StatusCode::FAILURE;
    }
  } else {
    return StatusCode::FAILURE;
  }
  // Check if there are any showers already present on the reconstructed side
  bool ShowersUsed;
  SmartRefVector<EvtRecTrack> Showers = (*DTTool_iter)->showers();
  if(Showers.size() == 0) {
    // If no showers, it's probably a Kpi, KKpipi or Kpipipi tag on the reconstructed side
    ShowersUsed = false;
  } else if(Showers.size() > 2) {
    // There shouldn't be more than 2 showers or just 1 shower
    return StatusCode::FAILURE;
  } else if(Showers.size() == 1) {
    log << MSG::ERROR << "Found a single shower when looking for KL" << endreq;
    return StatusCode::FAILURE;
  } else {
    // If 2 showers, it's probably a Kpipi0 tag on the reconstructed side
    ShowersUsed = true;
  }
  // Look for pi0
  for(EvtRecPi0Col::iterator Pi0_iter = evtRecPi0Col->begin(); Pi0_iter != evtRecPi0Col->end(); Pi0_iter++) {
    // Get photon tracks...?
    EvtRecTrack *HighEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Pi0_iter)->hiEnGamma());
    EvtRecTrack *LowEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Pi0_iter)->loEnGamma());
    // Get EM shower four-momenta of photons
    RecEmcShower *HighEPhotonShower = HighEnergyPhotonTrack->emcShower();
    RecEmcShower *LowEPhotonShower = LowEnergyPhotonTrack->emcShower();
    // Get photon track ID
    int HighEnergyPhotonTrackID = HighEnergyPhotonTrack->trackId();
    int LowEnergyPhotonTrackID = LowEnergyPhotonTrack->trackId();
    // If reconstructed side constains a shower, make sure this is not the same shower
    if(ShowersUsed) {
      if(HighEnergyPhotonTrackID == Showers[0]->trackId() || HighEnergyPhotonTrackID == Showers[1]->trackId()) {
	continue;
      }
      if(LowEnergyPhotonTrackID == Showers[0]->trackId() || LowEnergyPhotonTrackID == Showers[1]->trackId()) {
	continue;
      }
    }
    // Get photon momenta
    m_Pi0HighEPhotonP.push_back(KKpipiUtilities::GetPhoton4Vector(HighEPhotonShower->energy(), HighEPhotonShower->theta(), HighEPhotonShower->phi()));
    m_Pi0LowEPhotonP.push_back(KKpipiUtilities::GetPhoton4Vector(LowEPhotonShower->energy(), LowEPhotonShower->theta(), LowEPhotonShower->phi()));
    // Get kinematically constrained four-momenta of photons
    m_Pi0HighEPhotonPConstrained.push_back((*Pi0_iter)->hiPfit());
    m_Pi0LowEPhotonPConstrained.push_back((*Pi0_iter)->loPfit());
    m_Pi0Chi2Fit.push_back((*Pi0_iter)->chisq());
    m_NumberPi0++;
  }
  // If there are no pi0, reject event
  if(m_NumberPi0 == 0) {
    return StatusCode::FAILURE;
  }
  // Look for eta
  for(EvtRecEtaToGGCol::iterator Eta_iter = evtRecEtaToGGCol->begin(); Eta_iter != evtRecEtaToGGCol->end(); Eta_iter++) {
    // Get photon tracks...?
    EvtRecTrack *HighEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Eta_iter)->hiEnGamma());
    EvtRecTrack *LowEnergyPhotonTrack = const_cast<EvtRecTrack*>((*Eta_iter)->loEnGamma());
    // Get EM shower four-momenta of photons
    RecEmcShower *HighEPhotonShower = HighEnergyPhotonTrack->emcShower();
    RecEmcShower *LowEPhotonShower = LowEnergyPhotonTrack->emcShower();
    // Get photon track ID
    int HighEnergyPhotonTrackID = HighEnergyPhotonTrack->trackId();
    int LowEnergyPhotonTrackID = LowEnergyPhotonTrack->trackId();
    // If reconstructed side constains a shower, make sure this is not the same shower
    if(ShowersUsed) {
      if(HighEnergyPhotonTrackID == Showers[0]->trackId() || HighEnergyPhotonTrackID == Showers[1]->trackId()) {
	continue;
      }
      if(LowEnergyPhotonTrackID == Showers[0]->trackId() || LowEnergyPhotonTrackID == Showers[1]->trackId()) {
	continue;
      }
    }
    // Get photon momenta
    m_EtaHighEPhotonP.push_back(KKpipiUtilities::GetPhoton4Vector(HighEPhotonShower->energy(), HighEPhotonShower->theta(), HighEPhotonShower->phi()));
    m_EtaLowEPhotonP.push_back(KKpipiUtilities::GetPhoton4Vector(LowEPhotonShower->energy(), LowEPhotonShower->theta(), LowEPhotonShower->phi()));
    // Get kinematically constrained four-momenta of photons
    m_EtaHighEPhotonPConstrained.push_back((*Eta_iter)->hiPfit());
    m_EtaLowEPhotonPConstrained.push_back((*Eta_iter)->loPfit());
    m_EtaChi2Fit.push_back((*Eta_iter)->chisq());
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
    double Theta = 2*TMath::Pi();
    double Phi = 2*TMath::Pi();
    double Angle = 2*TMath::Pi();
    // Loop over all charged tracks
    for(int j = 0; j < evtRecEvent->totalCharged(); j++) {
      EvtRecTrackIterator Track_iter = evtRecTrackCol->begin() + j;
      // Check if track is valid
      if(!(*Track_iter)->isExtTrackValid()) {
	continue;
      }
      // Get track extrapolated to the EMC
      RecExtTrack *ExternalTrack = (*Track_iter)->extTrack();
      // Check if external track is valid
      if(ExternalTrack->emcVolumeNumber() == -1) {
	continue;
      }
      // Get position of external track
      CLHEP::Hep3Vector ExternalPosition = ExternalTrack->emcPosition();
      // Find angle between track and shower
      double DeltaAngle = ExternalPosition.angle(EMCPosition);
      // Find polar angle between track and shower
      double DeltaTheta = ExternalPosition.theta() - EMCPosition.theta();
      // Find azimuthal angle between track and shower
      double DeltaPhi = ExternalPosition.deltaPhi(EMCPosition);
      if(DeltaAngle < Angle) {
	Theta = DeltaTheta;
	Phi = DeltaPhi;
	Angle = DeltaAngle;
      }
    }
    if(Angle == 2*TMath::Pi()) {
      log << MSG::ERROR << "No charged tracks found when looking for photons";
    }
    m_PhotonEnergy.push_back(EMCShower->energy());
    m_PhotonAngleSeparation.push_back(Angle);
    m_PhotonThetaSeparation.push_back(Theta);
    m_PhotonPhiSeparation.push_back(Phi);
    m_NumberGamma++;
  }
  return StatusCode::SUCCESS;
}

int FindKL::GetFoundPionPair() const {
  return m_FoundPionPair;
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

double FindKL::GetPi0HighEPhotonPConstrained(int i, int j) const{
  return m_Pi0HighEPhotonPConstrained[j][i];
}

double FindKL::GetPi0LowEPhotonPConstrained(int i, int j) const{
  return m_Pi0LowEPhotonPConstrained[j][i];
}

double FindKL::GetPi0Chi2Fit(int j) const {
  return m_Pi0Chi2Fit[j];
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

double FindKL::GetEtaHighEPhotonPConstrained(int i, int j) const{
  return m_EtaHighEPhotonPConstrained[j][i];
}

double FindKL::GetEtaLowEPhotonPConstrained(int i, int j) const{
  return m_EtaLowEPhotonPConstrained[j][i];
}

double FindKL::GetEtaChi2Fit(int j) const {
  return m_EtaChi2Fit[j];
}

int FindKL::GetNumberEta() const {
  return m_NumberEta;
}

double FindKL::GetPhotonEnergy(int j) const {
  return m_PhotonEnergy[j];
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

int FindKL::GetNumberGamma() const {
  return m_NumberGamma;
}

std::vector<int> FindKL::GetDaughterTrackID() const {
  return m_DaughterTrackID;
}
