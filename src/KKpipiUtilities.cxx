// Martin Duy Tat 26th March 2021

// KKpipi
#include "KKpipi/KKpipiUtilities.h"
// Gaudi
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
// Event information
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "ExtEvent/RecExtTrack.h"
// ROOT
#include "TMath.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"

CLHEP::HepLorentzVector KKpipiUtilities::GetPhoton4Vector(double Energy, double Theta, double Phi) {
  double Px = Energy*TMath::Sin(Theta)*TMath::Cos(Phi);
  double Py = Energy*TMath::Sin(Theta)*TMath::Sin(Phi);
  double Pz = Energy*TMath::Cos(Theta);
  return CLHEP::HepLorentzVector(Px, Py, Pz, Energy);
}

bool KKpipiUtilities::GetPhotonAngularSeparation(const CLHEP::Hep3Vector &EMCPosition, double &Angle, double &Theta, double &Phi) {
  Theta = 2*TMath::Pi();
  Phi = 2*TMath::Pi();
  Angle = 2*TMath::Pi();
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
    return false;
  } else {
    return true;
  }
}
