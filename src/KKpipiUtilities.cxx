// Martin Duy Tat 26th March 2021

//STL
#include <cmath>
// KKpipi
#include "KKpipi/KKpipiUtilities.h"
#include "KKpipi/ParticleMasses.h"
// Gaudi
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ISvcLocator.h"
// Event information
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "ExtEvent/RecExtTrack.h"
// ROOT
#include "TMath.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Geometry/Point3D.h"
// Boss
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/Helix.h"

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
  // Prepare message service
  IMessageSvc *msgSvc;
  Gaudi::svcLocator()->service("MessageSvc", msgSvc);
  MsgStream log(msgSvc, "GetPhotonAngularSeparation");
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

CLHEP::HepLorentzVector KKpipiUtilities::GetMissingMomentum(CLHEP::HepLorentzVector P_D, CLHEP::HepLorentzVector P_X, double BeamE) {
  // Need to boost because the beams have a crossing angle of 11 mrad
  CLHEP::Hep3Vector CrossingAngle(-0.011, 0.0, 0.0);
  P_D.boost(CrossingAngle);
  P_X.boost(CrossingAngle);
  // Get the direction unit vector of the D meson three-momentum
  CLHEP::Hep3Vector P_Dunit = P_D.vect().unit();
  // Get the constrained magnitude of the D meson momentum
  double DMomentum = TMath::Sqrt(BeamE*BeamE - MASS::D_MASS*MASS::D_MASS);
  CLHEP::HepLorentzVector P_Dconstrained(DMomentum*P_Dunit, BeamE);
  // Use conservation of four-momentum to find missing four-momentum
  return CLHEP::HepLorentzVector(0.0, 0.0, 0.0, MASS::JPSI_MASS) - P_Dconstrained - P_X;
}

void KKpipiUtilities::GetIP(RecMdcKalTrack *MDCKalTrack, double &IP_Vxy, double &IP_Vz) {
  Hep3Vector xorigin(0,0,0);
  IVertexDbSvc*  vtxsvc;
  Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
  if(vtxsvc->isVertexValid()){
    double* dbv = vtxsvc->PrimaryVertex();
    xorigin.setX(dbv[0]);
    xorigin.setY(dbv[1]);
    xorigin.setZ(dbv[2]);
  }
  HepVector a = MDCKalTrack->getZHelixK();
  HepSymMatrix Ea = MDCKalTrack->getZErrorK();
  HepPoint3D point0(0.,0.,0.);
  HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]);
  VFHelix helixip3(point0,a,Ea);
  helixip3.pivot(IP);
  HepVector  vecipa = helixip3.a();
  IP_Vxy = fabs(vecipa[0]);
  IP_Vz = fabs(vecipa[3]);
}
