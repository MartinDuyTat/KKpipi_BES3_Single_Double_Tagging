// Martin Duy Tat 25th March 2021, based on code by Yu Zhang

// Header file
#include "KKpipi/FindKL.h"
#include "KKpipi/ParticleMasses.h"
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
#include "EvtRecEvent/EvtRecVeeVertex.h"
// CLHEP
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector/LorentzVector.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/VertexParameter.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/WTrackParameter.h"
// STL
#include<algorithm>
#include<vector>
// Particle masses
#include "KKpipi/ParticleMasses.h"

FindKL::FindKL(): FoundPionPair(0) {
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
  SmartDataPtr<EvtRecTracCol> evtRecTrkCol(EventDataService, "/Event/EvtRec/EvtRecTrackCol");
  if(!evtRecTrkCol) {
    log << MSG::ERROR << "EvtRecTrackCol not found" << endreq;
  }
  // Prepare pi0 service
  SmartDataPtr<EvtRecPi0Col> evtRecPi0Col(EventDataService, "/Event/EvtRec/EvtRecPi0Col");
  if(!evtRecPi0Col) {
    log << MSG::ERROR << "EvtRecPi0Col not found" << endreq;
  }
  // Prepare eta service
  SmartDataPtr<EvtRecEtaToGGCol> evtRecEtaToGG(EventDataService, "/Event/EvtRec/EvtRecEtaToGGCol");
  if(!evtRecEtaToGG) {
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
    if(DTTool.isPion(*Track_iter)) {
      RecMdcKalTrack *MDCKalTrack = (*Track_iter)->mdcKalTrack();
      MDCKalTrack->setPidType(RecMdcKalTrack::pion);
      if(MDCKalTrack->charge = +1) {
	NumberPiPlusTracks++;
	m_PiPlusP = MDCKalTrack->p4(MASS::PI_MASS);
      } else if(MDCKalTrack->charge = -1) {
	NumberPiMinusTracks++;
	m_PiMinusP = MDCKalTrack->p4(MASS::PI_MASS);
      }
    }
  }
  // If no pions are found, or a pi+ pi- pair is found, keep going, otherwise reject event
  if(NumberPiPlusTracks == NumberPiMinusTracks) {
    if(NumberPiPlusTracks == 1) {
      FoundPionPair == 1;
    } else if (NumberPiPlusTracks != 0) {
      return StatusCode::FAILURE;
    }
  } else {
    return Statuscode::FAILURE;
  }
}
