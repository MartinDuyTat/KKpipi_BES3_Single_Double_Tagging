// Martin Duy Tat 5th February 2021, based on code by Yu Zhang

// Header file
#include "KKpipi/FindKS.h"
// Gaudi
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"
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
// Boss
#include "DTagTool/DTagTool.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/VertexParameter.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/WTrackParameter.h"
// STL
#include<vector>
// Particle masses
#include "KKpipi/ParticleMasses.h"

#include<iostream>// remove

FindKS::FindKS(): m_DecayLengthVeeVertex(0.0), m_Chi2VeeVertex(0.0), m_KSMassVeeVertex(0.0), m_DecayLengthFit(0.0), m_DecayLengthErrorFit(0.0), m_Chi2PrimaryVertexFit(0.0), m_Chi2SecondaryVertexFit(0.0), m_KSMassFit(0.0) {
}

FindKS::~FindKS() {
}

StatusCode FindKS::findKS(DTagToolIterator &DTTool_iter, const std::vector<int> &PiTrackIndex) {
  IMessageSvc *msgSvc;
  Gaudi::svcLocator()->service("MessageSvc", msgSvc);
  MsgStream log(msgSvc, "FindKS");
  // Check if event has two pions
  if(PiTrackIndex.size() != 2) {
    log << MSG::ERROR << "Need two pions to reconstruct KS" << endreq;
  }
  IDataProviderSvc *eventSvc = nullptr;
  Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
  SmartDataPtr<EvtRecVeeVertexCol> evtRecVeeVertexCol(eventSvc, "/Event/EvtRec/EvtRecVeeVertexCol");
  if(!evtRecVeeVertexCol) {
    log << MSG::ERROR << "EvtRecVeeVertexCol not found" << endreq;
  }
  // Get tracks in the event
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  // Get Kalman tracks and pion track IDs
  RecMdcKalTrack *MDCKalmanTrack1 = (*Tracks.begin() + PiTrackIndex[0])->mdcKalTrack();
  int PiTrackID1 = (*Tracks.begin() + PiTrackIndex[0])->trackId();
  RecMdcKalTrack *MDCKalmanTrack2 = (*Tracks.begin() + PiTrackIndex[1])->mdcKalTrack();
  int PiTrackID2 = (*Tracks.begin() + PiTrackIndex[1])->trackId();
  // Loop over KS in the event (should only be one)
  for(EvtRecVeeVertexCol::iterator KS_iter = evtRecVeeVertexCol->begin(); KS_iter != evtRecVeeVertexCol->end(); KS_iter++) {
    // Get KS daughter tracks
    EvtRecTrack *KSChildTrack1 = (*KS_iter)->daughter(0);
    EvtRecTrack *KSChildTrack2 = (*KS_iter)->daughter(1);
    // Get KS daughter track IDs
    int KSChildTrackID1 = KSChildTrack1->trackId();
    int KSChildTrackID2 = KSChildTrack2->trackId();
    // Check if KS daughter tracks are the same as the pion tracks
    if(!((KSChildTrackID1 == PiTrackID1 && KSChildTrackID2 == PiTrackID2) || (KSChildTrackID1 == PiTrackID2 && KSChildTrackID2 == PiTrackID1))) {
      continue;
    }
    // Get KS position vector from VeeVertexAlg
    CLHEP::Hep3Vector KS_PositionVector((*KS_iter)->w()[4], (*KS_iter)->w()[5], (*KS_iter)->w()[6]);
    // Calculate length of position vector, which is the decay length
    m_DecayLengthVeeVertex = KS_PositionVector.mag();
    // Get VeeVertexAlg chi2 value
    m_Chi2VeeVertex = (*KS_iter)->chi2();
    // Get VeeVertexAlg KS mass
    m_KSMassVeeVertex = (*KS_iter)->mass();
    // Set up initial guess for secondary vertex position and error and put into a VertexParameter object
    HepPoint3D SecondaryVertexPosition(0.0, 0.0, 0.0);
    CLHEP::HepSymMatrix SecondaryVertexError(3, 0);
    std::cout << "One!\n"; // remove
    SecondaryVertexError[0][0] = 1e6;
    SecondaryVertexError[1][1] = 1e6;
    SecondaryVertexError[2][2] = 1e6;
    std::cout << "Two!\n"; // remove
    VertexParameter SecondaryVertexParam;
    SecondaryVertexParam.setVx(SecondaryVertexPosition);
    SecondaryVertexParam.setEvx(SecondaryVertexError);
    // Get Kalman fitted pion tracks and their track parameters
    RecMdcKalTrack *KSChildKalmanTrack1 = KSChildTrack1->mdcKalTrack();
    RecMdcKalTrack *KSChildKalmanTrack2 = KSChildTrack2->mdcKalTrack();
    WTrackParameter WTrackPion1(MASS::PI_MASS, KSChildKalmanTrack1->helix(), KSChildKalmanTrack1->err());
    WTrackParameter WTrackPion2(MASS::PI_MASS, KSChildKalmanTrack2->helix(), KSChildKalmanTrack2->err());
    // Start fitting secondary vertex
    VertexFit *SecondaryVertexFit = VertexFit::instance();
    SecondaryVertexFit->AddTrack(0, WTrackPion1);
    SecondaryVertexFit->AddTrack(1, WTrackPion2);
    SecondaryVertexFit->AddVertex(0, SecondaryVertexParam, 0, 1);
    SecondaryVertexFit->Fit(0);
    SecondaryVertexFit->BuildVirtualParticle(0);
    // Save fitted track parameters of the KS
    WTrackParameter WTrackKS = SecondaryVertexFit->wVirtualTrack(0);
    // Get VertexDbSvc, which determines the average beam position for each run
    IVertexDbSvc *VertexService;
    Gaudi::svcLocator()->service("VertexDbSvc", VertexService);
    if(!VertexService->isVertexValid()) {
      continue;
    }
    // Get primary vertex position and error
    double *PVertex = VertexService->PrimaryVertex();
    double *SigmaPV = VertexService->SigmaPrimaryVertex();
    // Put parameters into a VertexParameter object for fitting
    HepPoint3D PrimaryVertex(PVertex[0], PVertex[1], PVertex[2]);
    CLHEP::HepSymMatrix PVError(3, 0);
    std::cout << "Three!\n"; // remove
    PVError[0][0] = SigmaPV[0]*SigmaPV[0];
    PVError[1][1] = SigmaPV[1]*SigmaPV[1];
    PVError[2][2] = SigmaPV[2]*SigmaPV[2];
    std::cout << "Four!\n"; // remove
    VertexParameter KSOrigin;
    KSOrigin.setVx(PrimaryVertex);
    KSOrigin.setEvx(PVError);
    // Start fitting primary vertex (confusing name)
    SecondVertexFit *PrimaryVertexFit = SecondVertexFit::instance();
    PrimaryVertexFit->init();
    PrimaryVertexFit->AddTrack(0, WTrackKS);
    PrimaryVertexFit->setVpar(SecondaryVertexFit->vpar(0));
    PrimaryVertexFit->setPrimaryVertex(KSOrigin);
    if(PrimaryVertexFit->Fit()) {
      m_DecayLengthFit = PrimaryVertexFit->decayLength();
      m_DecayLengthErrorFit = PrimaryVertexFit->decayLengthError();
      m_Chi2PrimaryVertexFit = PrimaryVertexFit->chisq();
      m_Chi2SecondaryVertexFit = SecondaryVertexFit->chisq();
      m_KSMassFit = PrimaryVertexFit->p4par().m();
    }
    return StatusCode::SUCCESS;
  }
  return StatusCode::FAILURE;
}

double FindKS::getDecayLengthVeeVertex() const {
  return m_DecayLengthVeeVertex;
}

double FindKS::getChi2VeeVertex() const {
  return m_Chi2VeeVertex;
}

double FindKS::getKSMassVeeVertex() const {
  return m_KSMassVeeVertex;
}

double FindKS::getDecayLengthFit() const {
  return m_DecayLengthFit;
}

double FindKS::getDecayLengthErrorFit() const {
  return m_DecayLengthErrorFit;
}

double FindKS::getChi2PrimaryVertexFit() const {
  return m_Chi2PrimaryVertexFit;
}

double FindKS::getChi2SecondaryVertexFit() const {
  return m_Chi2SecondaryVertexFit;
}

double FindKS::getKSMassFit() const {
  return m_KSMassFit;
}
