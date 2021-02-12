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

FindKS::FindKS(): m_DecayLengthVeeVertex(0.0), m_Chi2VeeVertex(0.0), m_KSMassVeeVertex(0.0), m_DecayLengthFit(0.0), m_DecayLengthErrorFit(0.0), m_Chi2Fit(0.0), m_KSMassFit(0.0) {
}

FindKS::~FindKS() {
}

StatusCode FindKS::findKS(DTagToolIterator &DTTool_iter, const std::vector<SmartRefVector<EvtRecTrack>::iterator> &PiTrack_iter) {
  IMessageSvc *msgSvc;
  Gaudi::svcLocator()->service("MessageSvc", msgSvc);
  MsgStream log(msgSvc, "FindKS");
  // Check if event has two pions
  if(PiTrack_iter.size() != 2 || PiTrack_iter.size() != 0) {
    log << MSG::ERROR << "Need two pions to reconstruct KS" << endreq;
    return StatusCode::FAILURE;
  }
  IDataProviderSvc *eventSvc = nullptr;
  Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
  SmartDataPtr<EvtRecVeeVertexCol> evtRecVeeVertexCol(eventSvc, "/Event/EvtRec/EvtRecVeeVertexCol");
  if(!evtRecVeeVertexCol) {
    log << MSG::ERROR << "EvtRecVeeVertexCol not found" << endreq;
    return StatusCode::FAILURE;
  }
  // Get tracks in the event
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  // Get Kalman tracks and pion track IDs
  int PiTrackID1 = (*PiTrack_iter[0])->trackId();
  int PiTrackID2 = (*PiTrack_iter[1])->trackId();
  // Loop over KS in the event (should only be one)
  for(EvtRecVeeVertexCol::iterator KS_iter = evtRecVeeVertexCol->begin(); KS_iter != evtRecVeeVertexCol->end(); KS_iter++) {
    // Check if the vertex is actually a KS
    if((*KS_iter)->vertexId() != 310) {
      continue;
    }
    // Get KS daughter tracks
    EvtRecTrack *KSChildTrack1 = (*KS_iter)->daughter(0);
    EvtRecTrack *KSChildTrack2 = (*KS_iter)->daughter(1);
    // Get KS daughter track IDs
    int KSChildTrackID1 = KSChildTrack1->trackId();
    int KSChildTrackID2 = KSChildTrack2->trackId();
    // Check if KS daughter tracks are the same as the pion tracks (if pion tracks are given)
    if(PiTrack_iter.size() != 0 && !((KSChildTrackID1 == PiTrackID1 && KSChildTrackID2 == PiTrackID2) || (KSChildTrackID1 == PiTrackID2 && KSChildTrackID2 == PiTrackID1))) {
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
    SecondaryVertexError[0][0] = 100.0;
    SecondaryVertexError[1][1] = 100.0;
    SecondaryVertexError[2][2] = 100.0;
    VertexParameter SecondaryVertexParam;
    SecondaryVertexParam.setVx(SecondaryVertexPosition);
    SecondaryVertexParam.setEvx(SecondaryVertexError);
    // Get Kalman fitted pion tracks and their track parameters
    RecMdcKalTrack *KSChildKalmanTrack1 = KSChildTrack1->mdcKalTrack();
    RecMdcKalTrack *KSChildKalmanTrack2 = KSChildTrack2->mdcKalTrack();
    KSChildKalmanTrack1->setPidType(RecMdcKalTrack::pion);
    KSChildKalmanTrack2->setPidType(RecMdcKalTrack::pion);
    WTrackParameter WTrackPion1(MASS::PI_MASS, KSChildKalmanTrack1->helix(), KSChildKalmanTrack1->err());
    WTrackParameter WTrackPion2(MASS::PI_MASS, KSChildKalmanTrack2->helix(), KSChildKalmanTrack2->err());
    // Store the four-momenta of daughter particles from the MDC tracks
    m_KSPiPlusP = KSChildKalmanTrack1->p4(MASS::PI_MASS);
    m_KSPiMinusP = KSChildKalmanTrack2->p4(MASS::PI_MASS);
    // Start fitting secondary vertex
    VertexFit *SecondaryVertexFit = VertexFit::instance();
    SecondaryVertexFit->init();
    SecondaryVertexFit->AddTrack(0, WTrackPion1);
    SecondaryVertexFit->AddTrack(1, WTrackPion2);
    SecondaryVertexFit->AddVertex(0, SecondaryVertexParam, 0, 1);
    SecondaryVertexFit->Fit(0);
    SecondaryVertexFit->BuildVirtualParticle(0);
    // Save fitted track parameters of the KS
    WTrackParameter WTrackKS = SecondaryVertexFit->wVirtualTrack(0);
    // Store the four-momenta of daughter particles from fit
    m_KSPiPlusPFit = SecondaryVertexFit->wTrackInfit(0).p();
    m_KSPiMinusPFit = SecondaryVertexFit->wTrackInfit(1).p();
    // Swap pions if charges are the other way around
    if(KSChildKalmanTrack1->charge() < 0) {
      swap(m_KSPiPlusP, m_KSPiMinusP);
      swap(m_KSPiPlusPFit, m_KSPiMinusPFit);
    }
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
    PVError[0][0] = SigmaPV[0]*SigmaPV[0];
    PVError[1][1] = SigmaPV[1]*SigmaPV[1];
    PVError[2][2] = SigmaPV[2]*SigmaPV[2];
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
      m_Chi2Fit = PrimaryVertexFit->chisq();
      m_KSMassFit = PrimaryVertexFit->p4par().m();
      return StatusCode::SUCCESS;
    }
  }
  return StatusCode::FAILURE;
}

double FindKS::GetDecayLengthVeeVertex() const {
  return m_DecayLengthVeeVertex;
}

double FindKS::GetChi2VeeVertex() const {
  return m_Chi2VeeVertex;
}

double FindKS::GetKSMassVeeVertex() const {
  return m_KSMassVeeVertex;
}

double FindKS::GetDecayLengthFit() const {
  return m_DecayLengthFit;
}

double FindKS::GetDecayLengthErrorFit() const {
  return m_DecayLengthErrorFit;
}

double FindKS::GetChi2Fit() const {
  return m_Chi2Fit;
}

double FindKS::GetKSMassFit() const {
  return m_KSMassFit;
}

double FindKS::GetKSPiPlusP(int i) const {
  return m_KSPiPlusP[i];
}

double FindKS::GetKSPiMinusP(int i) const {
  return m_KSPiMinusP[i];
}

double FindKS::GetKSPiPlusPFit(int i) const {
  return m_KSPiPlusPFit[i];
}

double FindKS::GetKSPiMinusPFit(int i) const {
  return m_KSPiMinusPFit[i];
}
