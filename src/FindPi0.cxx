// Martin Duy Tat 12th February 2021, based on code by Yu Zhang

// Header file
#include "KKpipi/FindPi0.h"
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
#include "EmcRecEventModel/RecEmcShower.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"
#include "EvtRecEvent/EvtRecDTag.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecPi0.h"
#include "EvtRecEvent/EvtRecTrack.h"
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
// ROOT
#include "TMath.h"
// Particle masses
#include "KKpipi/ParticleMasses.h"

FindPi0::FindPi0(int npi0): m_Chi2Fit(0.0), m_npi0(npi0) {
  m_HighEPhotonP.resize(m_npi0);
  m_LowEPhoton.resize(m_npi0);
  m_HighEPhotonPConstrained.resize(m_npi0);
  m_LowEPhotonConstrained.resize(m_npi0);
  m_Chi2Fit.resize(m_npi0);
}

FindPi0::~FindPi0() {
}

CLHEP::HepLorentzVector FindPi0::GetPhoton4Vector(double Energy, double Theta, double Phi) const {
  double Px = Energy*TMath::Sin(Theta)*TMath::Cos(Phi);
  double Py = Energy*TMath::Sin(Theta)*TMath::Sin(Phi);
  double Pz = Energy*TMath::Cos(Theta);
  return CLHEP::HepLorentzVector(Px, Py, Pz, Energy);
}

StatusCode FindPi0::findPi0(DTagToolIterator &DTTool_iter, DTagTool &DTTool) {
  IMessageSvc *msgSvc;
  Gaudi::svcLocator()->service("MessageSvc", msgSvc);
  MsgStream log(msgSvc, "FindPi0");
  IDataProviderSvc *EventService = nullptr;
  Gaudi::svcLocator()->service("EventDataSvc", EventService);
  SmartDataPtr<EvtRecPi0Col> RecPi0Col(EventService, "/Event/EvtRec/EvtRecPi0Col");
  if(!RecPi0Col) {
    log << "Could not find EvtRecPi0Col" << endreq;
    return StatusCode::FAILURE;
  }
  // Get the track ID of the tagged pi0 candidates
  std::vector<int> Pi0TrackIDs = DTTool.pi0Id(DTTool_iter, 1);
  // Loop over all pi0 candidates found by EvtRecPi0
  int i = 0;
  for(EvtRecPi0Col::iterator Pi0_iter = RecPi0Col->begin(); Pi0_iter != RecPi0Col->end(); Pi0_iter++) {
    // Check if pi0 candidate from the DTagTool is the same as the candidate found by EvtRecPi0
    if(Pi0_iter - RecPi0Col->begin() != Pi0TrackIDs[0]) {
      continue;
    }
    // Get EM shower four-momenta of photons
    RecEmcShower *HighEPhotonShower = const_cast<EvtRecTrack*>((*Pi0_iter)->hiEnGamma())->emcShower();
    RecEmcShower *LowEPhotonShower = const_cast<EvtRecTrack*>((*Pi0_iter)->loEnGamma())->emcShower();
    m_HighEPhotonP[i] = GetPhoton4Vector(HighEPhotonShower->energy(), HighEPhotonShower->theta(), HighEPhotonShower->phi());
    m_LowEPhotonP[i] = GetPhoton4Vector(LowEPhotonShower->energy(), LowEPhotonShower->theta(), LowEPhotonShower->phi());
    // Get kinematically constrained four-momenta of photons
    m_HighEPhotonPConstrained[i] = (*Pi0_iter)->hiPfit();
    m_LowEPhotonPConstrained[i] = (*Pi0_iter)->loPfit();
    m_Chi2Fit[i] = (*Pi0_iter)->chisq();
    i++;
    if(i == m_npi0) {
      return StatusCode::SUCCESS;
    }
  }
  log << "Could not find the requested number of pi0 candidates" << endreq;
  return StatusCode::FAILURE;
}

double FindPi0::GetHighEPhotonP(int i, int pi0_index) const {
  return m_HighEPhotonP[pi0_index][i];
}

double FindPi0::GetLowEPhotonP(int i, int pi0_index) const {
  return m_LowEPhotonP[pi0_index][i];
}

double FindPi0::GetHighEPhotonPConstrained(int i, int pi0_index) const {
  return m_HighEPhotonPConstrained[pi0_index][i];
}

double FindPi0::GetLowEPhotonPConstrained(int i, int pi0_index) const {
  return m_LowEPhotonPConstrained[pi0_index][i];
}

double FindPi0::GetChi2Fit(int pi0_index) const {
  return m_Chi2Fit[pi0_index];
}
