// Martin Duy Tat 28th January 2021, based on code by Yu Zhang

// KKpipi
#include "KKpipi/FindMCInfo.h"
// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/SmartRefVector.h"
// Event information
#include "EventModel/Event.h"
#include "EventModel/EventModel.h"
#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecDTag.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McTruth/McParticle.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
// STL
#include<vector>
// Particle masses
#include "KKpipi/ParticleMasses.h"

FindMCInfo::FindMCInfo(): m_NumberParticles(0), m_MCmode(0) {
}

FindMCInfo::~FindMCInfo() {
}

StatusCode CalculateMCInfo() {
  SmartDataPtr<Event::McParticleCol> MCParticleCol(eventSvc(), "/Event/MC/McParticleCol");
  if(!MCParticleCol) {
    log << MSG::FATAL << "Could not load McParticleCol" << endreq;
    return StatusCode::FAILURE;
  }
  std::vector<int> pdgID, MotherIndex;
  std::vector<double> TruePx, TruePy, TruePz, TrueEnergy;
  for(Event::McParticleCol::iterator MCParticleCol_iter = MCParticleCol->begin(); MCParticleCol_iter != MCParticleCol->end(); MCParticleCol_iter++) {
    if((*MCParticleCol_iter)->primaryParticle() || !(*MCParticleCol_iter)->decayFromGenerator()) {
      continue;
    }
    if((*MCParticleCol_iter)->particleProperty() == 30443) {
      
      IMcDecayModeSvc *IMcDecayModeService;
      StatusCode McDecayModeSVC_Status = service("McDecayModeSvc", IMcDecayModeService);
      if(McDecayModeSVC_Status.isFailure()) {
	log << MSG::FATAL << "Could not load McDecayModeSvc" << endreq;
	return McDecayModeSVC_Status;
      }
      McDecayModeSvc *McDecayModeService = dynamic_cast<McDecayModeSvc*>(IMcDecayModeService);
      m_MCmode = McDecayModeService->extract(*MCParticleCol_iter, pdgID, MotherIndex, TruePx, TruePy, TruePz, TrueEnergy);
    }
  }
  m_NumberParticles = pdgID.size();
  for(int i = 0; i < m_NumberParticles; i++) {
    m_pdgID[i] = pdgID[i];
    m_MotherIndex[i] = MotherIndex[i];
    m_TruePx[i] = TruePx[i];
    m_TruePy[i] = TruePy[i];
    m_TruePz[i] = TruePz[i];
    m_TrueEnergy[i] = TrueEnergy[i];
  }
  return StatusCode::SUCCESS;
}

int GetNumberParticles() {
  return m_NumberParticles;
}

std::vector<int> GetpdgID() {
  return m_pdgID;
}

std::vector<int> GetMotherIndex() {
  return m_MotherIndex;
}

int GetMCmode() {
  return m_TruePx;
}

std::vector<int> GetTruePx() {
  return m_TruePx;
}

std::vector<int> GetTruePy() {
  return m_True_Py;
}

std::vector<int> GetTruePz() {
  return m_TruePz;
}

std::vector<int> GetTrueEnergy() {
  return m_TrueEnergy;
}
