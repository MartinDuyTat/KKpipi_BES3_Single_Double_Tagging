// Martin Duy Tat 28th January 2021

// Header file
#include "KKpipi/KKpipiSingleTag.h"
// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartDataPtr.h"
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
#include "SimplePIDSvc/ISimplePIDSvc.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McTruth/McParticle.h"
// STL
#include<vector>
#include<string>

KKpipiSingleTag::KKpipiSingleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKpipiSingleTag::~KKpipiSingleTag() {
}

StatusCode KKpipiSingleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KKpipi Single Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/SingleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/SingleTag", CLID_ColumnWiseTuple, "Single tagged D->KKpipi events");
    if(m_tuple) {
      status = m_tuple->addItem("Run", m_RunNumber);
      status = m_tuple->addItem("Event", m_EventNumber);
      status = m_tuple->addItem("NumberOfParticles", m_NumberParticles, 0, 100);
      status = m_tuple->addIndexedItem("ParticleIDs", m_NumberParticles, m_pdgID);
      status = m_tuple->addIndexedItem("MotherIndex", m_NumberParticles, m_MotherIndex);
      status = m_tuple->addItem("MCmode", m_MCmode);
      status = m_tuple->addItem("GeneratorNumberOfParticles", m_GeneratorNumberParticles, 0, 100);
      status = m_tuple->addIndexedItem("GeneratorParticleIDs", m_GeneratorNumberParticles, m_GeneratorPDGID);
      status = m_tuple->addIndexedItem("GeneratorMotherID", m_GeneratorNumberParticles, m_MotherID);
      status = m_tuple->addIndexedItem("True_P", m_GeneratorNumberParticles, m_TrueMomentum);
      status = m_tuple->addIndexedItem("True_PT", m_GeneratorNumberParticles, m_TruePT);
      status = m_tuple->addIndexedItem("True_phi", m_GeneratorNumberParticles, m_TruePhi);
      status = m_tuple->addIndexedItem("True_theta", m_GeneratorNumberParticles, m_TrueTheta);
      status = m_tuple->addItem("Charm", m_Charm);
      status = m_tuple->addItem("DMass", m_DMass);
      status = m_tuple->addItem("MBC", m_MBC);
      status = m_tuple->addItem("DeltaE", m_DeltaE);
      status = m_tuple->addItem("BeamE", m_BeamE);
      status = m_tuple->addItem("Dpx", m_Dpx);
      status = m_tuple->addItem("Dpy", m_Dpy);
      status = m_tuple->addItem("Dpz", m_Dpz);
      status = m_tuple->addItem("Denergy", m_Denergy);
      status = m_tuple->addItem("PiPluspx", m_PiPluspx);
      status = m_tuple->addItem("PiPluspy", m_PiPluspy);
      status = m_tuple->addItem("PiPluspz", m_PiPluspz);
      status = m_tuple->addItem("PiPlusenergy", m_PiPlusenergy);
      status = m_tuple->addItem("PiMinuspx", m_PiMinuspx);
      status = m_tuple->addItem("PiMinuspy", m_PiMinuspy);
      status = m_tuple->addItem("PiMinuspz", m_PiMinuspz);
      status = m_tuple->addItem("PiMinusenergy", m_PiMinusenergy);
      status = m_tuple->addItem("KPluspx", m_KPluspx);
      status = m_tuple->addItem("KPluspy", m_KPluspy);
      status = m_tuple->addItem("KPluspz", m_KPluspz);
      status = m_tuple->addItem("KPlusenergy", m_KPlusenergy);
      status = m_tuple->addItem("KMinuspx", m_KMinuspx);
      status = m_tuple->addItem("KMinuspy", m_KMinuspy);
      status = m_tuple->addItem("KMinuspz", m_KMinuspz);
      status = m_tuple->addItem("KMinusenergy", m_KMinusenergy);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KKpipi Single Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KKpipiSingleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KKpipi Single Tag Algorithm" << endreq;
  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
  m_RunNumber = eventHeader->runNumber();
  m_EventNumber = eventHeader->eventNumber();
  if(m_RunNumber < 0) {
    SmartDataPtr<Event::McParticleCol> MCParticleCol(eventSvc(), "/Event/MC/McParticleCol");
    if(!MCParticleCol) {
      log << MSG::FATAL << "Could not load McParticleCol" << endreq;
      return StatusCode::FAILURE;
    }
    std::vector<int> pdgID, MotherIndex;
    int ParticleNumber = 0;
    for(Event::McParticleCol::iterator MCParticleCol_iter = MCParticleCol->begin(); MCParticleCol_iter != MCParticleCol->end(); MCParticleCol_iter++) {
      if((*MCParticleCol_iter)->primaryParticle()) {
	continue;
      } else if((*MCParticleCol_iter)->decayFromGenerator()) {
	continue;
      } else if((*MCParticleCol_iter)->particleProperty() == 30443) {
	CLHEP::HepLorentzVector initialP = (*MCParticleCol_iter)->initialFourMomentum();
	m_TrueMomentum[ParticleNumber] = initialP.mag();
	m_TruePT[ParticleNumber] = initialP.perp();
	m_TruePhi[ParticleNumber] = initialP.phi();
	m_TrueTheta[ParticleNumber] = initialP.cosTheta();
	m_GeneratorPDGID[ParticleNumber] = (*MCParticleCol_iter)->particleProperty();
	m_MotherID[ParticleNumber] = (*MCParticleCol_iter)->mother().particleProperty();
	++ParticleNumber;
	IMcDecayModeSvc *IMcDecayModeService;
	StatusCode McDecayModeSVC_Status = service("McDecayModeSvc", IMcDecayModeService);
	if(McDecayModeSVC_Status.isFailure()) {
	  log << MSG::FATAL << "Could not load McDecayModeSvc" << endreq;
	  return McDecayModeSVC_Status;
	}
	McDecayModeSvc *McDecayModeService = dynamic_cast<McDecayModeSvc*>(IMcDecayModeService);
	m_MCmode = McDecayModeService->extract(*MCParticleCol_iter, pdgID, MotherIndex);
      }
    }
    m_GeneratorNumberParticles = ParticleNumber;
    m_NumberParticles = pdgID.size();
    for(int i = 0; i < m_NumberParticles; i++) {
      m_pdgID[i] = pdgID[i];
      m_MotherIndex[i] = MotherIndex[i];
    }
  }
  DTagTool DTTool;
  DTTool.setPID(true);
  if(DTTool.isDTagListEmpty()) {
    log << MSG::DEBUG << "No D candidates found" << endreq;
    return StatusCode::SUCCESS;
  }
  if(!DTTool.cosmicandleptonVeto()) {
    log << MSG::DEBUG << "Cosmic and lepton veto" << endreq;
    return StatusCode::SUCCESS;
  }
  if(DTTool.findSTag(EvtRecDTag::kD0toKKPiPi)) {
    DTagToolIterator DTTool_iter = DTTool.stag();
    AssignTagInfo(DTTool_iter);
    m_tuple->write();
  }
  return StatusCode::SUCCESS;
}

StatusCode KKpipiSingleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KKpipi Single Tagging" << endreq;
  return StatusCode::SUCCESS;
}

void KKpipiSingleTag::AssignTagInfo(DTagToolIterator DTTool_iter) {
  m_Charm = (*DTTool_iter)->charm();
  m_DMass = (*DTTool_iter)->mass();
  m_MBC = (*DTTool_iter)->mBC();
  m_DeltaE = (*DTTool_iter)->deltaE();
  m_BeamE = (*DTTool_iter)->beamE();
  m_Dpx = (*DTTool_iter)->p4().x();
  m_Dpy = (*DTTool_iter)->p4().y();
  m_Dpz = (*DTTool_iter)->p4().z();
  m_Denergy = (*DTTool_iter)->p4().t();
  m_Charm = (*DTTool_iter)->charm();
  m_Charm = (*DTTool_iter)->charm();
  m_Charm = (*DTTool_iter)->charm();
}
