// Martin Duy Tat 12th March 2021

// KKpipi
#include "KKpipi/KKSingleTag.h"
#include "KKpipi/FindhhTagInfo.h"
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
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McTruth/McParticle.h"
// STL
#include<vector>
#include<string>

KKSingleTag::KKSingleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKSingleTag::~KKSingleTag() {
}

StatusCode KKSingleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KK Single Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KKSingleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KKSingleTag", CLID_ColumnWiseTuple, "Single tagged D->KK events");
    if(m_tuple) {
      status = m_tuple->addItem("Run", m_RunNumber);
      status = m_tuple->addItem("Event", m_EventNumber);
      status = m_tuple->addItem("NumberOfParticles", m_NumberParticles, 0, 100);
      status = m_tuple->addIndexedItem("ParticleIDs", m_NumberParticles, m_pdgID);
      status = m_tuple->addIndexedItem("MotherIndex", m_NumberParticles, m_MotherIndex);
      status = m_tuple->addItem("NumberOfParticlesStripped", m_NumberParticlesStripped, 0, 100);
      status = m_tuple->addIndexedItem("ParticleIDsStripped", m_NumberParticlesStripped, m_pdgIDStripped);
      status = m_tuple->addIndexedItem("MotherIndexStripped", m_NumberParticlesStripped, m_MotherIndexStripped);
      status = m_tuple->addItem("MCmode", m_MCmode);
      status = m_tuple->addIndexedItem("True_Px", m_NumberParticles, m_TruePx);
      status = m_tuple->addIndexedItem("True_Py", m_NumberParticles, m_TruePy);
      status = m_tuple->addIndexedItem("True_Pz", m_NumberParticles, m_TruePz);
      status = m_tuple->addIndexedItem("True_Energy", m_NumberParticles, m_TrueEnergy);
      status = m_tuple->addItem("DMass", m_DMass);
      status = m_tuple->addItem("MBC", m_MBC);
      status = m_tuple->addItem("DeltaE", m_DeltaE);
      status = m_tuple->addItem("BeamE", m_BeamE);
      status = m_tuple->addItem("Dpx", m_Dpx);
      status = m_tuple->addItem("Dpy", m_Dpy);
      status = m_tuple->addItem("Dpz", m_Dpz);
      status = m_tuple->addItem("Denergy", m_Denergy);
      status = m_tuple->addItem("KPluspx", m_KPluspx);
      status = m_tuple->addItem("KPluspy", m_KPluspy);
      status = m_tuple->addItem("KPluspz", m_KPluspz);
      status = m_tuple->addItem("KPlusenergy", m_KPlusenergy);
      status = m_tuple->addItem("KMinuspx", m_KMinuspx);
      status = m_tuple->addItem("KMinuspy", m_KMinuspy);
      status = m_tuple->addItem("KMinuspz", m_KMinuspz);
      status = m_tuple->addItem("KMinusenergy", m_KMinusenergy);
      status = m_tuple->addItem("IsSameDMother", m_IsSameDMother);
      status = m_tuple->addItem("PIDTrue", m_PIDTrue);
      status = m_tuple->addItem("KPlusTrueID", m_KPlusTrueID);
      status = m_tuple->addItem("KMinusTrueID", m_KMinusTrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KK Single Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KKSingleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KK Single Tag Algorithm" << endreq;
  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
  m_RunNumber = eventHeader->runNumber();
  m_EventNumber = eventHeader->eventNumber();
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
  if(DTTool.findSTag(EvtRecDTag::kD0toKK)) {
    DTagToolIterator DTTool_iter = DTTool.stag();
    StatusCode FillTupleStatus = FillTuple(DTTool_iter, DTTool);
    if(FillTupleStatus != StatusCode::SUCCESS) {
      log << MSG::FATAL << "Assigning tuple info failed" << endreq;
      return StatusCode::FAILURE;
    }
    m_tuple->write();
  }
  return StatusCode::SUCCESS;
}

StatusCode KKSingleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KK Single Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KKSingleTag::FillTuple(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
  if(m_RunNumber < 0) {
    SmartDataPtr<Event::McParticleCol> MCParticleCol(eventSvc(), "/Event/MC/McParticleCol");
    if(!MCParticleCol) {
      return StatusCode::FAILURE;
    }
    IMcDecayModeSvc *IMcDecayModeService;
    StatusCode McDecayModeSVC_Status = service("McDecayModeSvc", IMcDecayModeService);
    if(McDecayModeSVC_Status.isFailure()) {
      return StatusCode::FAILURE;
    }
    FindMCInfo findMCInfo;
    StatusCode MCStatus = findMCInfo.CalculateMCInfo(MCParticleCol, IMcDecayModeService);
    if(MCStatus != StatusCode::SUCCESS) {
      return MCStatus;
    }
    m_NumberParticles = findMCInfo.GetNumberParticles();
    m_MCmode = findMCInfo.GetMCmode();
    for(int i = 0; i < m_NumberParticles; i++) {
      m_pdgID[i] = findMCInfo.GetpdgID(i);
      m_MotherIndex[i] = findMCInfo.GetMotherIndex(i);
      m_TruePx[i] = findMCInfo.GetTruePx(i);
      m_TruePy[i] = findMCInfo.GetTruePy(i);
      m_TruePz[i] = findMCInfo.GetTruePz(i);
      m_TrueEnergy[i] = findMCInfo.GetTrueEnergy(i);
    }
    m_NumberParticlesStripped = findMCInfo.GetNumberParticlesStripped();
    for(int i = 0; i < m_NumberParticlesStripped; i++) {
      m_pdgIDStripped[i] = findMCInfo.GetpdgIDStripped(i);
      m_MotherIndexStripped[i] = findMCInfo.GetMotherIndexStripped(i);
    }
  }
  m_DMass = (*DTTool_iter)->mass();
  m_MBC = (*DTTool_iter)->mBC();
  m_DeltaE = (*DTTool_iter)->deltaE();
  m_BeamE = (*DTTool_iter)->beamE();
  m_Dpx = (*DTTool_iter)->p4().x();
  m_Dpy = (*DTTool_iter)->p4().y();
  m_Dpz = (*DTTool_iter)->p4().z();
  m_Denergy = (*DTTool_iter)->p4().t();
  FindhhTagInfo findKKTagInfo("KK");
  StatusCode status = findKKTagInfo.CalculateTagInfo(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_KPluspx = findKKTagInfo.GethPlusP(0);
  m_KPluspy = findKKTagInfo.GethPlusP(1);
  m_KPluspz = findKKTagInfo.GethPlusP(2);
  m_KPlusenergy = findKKTagInfo.GethPlusP(3);
  m_KMinuspx = findKKTagInfo.GethMinusP(0);
  m_KMinuspy = findKKTagInfo.GethMinusP(1);
  m_KMinuspz = findKKTagInfo.GethMinusP(2);
  m_KMinusenergy = findKKTagInfo.GethMinusP(3);
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKKTagInfo.GetDaughterTrackID(), this);
    m_IsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[2] = {321, -321};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 2);
    m_PIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_KPlusTrueID = ReconstructedPID[0];
    m_KMinusTrueID = ReconstructedPID[1];
  }
  return StatusCode::SUCCESS;
}
