// Martin Duy Tat 28th January 2021, based on code by Yu Zhang

// KKpipi
#include "KKpipi/KKpipiSingleTag.h"
#include "KKpipi/FindKKpipiTagInfo.h"
#include "KKpipi/FindMCInfo.h"
#include "KKpipi/PIDTruth.h"
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

KKpipiSingleTag::KKpipiSingleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKpipiSingleTag::~KKpipiSingleTag() {
}

StatusCode KKpipiSingleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KKpipi Single Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KKpipiSingleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KKpipiSingleTag", CLID_ColumnWiseTuple, "Single tagged D->KKpipi events");
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
      status = m_tuple->addItem("KalmanFitSuccess", m_KalmanFitSuccess);
      status = m_tuple->addItem("KalmanFitChi2", m_KalmanFitChi2);
      status = m_tuple->addItem("PiPluspxKalmanFit", m_PiPluspxKalmanFit);
      status = m_tuple->addItem("PiPluspyKalmanFit", m_PiPluspyKalmanFit);
      status = m_tuple->addItem("PiPluspzKalmanFit", m_PiPluspzKalmanFit);
      status = m_tuple->addItem("PiPlusenergyKalmanFit", m_PiPlusenergyKalmanFit);
      status = m_tuple->addItem("PiMinuspxKalmanFit", m_PiMinuspxKalmanFit);
      status = m_tuple->addItem("PiMinuspyKalmanFit", m_PiMinuspyKalmanFit);
      status = m_tuple->addItem("PiMinuspzKalmanFit", m_PiMinuspzKalmanFit);
      status = m_tuple->addItem("PiMinusenergyKalmanFit", m_PiMinusenergyKalmanFit);
      status = m_tuple->addItem("KPluspxKalmanFit", m_KPluspxKalmanFit);
      status = m_tuple->addItem("KPluspyKalmanFit", m_KPluspyKalmanFit);
      status = m_tuple->addItem("KPluspzKalmanFit", m_KPluspzKalmanFit);
      status = m_tuple->addItem("KPlusenergyKalmanFit", m_KPlusenergyKalmanFit);
      status = m_tuple->addItem("KMinuspxKalmanFit", m_KMinuspxKalmanFit);
      status = m_tuple->addItem("KMinuspyKalmanFit", m_KMinuspyKalmanFit);
      status = m_tuple->addItem("KMinuspzKalmanFit", m_KMinuspzKalmanFit);
      status = m_tuple->addItem("KMinusenergyKalmanFit", m_KMinusenergyKalmanFit);
      status = m_tuple->addItem("Mpipi", m_Mpipi);
      status = m_tuple->addItem("KSFitSuccess", m_KSFitSuccess);
      status = m_tuple->addItem("KSDecayLengthVeeVertex", m_DecayLengthVeeVertex);
      status = m_tuple->addItem("KSChi2VeeVertex", m_Chi2VeeVertex);
      status = m_tuple->addItem("KSMassVeeVertex", m_KSMassVeeVertex);
      status = m_tuple->addItem("KSDecayLengthFit", m_DecayLengthFit);
      status = m_tuple->addItem("KSDecayLengthErrorFit", m_DecayLengthErrorFit);
      status = m_tuple->addItem("KSChi2Fit", m_Chi2Fit);
      status = m_tuple->addItem("IsSameDMother", m_IsSameDMother);
      status = m_tuple->addItem("PIDTrue", m_PIDTrue);
      status = m_tuple->addItem("KPlusTrueID", m_KPlusTrueID);
      status = m_tuple->addItem("KMinusTrueID", m_KMinusTrueID);
      status = m_tuple->addItem("PiPlusTrueID", m_PiPlusTrueID);
      status = m_tuple->addItem("PiMinusTrueID", m_PiMinusTrueID);
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
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
    m_RunNumber = eventHeader->runNumber();
    m_EventNumber = eventHeader->eventNumber();
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

StatusCode KKpipiSingleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KKpipi Single Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KKpipiSingleTag::FillTuple(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
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
  FindKKpipiTagInfo findKKpipiTagInfo;
  StatusCode status = findKKpipiTagInfo.CalculateTagInfo(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_KPluspx = findKKpipiTagInfo.GetKPlusP(0);
  m_KPluspy = findKKpipiTagInfo.GetKPlusP(1);
  m_KPluspz = findKKpipiTagInfo.GetKPlusP(2);
  m_KPlusenergy = findKKpipiTagInfo.GetKPlusP(3);
  m_KMinuspx = findKKpipiTagInfo.GetKMinusP(0);
  m_KMinuspy = findKKpipiTagInfo.GetKMinusP(1);
  m_KMinuspz = findKKpipiTagInfo.GetKMinusP(2);
  m_KMinusenergy = findKKpipiTagInfo.GetKMinusP(3);
  m_PiPluspx = findKKpipiTagInfo.GetPiPlusP(0);
  m_PiPluspy = findKKpipiTagInfo.GetPiPlusP(1);
  m_PiPluspz = findKKpipiTagInfo.GetPiPlusP(2);
  m_PiPlusenergy = findKKpipiTagInfo.GetPiPlusP(3);
  m_PiMinuspx = findKKpipiTagInfo.GetPiMinusP(0);
  m_PiMinuspy = findKKpipiTagInfo.GetPiMinusP(1);
  m_PiMinuspz = findKKpipiTagInfo.GetPiMinusP(2);
  m_PiMinusenergy = findKKpipiTagInfo.GetPiMinusP(3);
  m_KalmanFitSuccess = findKKpipiTagInfo.GetKalmanFitSuccess();
  m_KalmanFitChi2 = findKKpipiTagInfo.GetKalmanFitChi2();
  m_KPluspxKalmanFit = findKKpipiTagInfo.GetKPlusPKalmanFit(0);
  m_KPluspyKalmanFit = findKKpipiTagInfo.GetKPlusPKalmanFit(1);
  m_KPluspzKalmanFit = findKKpipiTagInfo.GetKPlusPKalmanFit(2);
  m_KPlusenergyKalmanFit = findKKpipiTagInfo.GetKPlusPKalmanFit(3);
  m_KMinuspxKalmanFit = findKKpipiTagInfo.GetKMinusPKalmanFit(0);
  m_KMinuspyKalmanFit = findKKpipiTagInfo.GetKMinusPKalmanFit(1);
  m_KMinuspzKalmanFit = findKKpipiTagInfo.GetKMinusPKalmanFit(2);
  m_KMinusenergyKalmanFit = findKKpipiTagInfo.GetKMinusPKalmanFit(3);
  m_PiPluspxKalmanFit = findKKpipiTagInfo.GetPiPlusPKalmanFit(0);
  m_PiPluspyKalmanFit = findKKpipiTagInfo.GetPiPlusPKalmanFit(1);
  m_PiPluspzKalmanFit = findKKpipiTagInfo.GetPiPlusPKalmanFit(2);
  m_PiPlusenergyKalmanFit = findKKpipiTagInfo.GetPiPlusPKalmanFit(3);
  m_PiMinuspxKalmanFit = findKKpipiTagInfo.GetPiMinusPKalmanFit(0);
  m_PiMinuspyKalmanFit = findKKpipiTagInfo.GetPiMinusPKalmanFit(1);
  m_PiMinuspzKalmanFit = findKKpipiTagInfo.GetPiMinusPKalmanFit(2);
  m_PiMinusenergyKalmanFit = findKKpipiTagInfo.GetPiMinusPKalmanFit(3);
  m_Mpipi = findKKpipiTagInfo.GetMpipi();
  m_KSFitSuccess = findKKpipiTagInfo.GetKSFitSuccess();
  m_DecayLengthVeeVertex = findKKpipiTagInfo.GetDecayLengthVeeVertex();
  m_Chi2VeeVertex = findKKpipiTagInfo.GetChi2VeeVertex();
  m_KSMassVeeVertex = findKKpipiTagInfo.GetKSMassVeeVertex();
  m_DecayLengthFit = findKKpipiTagInfo.GetDecayLengthFit();
  m_DecayLengthErrorFit = findKKpipiTagInfo.GetDecayLengthErrorFit();
  m_Chi2Fit = findKKpipiTagInfo.GetChi2Fit();
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKKpipiTagInfo.GetDaughterTrackID(), 4, this);
    m_IsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {321, -321, 211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_PIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_KPlusTrueID = ReconstructedPID[0];
    m_KMinusTrueID = ReconstructedPID[1];
    m_PiPlusTrueID = ReconstructedPID[2];
    m_PiMinusTrueID = ReconstructedPID[3];
  }
  return StatusCode::SUCCESS;
}
