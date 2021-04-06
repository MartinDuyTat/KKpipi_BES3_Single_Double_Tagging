// Martin Duy Tat 12th March 2021

// KKpipi
#include "KKpipi/KSpipiSingleTag.h"
#include "KKpipi/FindKSpipiTagInfo.h"
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

KSpipiSingleTag::KSpipiSingleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KSpipiSingleTag::~KSpipiSingleTag() {
}

StatusCode KSpipiSingleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KSpipi Single Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KSpipiSingleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KSpipiSingleTag", CLID_ColumnWiseTuple, "Single tagged D->KSpipi events");
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
      status = m_tuple->addItem("KSDecayLengthVeeVertex", m_DecayLengthVeeVertex);
      status = m_tuple->addItem("KSChi2VeeVertex", m_Chi2VeeVertex);
      status = m_tuple->addItem("KSMassVeeVertex", m_KSMassVeeVertex);
      status = m_tuple->addItem("KSDecayLengthFit", m_DecayLengthFit);
      status = m_tuple->addItem("KSDecayLengthErrorFit", m_DecayLengthErrorFit);
      status = m_tuple->addItem("KSChi2Fit", m_Chi2Fit);
      status = m_tuple->addItem("KSMassFit", m_KSMassFit);
      status = m_tuple->addItem("KSpx", m_KSpx);
      status = m_tuple->addItem("KSpy", m_KSpy);
      status = m_tuple->addItem("KSz", m_KSpz);
      status = m_tuple->addItem("KSenergy", m_KSenergy);
      status = m_tuple->addItem("KSPiPluspx", m_KSPiPluspx);
      status = m_tuple->addItem("KSPiPluspy", m_KSPiPluspy);
      status = m_tuple->addItem("KSPiPluspz", m_KSPiPluspz);
      status = m_tuple->addItem("KSPiPlusenergy", m_KSPiPlusenergy);
      status = m_tuple->addItem("KSPiMinuspx", m_KSPiMinuspx);
      status = m_tuple->addItem("KSPiMinuspy", m_KSPiMinuspy);
      status = m_tuple->addItem("KSPiMinuspz", m_KSPiMinuspz);
      status = m_tuple->addItem("KSPiMinusenergy", m_KSPiMinusenergy);
      status = m_tuple->addItem("PiPluspx", m_PiPluspx);
      status = m_tuple->addItem("PiPluspy", m_PiPluspy);
      status = m_tuple->addItem("PiPluspz", m_PiPluspz);
      status = m_tuple->addItem("PiPlusenergy", m_PiPlusenergy);
      status = m_tuple->addItem("PiMinuspx", m_PiMinuspx);
      status = m_tuple->addItem("PiMinuspy", m_PiMinuspy);
      status = m_tuple->addItem("PiMinuspz", m_PiMinuspz);
      status = m_tuple->addItem("PiMinusenergy", m_PiMinusenergy);
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
      status = m_tuple->addItem("KSpxKalmanFit", m_KSpxKalmanFit);
      status = m_tuple->addItem("KSpyKalmanFit", m_KSpyKalmanFit);
      status = m_tuple->addItem("KSpzKalmanFit", m_KSpzKalmanFit);
      status = m_tuple->addItem("KSenergyKalmanFit", m_KSenergyKalmanFit);
      status = m_tuple->addItem("pipiKSFitSuccess", m_pipiKSFitSuccess);
      status = m_tuple->addItem("pipiKSDecayLengthVeeVertex", m_pipiDecayLengthVeeVertex);
      status = m_tuple->addItem("pipiKSChi2VeeVertex", m_pipiChi2VeeVertex);
      status = m_tuple->addItem("pipiKSMassVeeVertex", m_pipiKSMassVeeVertex);
      status = m_tuple->addItem("pipiKSDecayLengthFit", m_pipiDecayLengthFit);
      status = m_tuple->addItem("pipiKSDecayLengthErrorFit", m_pipiDecayLengthErrorFit);
      status = m_tuple->addItem("pipiKSChi2Fit", m_pipiChi2Fit);
      status = m_tuple->addItem("pipiKSMassFit", m_pipiKSMassFit);
      status = m_tuple->addItem("IsSameDMother", m_IsSameDMother);
      status = m_tuple->addItem("PIDTrue", m_PIDTrue);
      status = m_tuple->addItem("KSPiPlusTrueID", m_KSPiPlusTrueID);
      status = m_tuple->addItem("KSPiMinusTrueID", m_KSPiMinusTrueID);
      status = m_tuple->addItem("PiPlusTrueID", m_PiPlusTrueID);
      status = m_tuple->addItem("PiMinusTrueID", m_PiMinusTrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KSpipi Single Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KSpipiSingleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KSpipi Single Tag Algorithm" << endreq;
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
  if(DTTool.findSTag(EvtRecDTag::kD0toKsPiPi)) {
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
    m_RunNumber = eventHeader->runNumber();
    m_EventNumber = eventHeader->eventNumber();
    DTagToolIterator DTTool_iter = DTTool.stag();
    StatusCode FillTupleStatus = FillTuple(DTTool_iter, DTTool);
    if(FillTupleStatus != StatusCode::SUCCESS) {
      if(FillTupleStatus == StatusCode::RECOVERABLE) {
	log << MSG::WARNING << "Vertex fit of KS failed, skipping event" << endreq;
	return StatusCode::SUCCESS;
      }
      log << MSG::FATAL << "Assigning tuple info failed" << endreq;
      return StatusCode::FAILURE;
    }
    m_tuple->write();
  }
  return StatusCode::SUCCESS;
}

StatusCode KSpipiSingleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KSpipi Single Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KSpipiSingleTag::FillTuple(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
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
  FindKSpipiTagInfo findKSpipiTagInfo;
  StatusCode status = findKSpipiTagInfo.CalculateTagInfo(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_DecayLengthVeeVertex = findKSpipiTagInfo.GetDecayLengthVeeVertex();
  m_Chi2VeeVertex = findKSpipiTagInfo.GetChi2VeeVertex();
  m_KSMassVeeVertex = findKSpipiTagInfo.GetKSMassVeeVertex();
  m_DecayLengthFit = findKSpipiTagInfo.GetDecayLengthFit();
  m_DecayLengthErrorFit = findKSpipiTagInfo.GetDecayLengthErrorFit();
  m_Chi2Fit = findKSpipiTagInfo.GetChi2Fit();
  m_KSMassFit = findKSpipiTagInfo.GetKSMassFit();
  m_KSPiPluspx = findKSpipiTagInfo.GetKSPiPlusP(0);
  m_KSPiPluspy = findKSpipiTagInfo.GetKSPiPlusP(1);
  m_KSPiPluspz = findKSpipiTagInfo.GetKSPiPlusP(2);
  m_KSPiPlusenergy = findKSpipiTagInfo.GetKSPiPlusP(3);
  m_KSPiMinuspx = findKSpipiTagInfo.GetKSPiMinusP(0);
  m_KSPiMinuspy = findKSpipiTagInfo.GetKSPiMinusP(1);
  m_KSPiMinuspz = findKSpipiTagInfo.GetKSPiMinusP(2);
  m_KSPiMinusenergy = findKSpipiTagInfo.GetKSPiMinusP(3);
  m_KSpx = findKSpipiTagInfo.GetKShortP(0);
  m_KSpy = findKSpipiTagInfo.GetKShortP(1);
  m_KSpz = findKSpipiTagInfo.GetKShortP(2);
  m_KSenergy = findKSpipiTagInfo.GetKShortP(3);
  m_PiPluspx = findKSpipiTagInfo.GetPiPlusP(0);
  m_PiPluspy = findKSpipiTagInfo.GetPiPlusP(1);
  m_PiPluspz = findKSpipiTagInfo.GetPiPlusP(2);
  m_PiPlusenergy = findKSpipiTagInfo.GetPiPlusP(3);
  m_PiMinuspx = findKSpipiTagInfo.GetPiMinusP(0);
  m_PiMinuspy = findKSpipiTagInfo.GetPiMinusP(1);
  m_PiMinuspz = findKSpipiTagInfo.GetPiMinusP(2);
  m_PiMinusenergy = findKSpipiTagInfo.GetPiMinusP(3);
  m_KalmanFitSuccess = findKSpipiTagInfo.GetKalmanFitSuccess();
  m_KalmanFitChi2 = findKSpipiTagInfo.GetKalmanFitChi2();
  m_KSpxKalmanFit = findKSpipiTagInfo.GetKShortPKalmanFit(0);
  m_KSpyKalmanFit = findKSpipiTagInfo.GetKShortPKalmanFit(1);
  m_KSpzKalmanFit = findKSpipiTagInfo.GetKShortPKalmanFit(2);
  m_KSenergyKalmanFit = findKSpipiTagInfo.GetKShortPKalmanFit(3);
  m_PiPluspxKalmanFit = findKSpipiTagInfo.GetPiPlusPKalmanFit(0);
  m_PiPluspyKalmanFit = findKSpipiTagInfo.GetPiPlusPKalmanFit(1);
  m_PiPluspzKalmanFit = findKSpipiTagInfo.GetPiPlusPKalmanFit(2);
  m_PiPlusenergyKalmanFit = findKSpipiTagInfo.GetPiPlusPKalmanFit(3);
  m_PiMinuspxKalmanFit = findKSpipiTagInfo.GetPiMinusPKalmanFit(0);
  m_PiMinuspyKalmanFit = findKSpipiTagInfo.GetPiMinusPKalmanFit(1);
  m_PiMinuspzKalmanFit = findKSpipiTagInfo.GetPiMinusPKalmanFit(2);
  m_PiMinusenergyKalmanFit = findKSpipiTagInfo.GetPiMinusPKalmanFit(3);
  m_pipiKSFitSuccess = findKSpipiTagInfo.GetpipiKSFitSuccess();
  m_pipiDecayLengthVeeVertex = findKSpipiTagInfo.GetpipiDecayLengthVeeVertex();
  m_pipiChi2VeeVertex = findKSpipiTagInfo.GetpipiChi2VeeVertex();
  m_pipiKSMassVeeVertex = findKSpipiTagInfo.GetpipiKSMassVeeVertex();
  m_pipiDecayLengthFit = findKSpipiTagInfo.GetpipiDecayLengthFit();
  m_pipiDecayLengthErrorFit = findKSpipiTagInfo.GetpipiDecayLengthErrorFit();
  m_pipiChi2Fit = findKSpipiTagInfo.GetpipiChi2Fit();
  m_pipiKSMassFit = findKSpipiTagInfo.GetpipiKSMassFit();
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKSpipiTagInfo.GetDaughterTrackID(), this);
    m_IsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {211, -211, 211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_PIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_KSPiPlusTrueID = ReconstructedPID[0];
    m_KSPiMinusTrueID = ReconstructedPID[1];
    m_PiPlusTrueID = ReconstructedPID[2];
    m_PiMinusTrueID = ReconstructedPID[3];
  }
  return StatusCode::SUCCESS;
}
