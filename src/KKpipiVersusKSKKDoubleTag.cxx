// Martin Duy Tat 4th March 2021

// KKpipi
#include "KKpipi/KKpipiVersusKSKKDoubleTag.h"
#include "KKpipi/FindKKpipiTagInfo.h"
#include "KKpipi/FindKS.h"
#include "KKpipi/FindKShhTagInfo.h"
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

KKpipiVersusKSKKDoubleTag::KKpipiVersusKSKKDoubleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKpipiVersusKSKKDoubleTag::~KKpipiVersusKSKKDoubleTag() {
}

StatusCode KKpipiVersusKSKKDoubleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KKpipi vs KSKK Double Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KSKKDoubleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KSKKDoubleTag", CLID_ColumnWiseTuple, "Double tagged D->KKpipi vs D->KSKK events");
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
      status = m_tuple->addItem("SignalDMass", m_SignalDMass);
      status = m_tuple->addItem("SignalMBC", m_SignalMBC);
      status = m_tuple->addItem("SignalDeltaE", m_SignalDeltaE);
      status = m_tuple->addItem("SignalBeamE", m_SignalBeamE);
      status = m_tuple->addItem("TagDMass", m_TagDMass);
      status = m_tuple->addItem("TagMBC", m_TagMBC);
      status = m_tuple->addItem("TagDeltaE", m_TagDeltaE);
      status = m_tuple->addItem("TagBeamE", m_TagBeamE);
      status = m_tuple->addItem("SignalDpx", m_SignalDpx);
      status = m_tuple->addItem("SignalDpy", m_SignalDpy);
      status = m_tuple->addItem("SignalDpz", m_SignalDpz);
      status = m_tuple->addItem("SignalDenergy", m_SignalDenergy);
      status = m_tuple->addItem("TagDpx", m_TagDpx);
      status = m_tuple->addItem("TagDpy", m_TagDpy);
      status = m_tuple->addItem("TagDpz", m_TagDpz);
      status = m_tuple->addItem("TagDenergy", m_TagDenergy);
      status = m_tuple->addItem("SignalPiPluspx", m_SignalPiPluspx);
      status = m_tuple->addItem("SignalPiPluspy", m_SignalPiPluspy);
      status = m_tuple->addItem("SignalPiPluspz", m_SignalPiPluspz);
      status = m_tuple->addItem("SignalPiPlusenergy", m_SignalPiPlusenergy);
      status = m_tuple->addItem("SignalPiMinuspx", m_SignalPiMinuspx);
      status = m_tuple->addItem("SignalPiMinuspy", m_SignalPiMinuspy);
      status = m_tuple->addItem("SignalPiMinuspz", m_SignalPiMinuspz);
      status = m_tuple->addItem("SignalPiMinusenergy", m_SignalPiMinusenergy);
      status = m_tuple->addItem("SignalKPluspx", m_SignalKPluspx);
      status = m_tuple->addItem("SignalKPluspy", m_SignalKPluspy);
      status = m_tuple->addItem("SignalKPluspz", m_SignalKPluspz);
      status = m_tuple->addItem("SignalKPlusenergy", m_SignalKPlusenergy);
      status = m_tuple->addItem("SignalKMinuspx", m_SignalKMinuspx);
      status = m_tuple->addItem("SignalKMinuspy", m_SignalKMinuspy);
      status = m_tuple->addItem("SignalKMinuspz", m_SignalKMinuspz);
      status = m_tuple->addItem("SignalKMinusenergy", m_SignalKMinusenergy);
      status = m_tuple->addItem("SignalKalmanFitSuccess", m_SignalKalmanFitSuccess);
      status = m_tuple->addItem("SignalKalmanFitChi2", m_SignalKalmanFitChi2);
      status = m_tuple->addItem("SignalPiPluspxKalmanFit", m_SignalPiPluspxKalmanFit);
      status = m_tuple->addItem("SignalPiPluspyKalmanFit", m_SignalPiPluspyKalmanFit);
      status = m_tuple->addItem("SignalPiPluspzKalmanFit", m_SignalPiPluspzKalmanFit);
      status = m_tuple->addItem("SignalPiPlusenergyKalmanFit", m_SignalPiPlusenergyKalmanFit);
      status = m_tuple->addItem("SignalPiMinuspxKalmanFit", m_SignalPiMinuspxKalmanFit);
      status = m_tuple->addItem("SignalPiMinuspyKalmanFit", m_SignalPiMinuspyKalmanFit);
      status = m_tuple->addItem("SignalPiMinuspzKalmanFit", m_SignalPiMinuspzKalmanFit);
      status = m_tuple->addItem("SignalPiMinusenergyKalmanFit", m_SignalPiMinusenergyKalmanFit);
      status = m_tuple->addItem("SignalKPluspxKalmanFit", m_SignalKPluspxKalmanFit);
      status = m_tuple->addItem("SignalKPluspyKalmanFit", m_SignalKPluspyKalmanFit);
      status = m_tuple->addItem("SignalKPluspzKalmanFit", m_SignalKPluspzKalmanFit);
      status = m_tuple->addItem("SignalKPlusenergyKalmanFit", m_SignalKPlusenergyKalmanFit);
      status = m_tuple->addItem("SignalKMinuspxKalmanFit", m_SignalKMinuspxKalmanFit);
      status = m_tuple->addItem("SignalKMinuspyKalmanFit", m_SignalKMinuspyKalmanFit);
      status = m_tuple->addItem("SignalKMinuspzKalmanFit", m_SignalKMinuspzKalmanFit);
      status = m_tuple->addItem("SignalKMinusenergyKalmanFit", m_SignalKMinusenergyKalmanFit);
      status = m_tuple->addItem("SignalMpipi", m_SignalMpipi);
      status = m_tuple->addItem("SignalKSFitSuccess", m_SignalKSFitSuccess);
      status = m_tuple->addItem("SignalKSDecayLengthVeeVertex", m_SignalDecayLengthVeeVertex);
      status = m_tuple->addItem("SignalKSChi2VeeVertex", m_SignalChi2VeeVertex);
      status = m_tuple->addItem("SignalKSMassVeeVertex", m_SignalKSMassVeeVertex);
      status = m_tuple->addItem("SignalKSDecayLengthFit", m_SignalDecayLengthFit);
      status = m_tuple->addItem("SignalKSDecayLengthErrorFit", m_SignalDecayLengthErrorFit);
      status = m_tuple->addItem("SignalKSChi2Fit", m_SignalChi2Fit);
      status = m_tuple->addItem("SignalIsSameDMother", m_SignalIsSameDMother);
      status = m_tuple->addItem("SignalPIDTrue", m_SignalPIDTrue);
      status = m_tuple->addItem("SignalKPlusTrueID", m_SignalKPlusTrueID);
      status = m_tuple->addItem("SignalKMinusTrueID", m_SignalKMinusTrueID);
      status = m_tuple->addItem("SignalPiPlusTrueID", m_SignalPiPlusTrueID);
      status = m_tuple->addItem("SignalPiMinusTrueID", m_SignalPiMinusTrueID);
      status = m_tuple->addItem("TagKSFitSuccess", m_TagKSFitSuccess);
      status = m_tuple->addItem("TagKSDecayLengthVeeVertex", m_TagDecayLengthVeeVertex);
      status = m_tuple->addItem("TagKSChi2VeeVertex", m_TagChi2VeeVertex);
      status = m_tuple->addItem("TagKSMassVeeVertex", m_TagKSMassVeeVertex);
      status = m_tuple->addItem("TagKSDecayLengthFit", m_TagDecayLengthFit);
      status = m_tuple->addItem("TagKSDecayLengthErrorFit", m_TagDecayLengthErrorFit);
      status = m_tuple->addItem("TagKSChi2Fit", m_TagChi2Fit);
      status = m_tuple->addItem("TagKSpx", m_TagKSpx);
      status = m_tuple->addItem("TagKSpy", m_TagKSpy);
      status = m_tuple->addItem("TagKSpz", m_TagKSpz);
      status = m_tuple->addItem("TagKSenergy", m_TagKSenergy);
      status = m_tuple->addItem("TagKSPiPluspx", m_TagKSPiPluspx);
      status = m_tuple->addItem("TagKSPiPluspy", m_TagKSPiPluspy);
      status = m_tuple->addItem("TagKSPiPluspz", m_TagKSPiPluspz);
      status = m_tuple->addItem("TagKSPiPlusenergy", m_TagKSPiPlusenergy);
      status = m_tuple->addItem("TagKSPiMinuspx", m_TagKSPiMinuspx);
      status = m_tuple->addItem("TagKSPiMinuspy", m_TagKSPiMinuspy);
      status = m_tuple->addItem("TagKSPiMinuspz", m_TagKSPiMinuspz);
      status = m_tuple->addItem("TagKSPiMinusenergy", m_TagKSPiMinusenergy);
      status = m_tuple->addItem("TagKPluspx", m_TagKPluspx);
      status = m_tuple->addItem("TagKPluspy", m_TagKPluspy);
      status = m_tuple->addItem("TagKPluspz", m_TagKPluspz);
      status = m_tuple->addItem("TagKPlusenergy", m_TagKPlusenergy);
      status = m_tuple->addItem("TagKMinuspx", m_TagKMinuspx);
      status = m_tuple->addItem("TagKMinuspy", m_TagKMinuspy);
      status = m_tuple->addItem("TagKMinuspz", m_TagKMinuspz);
      status = m_tuple->addItem("TagKMinusenergy", m_TagKMinusenergy);
      status = m_tuple->addItem("TagKalmanFitSuccess", m_TagKalmanFitSuccess);
      status = m_tuple->addItem("TagKalmanFitChi2", m_TagKalmanFitChi2);
      status = m_tuple->addItem("PiPluspxKalmanFit", m_TagKPluspxKalmanFit);
      status = m_tuple->addItem("PiPluspyKalmanFit", m_TagKPluspyKalmanFit);
      status = m_tuple->addItem("PiPluspzKalmanFit", m_TagKPluspzKalmanFit);
      status = m_tuple->addItem("PiPlusenergyKalmanFit", m_TagKPlusenergyKalmanFit);
      status = m_tuple->addItem("PiMinuspxKalmanFit", m_TagKMinuspxKalmanFit);
      status = m_tuple->addItem("PiMinuspyKalmanFit", m_TagKMinuspyKalmanFit);
      status = m_tuple->addItem("PiMinuspzKalmanFit", m_TagKMinuspzKalmanFit);
      status = m_tuple->addItem("PiMinusenergyKalmanFit", m_TagKMinusenergyKalmanFit);
      status = m_tuple->addItem("TagKSpxKalmanFit", m_TagKSpxKalmanFit);
      status = m_tuple->addItem("TagKSpyKalmanFit", m_TagKSpyKalmanFit);
      status = m_tuple->addItem("TagKSpzKalmanFit", m_TagKSpzKalmanFit);
      status = m_tuple->addItem("TagKSenergyKalmanFit", m_TagKSenergyKalmanFit);
      status = m_tuple->addItem("TagIsSameDMother", m_TagIsSameDMother);
      status = m_tuple->addItem("TagPIDTrue", m_TagPIDTrue);
      status = m_tuple->addItem("TagKSPiPlusTrueID", m_TagKSPiPlusTrueID);
      status = m_tuple->addItem("TagKSPiMinusTrueID", m_TagKSPiMinusTrueID);
      status = m_tuple->addItem("TagKPlusTrueID", m_TagKPlusTrueID);
      status = m_tuple->addItem("TagKMinusTrueID", m_TagKMinusTrueID);
      status = m_tuple->addItem("TagKSPiPlusMotherTrueID", m_TagKSPiPlusMotherTrueID);
      status = m_tuple->addItem("TagKSPiMinusMotherTrueID", m_TagKSPiMinusMotherTrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KKpipi vs KSKK Double Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KKpipiVersusKSKKDoubleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KKpipi vs KSKK Double Tag Algorithm" << endreq;
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
  if(DTTool.findDTag(EvtRecDTag::kD0toKKPiPi, EvtRecDTag::kD0toKsKK)) {
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
    m_RunNumber = eventHeader->runNumber();
    m_EventNumber = eventHeader->eventNumber();
    DTagToolIterator DTTool_Signal_iter = DTTool.dtag1();
    DTagToolIterator DTTool_Tag_iter = DTTool.dtag2();
    StatusCode FillTupleStatus = FillTuple(DTTool_Signal_iter, DTTool_Tag_iter, DTTool);
    if(FillTupleStatus != StatusCode::SUCCESS) {
      log << MSG::FATAL << "Assigning KSKK tuple info failed" << endreq;
      return StatusCode::FAILURE;
    }
    m_tuple->write();
  }
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKSKKDoubleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KKpipi vs KSKK Double Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKSKKDoubleTag::FillTuple(DTagToolIterator DTTool_Signal_iter, DTagToolIterator DTTool_Tag_iter, DTagTool &DTTool) {
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
    StatusCode MCStatus = findMCInfo.CalculateMCInfo(MCParticleCol, IMcDecayModeService, 333);
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
  m_SignalDMass = (*DTTool_Signal_iter)->mass();
  m_SignalMBC = (*DTTool_Signal_iter)->mBC();
  m_SignalDeltaE = (*DTTool_Signal_iter)->deltaE();
  m_SignalBeamE = (*DTTool_Signal_iter)->beamE();
  m_TagDMass = (*DTTool_Tag_iter)->mass();
  m_TagMBC = (*DTTool_Tag_iter)->mBC();
  m_TagDeltaE = (*DTTool_Tag_iter)->deltaE();
  m_TagBeamE = (*DTTool_Tag_iter)->beamE();
  m_SignalDpx = (*DTTool_Signal_iter)->p4().x();
  m_SignalDpy = (*DTTool_Signal_iter)->p4().y();
  m_SignalDpz = (*DTTool_Signal_iter)->p4().z();
  m_SignalDenergy = (*DTTool_Signal_iter)->p4().t();
  m_TagDpx = (*DTTool_Tag_iter)->p4().x();
  m_TagDpy = (*DTTool_Tag_iter)->p4().y();
  m_TagDpz = (*DTTool_Tag_iter)->p4().z();
  m_TagDenergy = (*DTTool_Tag_iter)->p4().t();
  FindKKpipiTagInfo findKKpipiTagInfo;
  StatusCode status = findKKpipiTagInfo.CalculateTagInfo(DTTool_Signal_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_SignalKPluspx = findKKpipiTagInfo.GetKPlusP(0);
  m_SignalKPluspy = findKKpipiTagInfo.GetKPlusP(1);
  m_SignalKPluspz = findKKpipiTagInfo.GetKPlusP(2);
  m_SignalKPlusenergy = findKKpipiTagInfo.GetKPlusP(3);
  m_SignalKMinuspx = findKKpipiTagInfo.GetKMinusP(0);
  m_SignalKMinuspy = findKKpipiTagInfo.GetKMinusP(1);
  m_SignalKMinuspz = findKKpipiTagInfo.GetKMinusP(2);
  m_SignalKMinusenergy = findKKpipiTagInfo.GetKMinusP(3);
  m_SignalPiPluspx = findKKpipiTagInfo.GetPiPlusP(0);
  m_SignalPiPluspy = findKKpipiTagInfo.GetPiPlusP(1);
  m_SignalPiPluspz = findKKpipiTagInfo.GetPiPlusP(2);
  m_SignalPiPlusenergy = findKKpipiTagInfo.GetPiPlusP(3);
  m_SignalPiMinuspx = findKKpipiTagInfo.GetPiMinusP(0);
  m_SignalPiMinuspy = findKKpipiTagInfo.GetPiMinusP(1);
  m_SignalPiMinuspz = findKKpipiTagInfo.GetPiMinusP(2);
  m_SignalPiMinusenergy = findKKpipiTagInfo.GetPiMinusP(3);
  m_SignalKalmanFitSuccess = findKKpipiTagInfo.GetKalmanFitSuccess();
  m_SignalKalmanFitChi2 = findKKpipiTagInfo.GetKalmanFitChi2();
  m_SignalKPluspxKalmanFit = findKKpipiTagInfo.GetKPlusPKalmanFit(0);
  m_SignalKPluspyKalmanFit = findKKpipiTagInfo.GetKPlusPKalmanFit(1);
  m_SignalKPluspzKalmanFit = findKKpipiTagInfo.GetKPlusPKalmanFit(2);
  m_SignalKPlusenergyKalmanFit = findKKpipiTagInfo.GetKPlusPKalmanFit(3);
  m_SignalKMinuspxKalmanFit = findKKpipiTagInfo.GetKMinusPKalmanFit(0);
  m_SignalKMinuspyKalmanFit = findKKpipiTagInfo.GetKMinusPKalmanFit(1);
  m_SignalKMinuspzKalmanFit = findKKpipiTagInfo.GetKMinusPKalmanFit(2);
  m_SignalKMinusenergyKalmanFit = findKKpipiTagInfo.GetKMinusPKalmanFit(3);
  m_SignalPiPluspxKalmanFit = findKKpipiTagInfo.GetPiPlusPKalmanFit(0);
  m_SignalPiPluspyKalmanFit = findKKpipiTagInfo.GetPiPlusPKalmanFit(1);
  m_SignalPiPluspzKalmanFit = findKKpipiTagInfo.GetPiPlusPKalmanFit(2);
  m_SignalPiPlusenergyKalmanFit = findKKpipiTagInfo.GetPiPlusPKalmanFit(3);
  m_SignalPiMinuspxKalmanFit = findKKpipiTagInfo.GetPiMinusPKalmanFit(0);
  m_SignalPiMinuspyKalmanFit = findKKpipiTagInfo.GetPiMinusPKalmanFit(1);
  m_SignalPiMinuspzKalmanFit = findKKpipiTagInfo.GetPiMinusPKalmanFit(2);
  m_SignalPiMinusenergyKalmanFit = findKKpipiTagInfo.GetPiMinusPKalmanFit(3);
  m_SignalMpipi = findKKpipiTagInfo.GetMpipi();
  m_SignalKSFitSuccess = findKKpipiTagInfo.GetKSFitSuccess();
  m_SignalDecayLengthVeeVertex = findKKpipiTagInfo.GetDecayLengthVeeVertex();
  m_SignalChi2VeeVertex = findKKpipiTagInfo.GetChi2VeeVertex();
  m_SignalKSMassVeeVertex = findKKpipiTagInfo.GetKSMassVeeVertex();
  m_SignalDecayLengthFit = findKKpipiTagInfo.GetDecayLengthFit();
  m_SignalDecayLengthErrorFit = findKKpipiTagInfo.GetDecayLengthErrorFit();
  m_SignalChi2Fit = findKKpipiTagInfo.GetChi2Fit();
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKKpipiTagInfo.GetDaughterTrackID(), 4, this);
    m_SignalIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {321, -321, 211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_SignalPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_SignalKPlusTrueID = ReconstructedPID[0];
    m_SignalKMinusTrueID = ReconstructedPID[1];
    m_SignalPiPlusTrueID = ReconstructedPID[2];
    m_SignalPiMinusTrueID = ReconstructedPID[3];
  }
  FindKShhTagInfo findKSKKTagInfo("KSKK");
  findKSKKTagInfo.CalculateTagInfo(DTTool_Tag_iter, DTTool);
  m_TagKSFitSuccess = findKSKKTagInfo.GetKSFitSuccess();
  m_TagDecayLengthVeeVertex = findKSKKTagInfo.GetDecayLengthVeeVertex();
  m_TagChi2VeeVertex = findKSKKTagInfo.GetChi2VeeVertex();
  m_TagKSMassVeeVertex = findKSKKTagInfo.GetKSMassVeeVertex();
  m_TagDecayLengthFit = findKSKKTagInfo.GetDecayLengthFit();
  m_TagDecayLengthErrorFit = findKSKKTagInfo.GetDecayLengthErrorFit();
  m_TagChi2Fit = findKSKKTagInfo.GetChi2Fit();
  m_TagKSPiPluspx = findKSKKTagInfo.GetKSPiPlusP(0);
  m_TagKSPiPluspy = findKSKKTagInfo.GetKSPiPlusP(1);
  m_TagKSPiPluspz = findKSKKTagInfo.GetKSPiPlusP(2);
  m_TagKSPiPlusenergy = findKSKKTagInfo.GetKSPiPlusP(3);
  m_TagKSPiMinuspx = findKSKKTagInfo.GetKSPiMinusP(0);
  m_TagKSPiMinuspy = findKSKKTagInfo.GetKSPiMinusP(1);
  m_TagKSPiMinuspz = findKSKKTagInfo.GetKSPiMinusP(2);
  m_TagKSPiMinusenergy = findKSKKTagInfo.GetKSPiMinusP(3);
  m_TagKSpx = findKSKKTagInfo.GetKShortP(0);
  m_TagKSpy = findKSKKTagInfo.GetKShortP(1);
  m_TagKSpz = findKSKKTagInfo.GetKShortP(2);
  m_TagKSenergy = findKSKKTagInfo.GetKShortP(3);
  m_TagKPluspx = findKSKKTagInfo.GethPlusP(0);
  m_TagKPluspy = findKSKKTagInfo.GethPlusP(1);
  m_TagKPluspz = findKSKKTagInfo.GethPlusP(2);
  m_TagKPlusenergy = findKSKKTagInfo.GethPlusP(3);
  m_TagKMinuspx = findKSKKTagInfo.GethMinusP(0);
  m_TagKMinuspy = findKSKKTagInfo.GethMinusP(1);
  m_TagKMinuspz = findKSKKTagInfo.GethMinusP(2);
  m_TagKMinusenergy = findKSKKTagInfo.GethMinusP(3);
  m_TagKalmanFitSuccess = findKSKKTagInfo.GetKalmanFitSuccess();
  m_TagKalmanFitChi2 = findKSKKTagInfo.GetKalmanFitChi2();
  m_TagKSpxKalmanFit = findKSKKTagInfo.GetKShortPKalmanFit(0);
  m_TagKSpyKalmanFit = findKSKKTagInfo.GetKShortPKalmanFit(1);
  m_TagKSpzKalmanFit = findKSKKTagInfo.GetKShortPKalmanFit(2);
  m_TagKSenergyKalmanFit = findKSKKTagInfo.GetKShortPKalmanFit(3);
  m_TagKPluspxKalmanFit = findKSKKTagInfo.GethPlusPKalmanFit(0);
  m_TagKPluspyKalmanFit = findKSKKTagInfo.GethPlusPKalmanFit(1);
  m_TagKPluspzKalmanFit = findKSKKTagInfo.GethPlusPKalmanFit(2);
  m_TagKPlusenergyKalmanFit = findKSKKTagInfo.GethPlusPKalmanFit(3);
  m_TagKMinuspxKalmanFit = findKSKKTagInfo.GethMinusPKalmanFit(0);
  m_TagKMinuspyKalmanFit = findKSKKTagInfo.GethMinusPKalmanFit(1);
  m_TagKMinuspzKalmanFit = findKSKKTagInfo.GethMinusPKalmanFit(2);
  m_TagKMinusenergyKalmanFit = findKSKKTagInfo.GethMinusPKalmanFit(3);
  if(m_RunNumber < 0) {
    std::vector<int> DaughterTrackIDs = findKSKKTagInfo.GetDaughterTrackID();
    PIDTruth PID_Truth(DaughterTrackIDs, 4, this);
    m_TagIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {211, -211, 321, -321};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_TagPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_TagKPlusTrueID = ReconstructedPID[0];
    m_TagKMinusTrueID = ReconstructedPID[1];
    m_TagKSPiPlusTrueID = ReconstructedPID[2];
    m_TagKSPiMinusTrueID = ReconstructedPID[3];
    m_TagKSPiPlusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[0], true);
    m_TagKSPiMinusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[1], true);
  }
  return StatusCode::SUCCESS;
}
