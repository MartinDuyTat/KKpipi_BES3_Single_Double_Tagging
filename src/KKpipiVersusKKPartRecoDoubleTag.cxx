// Martin Duy Tat 19th January 2023

// KKpipi
#include "KKpipi/KKpipiVersusKKPartRecoDoubleTag.h"
#include "KKpipi/FindKKpipiVersusKKPartRecoTagInfo.h"
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

KKpipiVersusKKPartRecoDoubleTag::KKpipiVersusKKPartRecoDoubleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKpipiVersusKKPartRecoDoubleTag::~KKpipiVersusKKPartRecoDoubleTag() {
}

StatusCode KKpipiVersusKKPartRecoDoubleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KKpipi vs KK PartReco Double Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KKPartRecoDoubleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KKPartRecoDoubleTag", CLID_ColumnWiseTuple, "Double tagged D->KKpipi vs D->KK events with missing K in KKpipi");
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
      status = m_tuple->addItem("TagDMass", m_TagDMass);
      status = m_tuple->addItem("TagMBC", m_TagMBC);
      status = m_tuple->addItem("TagDeltaE", m_TagDeltaE);
      status = m_tuple->addItem("TagBeamE", m_TagBeamE);
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
      status = m_tuple->addItem("SignalMissingMass2", m_SignalMissingMass2);
      status = m_tuple->addItem("SignalRecKCharge", m_SignalRecKCharge);
      status = m_tuple->addItem("SignalNumberPi0", m_SignalNumberPi0);
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
      status = m_tuple->addItem("SignalKTrueID", m_SignalKTrueID);
      status = m_tuple->addItem("SignalPiPlusTrueID", m_SignalPiPlusTrueID);
      status = m_tuple->addItem("SignalPiMinusTrueID", m_SignalPiMinusTrueID);
      status = m_tuple->addItem("TagKPluspx", m_TagKPluspx);
      status = m_tuple->addItem("TagKPluspy", m_TagKPluspy);
      status = m_tuple->addItem("TagKPluspz", m_TagKPluspz);
      status = m_tuple->addItem("TagKPlusenergy", m_TagKPlusenergy);
      status = m_tuple->addItem("TagKMinuspx", m_TagKMinuspx);
      status = m_tuple->addItem("TagKMinuspy", m_TagKMinuspy);
      status = m_tuple->addItem("TagKMinuspz", m_TagKMinuspz);
      status = m_tuple->addItem("TagKMinusenergy", m_TagKMinusenergy);
      status = m_tuple->addItem("TagIsSameDMother", m_TagIsSameDMother);
      status = m_tuple->addItem("TagPIDTrue", m_TagPIDTrue);
      status = m_tuple->addItem("TagKPlusTrueID", m_TagKPlusTrueID);
      status = m_tuple->addItem("TagKMinusTrueID", m_TagKMinusTrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KKpipi vs KK Part Reco Double Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KKpipiVersusKKPartRecoDoubleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KKpipi vs KK PartReco Double Tag Algorithm" << endreq;
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
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
    m_RunNumber = eventHeader->runNumber();
    m_EventNumber = eventHeader->eventNumber();
    DTagToolIterator DTTool_Tag_iter = DTTool.stag();
    StatusCode FillTupleStatus = FillTuple(DTTool_Tag_iter, DTTool);
    if(FillTupleStatus != StatusCode::SUCCESS) {
      return StatusCode::SUCCESS;
    }
    m_tuple->write();
  }
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKKPartRecoDoubleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KKpipi vs KK Double Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKKPartRecoDoubleTag::FillTuple(DTagToolIterator DTTool_Tag_iter, DTagTool &DTTool) {
  FindKKpipiVersusKKPartRecoTagInfo findKKpipiVersusKKPartRecoTagInfo;
  StatusCode status = findKKpipiVersusKKPartRecoTagInfo.CalculateTagInfo(DTTool_Tag_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
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
  m_TagDMass = (*DTTool_Tag_iter)->mass();
  m_TagMBC = (*DTTool_Tag_iter)->mBC();
  m_TagDeltaE = (*DTTool_Tag_iter)->deltaE();
  m_TagBeamE = (*DTTool_Tag_iter)->beamE();
  m_TagDpx = (*DTTool_Tag_iter)->p4().x();
  m_TagDpy = (*DTTool_Tag_iter)->p4().y();
  m_TagDpz = (*DTTool_Tag_iter)->p4().z();
  m_TagDenergy = (*DTTool_Tag_iter)->p4().t();
  m_SignalKPluspx = findKKpipiVersusKKPartRecoTagInfo.GetKPlusP(0);
  m_SignalKPluspy = findKKpipiVersusKKPartRecoTagInfo.GetKPlusP(1);
  m_SignalKPluspz = findKKpipiVersusKKPartRecoTagInfo.GetKPlusP(2);
  m_SignalKPlusenergy = findKKpipiVersusKKPartRecoTagInfo.GetKPlusP(3);
  m_SignalKMinuspx = findKKpipiVersusKKPartRecoTagInfo.GetKMinusP(0);
  m_SignalKMinuspy = findKKpipiVersusKKPartRecoTagInfo.GetKMinusP(1);
  m_SignalKMinuspz = findKKpipiVersusKKPartRecoTagInfo.GetKMinusP(2);
  m_SignalKMinusenergy = findKKpipiVersusKKPartRecoTagInfo.GetKMinusP(3);
  m_SignalPiPluspx = findKKpipiVersusKKPartRecoTagInfo.GetPiPlusP(0);
  m_SignalPiPluspy = findKKpipiVersusKKPartRecoTagInfo.GetPiPlusP(1);
  m_SignalPiPluspz = findKKpipiVersusKKPartRecoTagInfo.GetPiPlusP(2);
  m_SignalPiPlusenergy = findKKpipiVersusKKPartRecoTagInfo.GetPiPlusP(3);
  m_SignalPiMinuspx = findKKpipiVersusKKPartRecoTagInfo.GetPiMinusP(0);
  m_SignalPiMinuspy = findKKpipiVersusKKPartRecoTagInfo.GetPiMinusP(1);
  m_SignalPiMinuspz = findKKpipiVersusKKPartRecoTagInfo.GetPiMinusP(2);
  m_SignalPiMinusenergy = findKKpipiVersusKKPartRecoTagInfo.GetPiMinusP(3);
  m_SignalKalmanFitSuccess = findKKpipiVersusKKPartRecoTagInfo.GetKalmanFitSuccess();
  m_SignalKalmanFitChi2 = findKKpipiVersusKKPartRecoTagInfo.GetKalmanFitChi2();
  m_SignalKPluspxKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetKPlusPKalmanFit(0);
  m_SignalKPluspyKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetKPlusPKalmanFit(1);
  m_SignalKPluspzKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetKPlusPKalmanFit(2);
  m_SignalKPlusenergyKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetKPlusPKalmanFit(3);
  m_SignalKMinuspxKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetKMinusPKalmanFit(0);
  m_SignalKMinuspyKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetKMinusPKalmanFit(1);
  m_SignalKMinuspzKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetKMinusPKalmanFit(2);
  m_SignalKMinusenergyKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetKMinusPKalmanFit(3);
  m_SignalPiPluspxKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetPiPlusPKalmanFit(0);
  m_SignalPiPluspyKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetPiPlusPKalmanFit(1);
  m_SignalPiPluspzKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetPiPlusPKalmanFit(2);
  m_SignalPiPlusenergyKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetPiPlusPKalmanFit(3);
  m_SignalPiMinuspxKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetPiMinusPKalmanFit(0);
  m_SignalPiMinuspyKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetPiMinusPKalmanFit(1);
  m_SignalPiMinuspzKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetPiMinusPKalmanFit(2);
  m_SignalPiMinusenergyKalmanFit = findKKpipiVersusKKPartRecoTagInfo.GetPiMinusPKalmanFit(3);
  m_SignalMpipi = findKKpipiVersusKKPartRecoTagInfo.GetMpipi_KKpipi();
  m_SignalMissingMass2 = findKKpipiVersusKKPartRecoTagInfo.GetMissingMass2();
  m_SignalRecKCharge = findKKpipiVersusKKPartRecoTagInfo.GetRecKCharge();
  m_SignalNumberPi0 = findKKpipiVersusKKPartRecoTagInfo.GetNumberPi0();
  m_SignalKSFitSuccess = findKKpipiVersusKKPartRecoTagInfo.GetKSFitSuccess_KKpipi();
  m_SignalDecayLengthVeeVertex = findKKpipiVersusKKPartRecoTagInfo.GetDecayLengthVeeVertex_KKpipi();
  m_SignalChi2VeeVertex = findKKpipiVersusKKPartRecoTagInfo.GetChi2VeeVertex_KKpipi();
  m_SignalKSMassVeeVertex = findKKpipiVersusKKPartRecoTagInfo.GetKSMassVeeVertex_KKpipi();
  m_SignalDecayLengthFit = findKKpipiVersusKKPartRecoTagInfo.GetDecayLengthFit_KKpipi();
  m_SignalDecayLengthErrorFit = findKKpipiVersusKKPartRecoTagInfo.GetDecayLengthErrorFit_KKpipi();
  m_SignalChi2Fit = findKKpipiVersusKKPartRecoTagInfo.GetChi2Fit_KKpipi();
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKKpipiVersusKKPartRecoTagInfo.GetDaughterTrackID_KKpipi(), 3, this);
    m_SignalIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[3] = {m_SignalRecKCharge*321, 211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 3);
    m_SignalPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_SignalKTrueID = ReconstructedPID[0];
    m_SignalPiPlusTrueID = ReconstructedPID[1];
    m_SignalPiMinusTrueID = ReconstructedPID[2];
  }
  m_TagKPluspx = findKKpipiVersusKKPartRecoTagInfo.GetKKTagInfo().GethPlusP(0);
  m_TagKPluspy = findKKpipiVersusKKPartRecoTagInfo.GetKKTagInfo().GethPlusP(1);
  m_TagKPluspz = findKKpipiVersusKKPartRecoTagInfo.GetKKTagInfo().GethPlusP(2);
  m_TagKPlusenergy = findKKpipiVersusKKPartRecoTagInfo.GetKKTagInfo().GethPlusP(3);
  m_TagKMinuspx = findKKpipiVersusKKPartRecoTagInfo.GetKKTagInfo().GethMinusP(0);
  m_TagKMinuspy = findKKpipiVersusKKPartRecoTagInfo.GetKKTagInfo().GethMinusP(1);
  m_TagKMinuspz = findKKpipiVersusKKPartRecoTagInfo.GetKKTagInfo().GethMinusP(2);
  m_TagKMinusenergy = findKKpipiVersusKKPartRecoTagInfo.GetKKTagInfo().GethMinusP(3);
  if(m_RunNumber < 0) {
    std::vector<int> DaughterTrackIDs = findKKpipiVersusKKPartRecoTagInfo.GetKKTagInfo().GetDaughterTrackID();
    PIDTruth PID_Truth(DaughterTrackIDs, 2, this);
    m_TagIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {321, -321};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 2);
    m_TagPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_TagKPlusTrueID = ReconstructedPID[0];
    m_TagKMinusTrueID = ReconstructedPID[1];
  }
  return StatusCode::SUCCESS;
}
