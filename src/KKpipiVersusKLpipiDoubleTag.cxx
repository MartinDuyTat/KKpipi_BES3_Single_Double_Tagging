// Martin Duy Tat 2nd January 2022

// KKpipi
#include "KKpipi/KKpipiVersusKLpipiDoubleTag.h"
#include "KKpipi/FindKKpipiTagInfo.h"
#include "KKpipi/FindMCInfo.h"
#include "KKpipi/FindKL.h"
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
// ROOT
#include "TMath.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McTruth/McParticle.h"
// STL
#include<vector>
#include<string>

DECLARE_COMPONENT(KKpipiVersusKLpipiDoubleTag)

KKpipiVersusKLpipiDoubleTag::KKpipiVersusKLpipiDoubleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKpipiVersusKLpipiDoubleTag::~KKpipiVersusKLpipiDoubleTag() {
}

StatusCode KKpipiVersusKLpipiDoubleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KKpipi vs KLpipi Double Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KLpipiDoubleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KLpipiDoubleTag", CLID_ColumnWiseTuple, "Double tagged D->KKpipi vs D->KLpipi events");
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
      status = m_tuple->addItem("SignalDpx", m_SignalDpx);
      status = m_tuple->addItem("SignalDpy", m_SignalDpy);
      status = m_tuple->addItem("SignalDpz", m_SignalDpz);
      status = m_tuple->addItem("SignalDenergy", m_SignalDenergy);
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
      status = m_tuple->addItem("SignalDaughters", m_SignalDaughters, 0, 100);
      status = m_tuple->addIndexedItem("SignalDOrigin", m_SignalDaughters, m_SignalDOrigin);
      status = m_tuple->addItem("TagKLpx", m_TagKLpx);
      status = m_tuple->addItem("TagKLpy", m_TagKLpy);
      status = m_tuple->addItem("TagKLpz", m_TagKLpz);
      status = m_tuple->addItem("TagKLenergy", m_TagKLenergy);
      status = m_tuple->addItem("TagPiPluspx", m_TagPiPluspx);
      status = m_tuple->addItem("TagPiPluspy", m_TagPiPluspy);
      status = m_tuple->addItem("TagPiPluspz", m_TagPiPluspz);
      status = m_tuple->addItem("TagPiPlusenergy", m_TagPiPlusenergy);
      status = m_tuple->addItem("TagPiMinuspx", m_TagPiMinuspx);
      status = m_tuple->addItem("TagPiMinuspy", m_TagPiMinuspy);
      status = m_tuple->addItem("TagPiMinuspz", m_TagPiMinuspz);
      status = m_tuple->addItem("TagPiMinusenergy", m_TagPiMinusenergy);
      status = m_tuple->addItem("TagKalmanFitSuccess", m_TagKalmanFitSuccess);
      status = m_tuple->addItem("TagKalmanFitChi2", m_TagKalmanFitChi2);
      status = m_tuple->addItem("TagKLpxKalmanFit", m_TagKLpxKalmanFit);
      status = m_tuple->addItem("TagKLpyKalmanFit", m_TagKLpyKalmanFit);
      status = m_tuple->addItem("TagKLpzKalmanFit", m_TagKLpzKalmanFit);
      status = m_tuple->addItem("TagKLenergyKalmanFit", m_TagKLenergyKalmanFit);
      status = m_tuple->addItem("TagPiPluspxKalmanFit", m_TagPiPluspxKalmanFit);
      status = m_tuple->addItem("TagPiPluspyKalmanFit", m_TagPiPluspyKalmanFit);
      status = m_tuple->addItem("TagPiPluspzKalmanFit", m_TagPiPluspzKalmanFit);
      status = m_tuple->addItem("TagPiPlusenergyKalmanFit", m_TagPiPlusenergyKalmanFit);
      status = m_tuple->addItem("TagPiMinuspxKalmanFit", m_TagPiMinuspxKalmanFit);
      status = m_tuple->addItem("TagPiMinuspyKalmanFit", m_TagPiMinuspyKalmanFit);
      status = m_tuple->addItem("TagPiMinuspzKalmanFit", m_TagPiMinuspzKalmanFit);
      status = m_tuple->addItem("TagPiMinusenergyKalmanFit", m_TagPiMinusenergyKalmanFit);
      status = m_tuple->addItem("TagNumberGamma", m_TagNumberGamma, 0, 100);
      status = m_tuple->addIndexedItem("TagPhotonEnergy", m_TagNumberGamma, m_TagPhotonEnergy);
      status = m_tuple->addIndexedItem("TagPhotonPx", m_TagNumberGamma, m_TagPhotonPx);
      status = m_tuple->addIndexedItem("TagPhotonPy", m_TagNumberGamma, m_TagPhotonPy);
      status = m_tuple->addIndexedItem("TagPhotonPz", m_TagNumberGamma, m_TagPhotonPz);
      status = m_tuple->addIndexedItem("TagPhotonAngleSeparation", m_TagNumberGamma, m_TagPhotonAngleSeparation);
      status = m_tuple->addIndexedItem("TagPhotonThetaSeparation", m_TagNumberGamma, m_TagPhotonThetaSeparation);
      status = m_tuple->addIndexedItem("TagPhotonPhiSeparation", m_TagNumberGamma, m_TagPhotonPhiSeparation);
      status = m_tuple->addItem("TagMissingMass2", m_TagMissingMass2);
      status = m_tuple->addItem("TagIsSameDMother", m_TagIsSameDMother);
      status = m_tuple->addItem("TagPIDTrue", m_TagPIDTrue);
      status = m_tuple->addItem("TagPiPlusTrueID", m_TagPiPlusTrueID);
      status = m_tuple->addItem("TagPiMinusTrueID", m_TagPiMinusTrueID);
      status = m_tuple->addItem("TagDaughters", m_TagDaughters, 0, 100);
      status = m_tuple->addIndexedItem("TagDOrigin", m_TagDaughters, m_TagDOrigin);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KKpipi vs KLpipi Double Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KKpipiVersusKLpipiDoubleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KKpipi vs KLpipi Double Tag Algorithm" << endreq;
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
    DTagToolIterator DTTool_Signal_iter = DTTool.stag();
    StatusCode FillTupleStatus = FillTuple(DTTool_Signal_iter, DTTool);
    if(FillTupleStatus == StatusCode::RECOVERABLE) {
      m_RunNumber = 0;
      m_EventNumber = 0;
      return StatusCode::SUCCESS;
    } else if(FillTupleStatus == StatusCode::FAILURE) {
      log << MSG::FATAL << "Assigning KL tuple info failed" << endreq;
      return StatusCode::FAILURE;
    } else {
      m_tuple->write();
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKLpipiDoubleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KKpipi vs KLpipi Double Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKLpipiDoubleTag::FillTuple(DTagToolIterator DTTool_Signal_iter, DTagTool &DTTool) {
  // First check if there are any KLpipi candidates, otherwise no point in saving all the other stuff
  FindKL findKL;
  StatusCode FoundKL = findKL.findKL(DTTool_Signal_iter, DTTool);
  if(FoundKL == StatusCode::FAILURE) {
    return StatusCode::RECOVERABLE;
  }
  if(!findKL.FoundKLpipiTag()) {
    return StatusCode::RECOVERABLE;
  }
  // Save MC truth information
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
  // Save information about signal side
  m_SignalDMass = (*DTTool_Signal_iter)->mass();
  m_SignalMBC = (*DTTool_Signal_iter)->mBC();
  m_SignalDeltaE = (*DTTool_Signal_iter)->deltaE();
  m_SignalBeamE = (*DTTool_Signal_iter)->beamE();
  m_SignalDpx = (*DTTool_Signal_iter)->p4().x();
  m_SignalDpy = (*DTTool_Signal_iter)->p4().y();
  m_SignalDpz = (*DTTool_Signal_iter)->p4().z();
  m_SignalDenergy = (*DTTool_Signal_iter)->p4().t();
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
    m_SignalDaughters = 4;
    for(int i = 0; i < 4; i++) {
      m_SignalDOrigin[i] = PID_Truth.FindDOrigin(findKKpipiTagInfo.GetDaughterTrackID()[i], true);
    }
  }
  // Save information about tag side
  m_TagKLpx = findKL.GetKLongP(0);
  m_TagKLpy = findKL.GetKLongP(1);
  m_TagKLpz = findKL.GetKLongP(2);
  m_TagKLenergy = findKL.GetKLongP(3);
  m_TagPiPluspx = findKL.GethPlusP(0);
  m_TagPiPluspy = findKL.GethPlusP(1);
  m_TagPiPluspz = findKL.GethPlusP(2);
  m_TagPiPlusenergy = findKL.GethPlusP(3);
  m_TagPiMinuspx = findKL.GethMinusP(0);
  m_TagPiMinuspy = findKL.GethMinusP(1);
  m_TagPiMinuspz = findKL.GethMinusP(2);
  m_TagPiMinusenergy = findKL.GethMinusP(3);
  m_TagKalmanFitSuccess = findKL.GetKalmanFitSuccess();
  m_TagKalmanFitChi2 = findKL.GetKalmanFitChi2();
  m_TagKLpxKalmanFit = findKL.GetKLongPKalmanFit(0);
  m_TagKLpyKalmanFit = findKL.GetKLongPKalmanFit(1);
  m_TagKLpzKalmanFit = findKL.GetKLongPKalmanFit(2);
  m_TagKLenergyKalmanFit = findKL.GetKLongPKalmanFit(3);
  m_TagPiPluspxKalmanFit = findKL.GethPlusPKalmanFit(0);
  m_TagPiPluspyKalmanFit = findKL.GethPlusPKalmanFit(1);
  m_TagPiPluspzKalmanFit = findKL.GethPlusPKalmanFit(2);
  m_TagPiPlusenergyKalmanFit = findKL.GethPlusPKalmanFit(3);
  m_TagPiMinuspxKalmanFit = findKL.GethMinusPKalmanFit(0);
  m_TagPiMinuspyKalmanFit = findKL.GethMinusPKalmanFit(1);
  m_TagPiMinuspzKalmanFit = findKL.GethMinusPKalmanFit(2);
  m_TagPiMinusenergyKalmanFit = findKL.GethMinusPKalmanFit(3);
  m_TagMissingMass2 = findKL.GetMissingMass2();
  m_TagNumberGamma = findKL.GetNumberGamma();
  for(int j = 0; j < m_TagNumberGamma; j++) {
    m_TagPhotonPx[j] = findKL.GetPhotonP(0, j);
    m_TagPhotonPy[j] = findKL.GetPhotonP(1, j);
    m_TagPhotonPz[j] = findKL.GetPhotonP(2, j);
    m_TagPhotonEnergy[j] = findKL.GetPhotonP(3, j);
    m_TagPhotonAngleSeparation[j] = findKL.GetPhotonAngleSeparation(j);
    m_TagPhotonThetaSeparation[j] = findKL.GetPhotonThetaSeparation(j);
    m_TagPhotonPhiSeparation[j] = findKL.GetPhotonPhiSeparation(j);
  }
  if(m_RunNumber < 0) {
    std::vector<int> DaughterTrackIDs = findKL.GetDaughterTrackID();
    PIDTruth PID_Truth(DaughterTrackIDs, 2, this);
    m_TagIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[2] = {211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 2);
    m_TagPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_TagPiPlusTrueID = ReconstructedPID[0];
    m_TagPiMinusTrueID = ReconstructedPID[1];
    m_TagDaughters = 2;
    m_TagDOrigin[0] = PID_Truth.FindDOrigin(DaughterTrackIDs[0], true);
    m_TagDOrigin[1] = PID_Truth.FindDOrigin(DaughterTrackIDs[1], true);
  }
  return StatusCode::SUCCESS;
}
