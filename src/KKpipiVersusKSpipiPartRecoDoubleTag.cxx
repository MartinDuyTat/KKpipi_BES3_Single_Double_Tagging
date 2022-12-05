// Martin Duy Tat 2nd February 2022

// KKpipi
#include "KKpipi/KKpipiVersusKSpipiPartRecoDoubleTag.h"
#include "KKpipi/FindKKpipiVersusKSpipiPartRecoTagInfo.h"
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

KKpipiVersusKSpipiPartRecoDoubleTag::KKpipiVersusKSpipiPartRecoDoubleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKpipiVersusKSpipiPartRecoDoubleTag::~KKpipiVersusKSpipiPartRecoDoubleTag() {
}

StatusCode KKpipiVersusKSpipiPartRecoDoubleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KKpipi vs KSpipi PartReco Double Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KSpipiPartRecoDoubleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KSpipiPartRecoDoubleTag", CLID_ColumnWiseTuple, "Double tagged D->KKpipi vs D->KSpipi events with missing K in KKpipi");
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
      status = m_tuple->addItem("TagNumberGamma", m_TagNumberGamma, 0, 100);
      status = m_tuple->addIndexedItem("TagExtraShowerEnergy", m_TagNumberGamma, m_TagExtraShowerEnergy);
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
      status = m_tuple->addItem("TagPiPluspx", m_TagPiPluspx);
      status = m_tuple->addItem("TagPiPluspy", m_TagPiPluspy);
      status = m_tuple->addItem("TagPiPluspz", m_TagPiPluspz);
      status = m_tuple->addItem("TagPiPlusenergy", m_TagPiPlusenergy);
      status = m_tuple->addItem("TagPiMinuspx", m_TagPiMinuspx);
      status = m_tuple->addItem("TagPiMinuspy", m_TagPiMinuspy);
      status = m_tuple->addItem("TagPiMinuspz", m_TagPiMinuspz);
      status = m_tuple->addItem("TagPiMinusenergy", m_TagPiMinusenergy);
      status = m_tuple->addItem("TagPiPluspxKalmanFit", m_TagPiPluspxKalmanFit);
      status = m_tuple->addItem("TagPiPluspyKalmanFit", m_TagPiPluspyKalmanFit);
      status = m_tuple->addItem("TagPiPluspzKalmanFit", m_TagPiPluspzKalmanFit);
      status = m_tuple->addItem("TagPiPlusenergyKalmanFit", m_TagPiPlusenergyKalmanFit);
      status = m_tuple->addItem("TagPiMinuspxKalmanFit", m_TagPiMinuspxKalmanFit);
      status = m_tuple->addItem("TagPiMinuspyKalmanFit", m_TagPiMinuspyKalmanFit);
      status = m_tuple->addItem("TagPiMinuspzKalmanFit", m_TagPiMinuspzKalmanFit);
      status = m_tuple->addItem("TagPiMinusenergyKalmanFit", m_TagPiMinusenergyKalmanFit);
      status = m_tuple->addItem("TagKSpxKalmanFit", m_TagKSpxKalmanFit);
      status = m_tuple->addItem("TagKSpyKalmanFit", m_TagKSpyKalmanFit);
      status = m_tuple->addItem("TagKSpzKalmanFit", m_TagKSpzKalmanFit);
      status = m_tuple->addItem("TagKSenergyKalmanFit", m_TagKSenergyKalmanFit);
      status = m_tuple->addItem("TagIsSameDMother", m_TagIsSameDMother);
      status = m_tuple->addItem("TagPIDTrue", m_TagPIDTrue);
      status = m_tuple->addItem("TagKSPiPlusTrueID", m_TagKSPiPlusTrueID);
      status = m_tuple->addItem("TagKSPiMinusTrueID", m_TagKSPiMinusTrueID);
      status = m_tuple->addItem("TagPiPlusTrueID", m_TagPiPlusTrueID);
      status = m_tuple->addItem("TagPiMinusTrueID", m_TagPiMinusTrueID);
      status = m_tuple->addItem("TagKSPiPlusMotherTrueID", m_TagKSPiPlusMotherTrueID);
      status = m_tuple->addItem("TagKSPiMinusMotherTrueID", m_TagKSPiMinusMotherTrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KKpipi vs KSpipi Part Reco Double Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KKpipiVersusKSpipiPartRecoDoubleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KKpipi vs KSpipi PartReco Double Tag Algorithm" << endreq;
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
    DTagToolIterator DTTool_Tag_iter = DTTool.stag();
    StatusCode FillTupleStatus = FillTuple(DTTool_Tag_iter, DTTool);
    if(FillTupleStatus != StatusCode::SUCCESS) {
      return StatusCode::SUCCESS;
    }
    m_tuple->write();
  }
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKSpipiPartRecoDoubleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KKpipi vs KSpipi Double Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKSpipiPartRecoDoubleTag::FillTuple(DTagToolIterator DTTool_Tag_iter, DTagTool &DTTool) {
  FindKKpipiVersusKSpipiPartRecoTagInfo findKKpipiVersusKSpipiPartRecoTagInfo;
  StatusCode status = findKKpipiVersusKSpipiPartRecoTagInfo.CalculateTagInfo(DTTool_Tag_iter, DTTool);
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
  m_SignalKPluspx = findKKpipiVersusKSpipiPartRecoTagInfo.GetKPlusP(0);
  m_SignalKPluspy = findKKpipiVersusKSpipiPartRecoTagInfo.GetKPlusP(1);
  m_SignalKPluspz = findKKpipiVersusKSpipiPartRecoTagInfo.GetKPlusP(2);
  m_SignalKPlusenergy = findKKpipiVersusKSpipiPartRecoTagInfo.GetKPlusP(3);
  m_SignalKMinuspx = findKKpipiVersusKSpipiPartRecoTagInfo.GetKMinusP(0);
  m_SignalKMinuspy = findKKpipiVersusKSpipiPartRecoTagInfo.GetKMinusP(1);
  m_SignalKMinuspz = findKKpipiVersusKSpipiPartRecoTagInfo.GetKMinusP(2);
  m_SignalKMinusenergy = findKKpipiVersusKSpipiPartRecoTagInfo.GetKMinusP(3);
  m_SignalPiPluspx = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiPlusP(0);
  m_SignalPiPluspy = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiPlusP(1);
  m_SignalPiPluspz = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiPlusP(2);
  m_SignalPiPlusenergy = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiPlusP(3);
  m_SignalPiMinuspx = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiMinusP(0);
  m_SignalPiMinuspy = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiMinusP(1);
  m_SignalPiMinuspz = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiMinusP(2);
  m_SignalPiMinusenergy = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiMinusP(3);
  m_SignalKalmanFitSuccess = findKKpipiVersusKSpipiPartRecoTagInfo.GetKalmanFitSuccess();
  m_SignalKalmanFitChi2 = findKKpipiVersusKSpipiPartRecoTagInfo.GetKalmanFitChi2();
  m_SignalKPluspxKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKPlusPKalmanFit(0);
  m_SignalKPluspyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKPlusPKalmanFit(1);
  m_SignalKPluspzKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKPlusPKalmanFit(2);
  m_SignalKPlusenergyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKPlusPKalmanFit(3);
  m_SignalKMinuspxKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKMinusPKalmanFit(0);
  m_SignalKMinuspyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKMinusPKalmanFit(1);
  m_SignalKMinuspzKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKMinusPKalmanFit(2);
  m_SignalKMinusenergyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKMinusPKalmanFit(3);
  m_SignalPiPluspxKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiPlusPKalmanFit(0);
  m_SignalPiPluspyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiPlusPKalmanFit(1);
  m_SignalPiPluspzKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiPlusPKalmanFit(2);
  m_SignalPiPlusenergyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiPlusPKalmanFit(3);
  m_SignalPiMinuspxKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiMinusPKalmanFit(0);
  m_SignalPiMinuspyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiMinusPKalmanFit(1);
  m_SignalPiMinuspzKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiMinusPKalmanFit(2);
  m_SignalPiMinusenergyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetPiMinusPKalmanFit(3);
  m_SignalMpipi = findKKpipiVersusKSpipiPartRecoTagInfo.GetMpipi_KKpipi();
  m_SignalMissingMass2 = findKKpipiVersusKSpipiPartRecoTagInfo.GetMissingMass2();
  m_SignalRecKCharge = findKKpipiVersusKSpipiPartRecoTagInfo.GetRecKCharge();
  m_SignalNumberPi0 = findKKpipiVersusKSpipiPartRecoTagInfo.GetNumberPi0();
  m_TagNumberGamma = findKKpipiVersusKSpipiPartRecoTagInfo.GetNumberGamma();
  for(int j = 0; j < m_TagNumberGamma; j++) {
    m_TagExtraShowerEnergy[j] = findKKpipiVersusKSpipiPartRecoTagInfo.GetExtraShowerEnergy(j);
  }
  m_SignalKSFitSuccess = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSFitSuccess_KKpipi();
  m_SignalDecayLengthVeeVertex = findKKpipiVersusKSpipiPartRecoTagInfo.GetDecayLengthVeeVertex_KKpipi();
  m_SignalChi2VeeVertex = findKKpipiVersusKSpipiPartRecoTagInfo.GetChi2VeeVertex_KKpipi();
  m_SignalKSMassVeeVertex = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSMassVeeVertex_KKpipi();
  m_SignalDecayLengthFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetDecayLengthFit_KKpipi();
  m_SignalDecayLengthErrorFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetDecayLengthErrorFit_KKpipi();
  m_SignalChi2Fit = findKKpipiVersusKSpipiPartRecoTagInfo.GetChi2Fit_KKpipi();
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKKpipiVersusKSpipiPartRecoTagInfo.GetDaughterTrackID_KKpipi(), 3, this);
    m_SignalIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[3] = {m_SignalRecKCharge*321, 211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 3);
    m_SignalPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_SignalKTrueID = ReconstructedPID[0];
    m_SignalPiPlusTrueID = ReconstructedPID[1];
    m_SignalPiMinusTrueID = ReconstructedPID[2];
  }
  m_TagKSFitSuccess = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSFitSuccess_KSpipi();
  m_TagDecayLengthVeeVertex = findKKpipiVersusKSpipiPartRecoTagInfo.GetDecayLengthVeeVertex_KSpipi();
  m_TagChi2VeeVertex = findKKpipiVersusKSpipiPartRecoTagInfo.GetChi2VeeVertex_KSpipi();
  m_TagKSMassVeeVertex = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSMassVeeVertex_KSpipi();
  m_TagDecayLengthFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetDecayLengthFit_KSpipi();
  m_TagDecayLengthErrorFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetDecayLengthErrorFit_KSpipi();
  m_TagChi2Fit = findKKpipiVersusKSpipiPartRecoTagInfo.GetChi2Fit_KSpipi();
  m_TagKSPiPluspx = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSPiPlusP(0);
  m_TagKSPiPluspy = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSPiPlusP(1);
  m_TagKSPiPluspz = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSPiPlusP(2);
  m_TagKSPiPlusenergy = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSPiPlusP(3);
  m_TagKSPiMinuspx = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSPiMinusP(0);
  m_TagKSPiMinuspy = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSPiMinusP(1);
  m_TagKSPiMinuspz = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSPiMinusP(2);
  m_TagKSPiMinusenergy = findKKpipiVersusKSpipiPartRecoTagInfo.GetKSPiMinusP(3);
  m_TagKSpx = findKKpipiVersusKSpipiPartRecoTagInfo.GetKShortP(0);
  m_TagKSpy = findKKpipiVersusKSpipiPartRecoTagInfo.GetKShortP(1);
  m_TagKSpz = findKKpipiVersusKSpipiPartRecoTagInfo.GetKShortP(2);
  m_TagKSenergy = findKKpipiVersusKSpipiPartRecoTagInfo.GetKShortP(3);
  m_TagPiPluspx = findKKpipiVersusKSpipiPartRecoTagInfo.GethPlusP(0);
  m_TagPiPluspy = findKKpipiVersusKSpipiPartRecoTagInfo.GethPlusP(1);
  m_TagPiPluspz = findKKpipiVersusKSpipiPartRecoTagInfo.GethPlusP(2);
  m_TagPiPlusenergy = findKKpipiVersusKSpipiPartRecoTagInfo.GethPlusP(3);
  m_TagPiMinuspx = findKKpipiVersusKSpipiPartRecoTagInfo.GethMinusP(0);
  m_TagPiMinuspy = findKKpipiVersusKSpipiPartRecoTagInfo.GethMinusP(1);
  m_TagPiMinuspz = findKKpipiVersusKSpipiPartRecoTagInfo.GethMinusP(2);
  m_TagPiMinusenergy = findKKpipiVersusKSpipiPartRecoTagInfo.GethMinusP(3);
  m_TagKSpxKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKShortPKalmanFit(0);
  m_TagKSpyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKShortPKalmanFit(1);
  m_TagKSpzKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKShortPKalmanFit(2);
  m_TagKSenergyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GetKShortPKalmanFit(3);
  m_TagPiPluspxKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GethPlusPKalmanFit(0);
  m_TagPiPluspyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GethPlusPKalmanFit(1);
  m_TagPiPluspzKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GethPlusPKalmanFit(2);
  m_TagPiPlusenergyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GethPlusPKalmanFit(3);
  m_TagPiMinuspxKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GethMinusPKalmanFit(0);
  m_TagPiMinuspyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GethMinusPKalmanFit(1);
  m_TagPiMinuspzKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GethMinusPKalmanFit(2);
  m_TagPiMinusenergyKalmanFit = findKKpipiVersusKSpipiPartRecoTagInfo.GethMinusPKalmanFit(3);
  if(m_RunNumber < 0) {
    std::vector<int> DaughterTrackIDs = findKKpipiVersusKSpipiPartRecoTagInfo.GetDaughterTrackID_KSpipi();
    PIDTruth PID_Truth(DaughterTrackIDs, 4, this);
    m_TagIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {211, -211, 211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_TagPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_TagKSPiPlusTrueID = ReconstructedPID[0];
    m_TagKSPiMinusTrueID = ReconstructedPID[1];
    m_TagPiPlusTrueID = ReconstructedPID[2];
    m_TagPiMinusTrueID = ReconstructedPID[3];
    m_TagKSPiPlusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[0], true);
    m_TagKSPiMinusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[1], true);
  }
  return StatusCode::SUCCESS;
}
