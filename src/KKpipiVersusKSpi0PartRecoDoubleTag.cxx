// Martin Duy Tat 18th January 2023

// KKpipi
#include "KKpipi/KKpipiVersusKSpi0PartRecoDoubleTag.h"
#include "KKpipi/FindKKpipiVersusKSpi0PartRecoTagInfo.h"
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

KKpipiVersusKSpi0PartRecoDoubleTag::KKpipiVersusKSpi0PartRecoDoubleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKpipiVersusKSpi0PartRecoDoubleTag::~KKpipiVersusKSpi0PartRecoDoubleTag() {
}

StatusCode KKpipiVersusKSpi0PartRecoDoubleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KKpipi vs KSpi0 PartReco Double Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KSpi0PartRecoDoubleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KSpi0PartRecoDoubleTag", CLID_ColumnWiseTuple, "Double tagged D->KKpipi vs D->KSpi0 events with missing K in KKpipi");
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
      status = m_tuple->addItem("TagHighEPi0px", m_TagHighEPi0px);
      status = m_tuple->addItem("TagHighEPi0py", m_TagHighEPi0py);
      status = m_tuple->addItem("TagHighEPi0pz", m_TagHighEPi0pz);
      status = m_tuple->addItem("TagHighEPi0energy", m_TagHighEPi0energy);
      status = m_tuple->addItem("TagLowEPi0px", m_TagLowEPi0px);
      status = m_tuple->addItem("TagLowEPi0py", m_TagLowEPi0py);
      status = m_tuple->addItem("TagLowEPi0pz", m_TagLowEPi0pz);
      status = m_tuple->addItem("TagLowEPi0energy", m_TagLowEPi0energy);
      status = m_tuple->addItem("TagMgammagamma", m_TagMgammagamma);
      status = m_tuple->addItem("TagHighEPi0Constrainedpx", m_TagHighEPi0Constrainedpx);
      status = m_tuple->addItem("TagHighEPi0Constrainedpy", m_TagHighEPi0Constrainedpy);
      status = m_tuple->addItem("TagHighEPi0Constrainedpz", m_TagHighEPi0Constrainedpz);
      status = m_tuple->addItem("TagHighEPi0Constrainedenergy", m_TagHighEPi0Constrainedenergy);
      status = m_tuple->addItem("TagLowEPi0Constrainedpx", m_TagLowEPi0Constrainedpx);
      status = m_tuple->addItem("TagLowEPi0Constrainedpy", m_TagLowEPi0Constrainedpy);
      status = m_tuple->addItem("TagLowEPi0Constrainedpz", m_TagLowEPi0Constrainedpz);
      status = m_tuple->addItem("TagLowEPi0Constrainedenergy", m_TagLowEPi0Constrainedenergy);
      status = m_tuple->addItem("TagPi0Chi2Fit", m_Pi0Chi2Fit);
      status = m_tuple->addItem("TagIsSameDMother", m_TagIsSameDMother);
      status = m_tuple->addItem("TagIsSameDMotherAll", m_TagIsSameDMotherAll);
      status = m_tuple->addItem("TagPIDTrue", m_TagPIDTrue);
      status = m_tuple->addItem("TagKSPiPlusTrueID", m_TagKSPiPlusTrueID);
      status = m_tuple->addItem("TagKSPiMinusTrueID", m_TagKSPiMinusTrueID);
      status = m_tuple->addItem("TagHighEPi0PhotonTrueID", m_TagHighEPi0PhotonTrueID);
      status = m_tuple->addItem("TagLowEPi0PhotonTrueID", m_TagLowEPi0PhotonTrueID);
      status = m_tuple->addItem("TagKSPiPlusMotherTrueID", m_TagKSPiPlusMotherTrueID);
      status = m_tuple->addItem("TagKSPiMinusMotherTrueID", m_TagKSPiMinusMotherTrueID);
      status = m_tuple->addItem("TagHighEPi0PhotonMotherTrueID", m_TagHighEPi0PhotonMotherTrueID);
      status = m_tuple->addItem("TagLowEPi0PhotonMotherTrueID", m_TagLowEPi0PhotonMotherTrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KKpipi vs KSpi0 Part Reco Double Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KKpipiVersusKSpi0PartRecoDoubleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KKpipi vs KSpi0 PartReco Double Tag Algorithm" << endreq;
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
  if(DTTool.findSTag(EvtRecDTag::kD0toKsPi0)) {
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

StatusCode KKpipiVersusKSpi0PartRecoDoubleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KKpipi vs KSpi0 Double Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKSpi0PartRecoDoubleTag::FillTuple(DTagToolIterator DTTool_Tag_iter, DTagTool &DTTool) {
  FindKKpipiVersusKSpi0PartRecoTagInfo findKKpipiVersusKSpi0PartRecoTagInfo;
  StatusCode status = findKKpipiVersusKSpi0PartRecoTagInfo.CalculateTagInfo(DTTool_Tag_iter, DTTool);
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
  m_SignalKPluspx = findKKpipiVersusKSpi0PartRecoTagInfo.GetKPlusP(0);
  m_SignalKPluspy = findKKpipiVersusKSpi0PartRecoTagInfo.GetKPlusP(1);
  m_SignalKPluspz = findKKpipiVersusKSpi0PartRecoTagInfo.GetKPlusP(2);
  m_SignalKPlusenergy = findKKpipiVersusKSpi0PartRecoTagInfo.GetKPlusP(3);
  m_SignalKMinuspx = findKKpipiVersusKSpi0PartRecoTagInfo.GetKMinusP(0);
  m_SignalKMinuspy = findKKpipiVersusKSpi0PartRecoTagInfo.GetKMinusP(1);
  m_SignalKMinuspz = findKKpipiVersusKSpi0PartRecoTagInfo.GetKMinusP(2);
  m_SignalKMinusenergy = findKKpipiVersusKSpi0PartRecoTagInfo.GetKMinusP(3);
  m_SignalPiPluspx = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiPlusP(0);
  m_SignalPiPluspy = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiPlusP(1);
  m_SignalPiPluspz = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiPlusP(2);
  m_SignalPiPlusenergy = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiPlusP(3);
  m_SignalPiMinuspx = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiMinusP(0);
  m_SignalPiMinuspy = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiMinusP(1);
  m_SignalPiMinuspz = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiMinusP(2);
  m_SignalPiMinusenergy = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiMinusP(3);
  m_SignalKalmanFitSuccess = findKKpipiVersusKSpi0PartRecoTagInfo.GetKalmanFitSuccess();
  m_SignalKalmanFitChi2 = findKKpipiVersusKSpi0PartRecoTagInfo.GetKalmanFitChi2();
  m_SignalKPluspxKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetKPlusPKalmanFit(0);
  m_SignalKPluspyKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetKPlusPKalmanFit(1);
  m_SignalKPluspzKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetKPlusPKalmanFit(2);
  m_SignalKPlusenergyKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetKPlusPKalmanFit(3);
  m_SignalKMinuspxKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetKMinusPKalmanFit(0);
  m_SignalKMinuspyKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetKMinusPKalmanFit(1);
  m_SignalKMinuspzKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetKMinusPKalmanFit(2);
  m_SignalKMinusenergyKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetKMinusPKalmanFit(3);
  m_SignalPiPluspxKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiPlusPKalmanFit(0);
  m_SignalPiPluspyKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiPlusPKalmanFit(1);
  m_SignalPiPluspzKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiPlusPKalmanFit(2);
  m_SignalPiPlusenergyKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiPlusPKalmanFit(3);
  m_SignalPiMinuspxKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiMinusPKalmanFit(0);
  m_SignalPiMinuspyKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiMinusPKalmanFit(1);
  m_SignalPiMinuspzKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiMinusPKalmanFit(2);
  m_SignalPiMinusenergyKalmanFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetPiMinusPKalmanFit(3);
  m_SignalMpipi = findKKpipiVersusKSpi0PartRecoTagInfo.GetMpipi_KKpipi();
  m_SignalMissingMass2 = findKKpipiVersusKSpi0PartRecoTagInfo.GetMissingMass2();
  m_SignalRecKCharge = findKKpipiVersusKSpi0PartRecoTagInfo.GetRecKCharge();
  m_SignalNumberPi0 = findKKpipiVersusKSpi0PartRecoTagInfo.GetNumberPi0();
  m_SignalKSFitSuccess = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSFitSuccess_KKpipi();
  m_SignalDecayLengthVeeVertex = findKKpipiVersusKSpi0PartRecoTagInfo.GetDecayLengthVeeVertex_KKpipi();
  m_SignalChi2VeeVertex = findKKpipiVersusKSpi0PartRecoTagInfo.GetChi2VeeVertex_KKpipi();
  m_SignalKSMassVeeVertex = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSMassVeeVertex_KKpipi();
  m_SignalDecayLengthFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetDecayLengthFit_KKpipi();
  m_SignalDecayLengthErrorFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetDecayLengthErrorFit_KKpipi();
  m_SignalChi2Fit = findKKpipiVersusKSpi0PartRecoTagInfo.GetChi2Fit_KKpipi();
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKKpipiVersusKSpi0PartRecoTagInfo.GetDaughterTrackID_KKpipi(), 3, this);
    m_SignalIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[3] = {m_SignalRecKCharge*321, 211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 3);
    m_SignalPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_SignalKTrueID = ReconstructedPID[0];
    m_SignalPiPlusTrueID = ReconstructedPID[1];
    m_SignalPiMinusTrueID = ReconstructedPID[2];
  }
  m_TagKSFitSuccess = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSFitSuccess_KSpi0();
  m_TagDecayLengthVeeVertex = findKKpipiVersusKSpi0PartRecoTagInfo.GetDecayLengthVeeVertex_KSpi0();
  m_TagChi2VeeVertex = findKKpipiVersusKSpi0PartRecoTagInfo.GetChi2VeeVertex_KSpi0();
  m_TagKSMassVeeVertex = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSMassVeeVertex_KSpi0();
  m_TagDecayLengthFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetDecayLengthFit_KSpi0();
  m_TagDecayLengthErrorFit = findKKpipiVersusKSpi0PartRecoTagInfo.GetDecayLengthErrorFit_KSpi0();
  m_TagChi2Fit = findKKpipiVersusKSpi0PartRecoTagInfo.GetChi2Fit_KSpi0();
  m_TagKSPiPluspx = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSPiPlusP(0);
  m_TagKSPiPluspy = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSPiPlusP(1);
  m_TagKSPiPluspz = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSPiPlusP(2);
  m_TagKSPiPlusenergy = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSPiPlusP(3);
  m_TagKSPiMinuspx = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSPiMinusP(0);
  m_TagKSPiMinuspy = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSPiMinusP(1);
  m_TagKSPiMinuspz = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSPiMinusP(2);
  m_TagKSPiMinusenergy = findKKpipiVersusKSpi0PartRecoTagInfo.GetKSPiMinusP(3);
  m_TagKSpx = findKKpipiVersusKSpi0PartRecoTagInfo.GetKShortP(0);
  m_TagKSpy = findKKpipiVersusKSpi0PartRecoTagInfo.GetKShortP(1);
  m_TagKSpz = findKKpipiVersusKSpi0PartRecoTagInfo.GetKShortP(2);
  m_TagKSenergy = findKKpipiVersusKSpi0PartRecoTagInfo.GetKShortP(3);
  m_TagHighEPi0px = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetHighEPhotonP(0);
  m_TagHighEPi0py = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetHighEPhotonP(1);
  m_TagHighEPi0pz = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetHighEPhotonP(2);
  m_TagHighEPi0energy = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetHighEPhotonP(3);
  m_TagLowEPi0px = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetLowEPhotonP(0);
  m_TagLowEPi0py = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetLowEPhotonP(1);
  m_TagLowEPi0pz = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetLowEPhotonP(2);
  m_TagLowEPi0energy = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetLowEPhotonP(3);
  m_TagMgammagamma = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetMgammagamma();
  m_TagHighEPi0Constrainedpx = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetHighEPhotonPConstrained(0);
  m_TagHighEPi0Constrainedpy = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetHighEPhotonPConstrained(1);
  m_TagHighEPi0Constrainedpz = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetHighEPhotonPConstrained(2);
  m_TagHighEPi0Constrainedenergy = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetHighEPhotonPConstrained(3);
  m_TagLowEPi0Constrainedpx = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetLowEPhotonPConstrained(0);
  m_TagLowEPi0Constrainedpy = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetLowEPhotonPConstrained(1);
  m_TagLowEPi0Constrainedpz = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetLowEPhotonPConstrained(2);
  m_TagLowEPi0Constrainedenergy = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetLowEPhotonPConstrained(3);
  m_Pi0Chi2Fit = findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info.GetChi2Fit();
  if(m_RunNumber < 0) {
    std::vector<int> DaughterTrackIDs = findKKpipiVersusKSpi0PartRecoTagInfo.GetDaughterTrackID_KSpi0();
    PIDTruth PID_Truth(DaughterTrackIDs, 2, this);
    m_TagIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    std::vector<std::pair<int, int> > PhotonPairTrackID;
    PhotonPairTrackID.push_back(std::make_pair(findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info().GetHighEPhotonTrackID(), findKKpipiVersusKSpi0PartRecoTagInfo.GetPi0Info().GetLowEPhotonTrackID()));
    PID_Truth = PIDTruth(DaughterTrackIDs, 2, this, PhotonPairTrackID);
    m_TagIsSameDMotherAll = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {211, -211, 0, 0};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_TagPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_TagKSPiPlusTrueID = ReconstructedPID[0];
    m_TagKSPiMinusTrueID = ReconstructedPID[1];
    m_TagHighEPi0PhotonTrueID = ReconstructedPID[2];
    m_TagLowEPi0PhotonTrueID = ReconstructedPID[3];
    m_TagKSPiPlusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[0], true);
    m_TagKSPiMinusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[1], true);
    m_TagHighEPi0PhotonMotherTrueID = PID_Truth.GetTrueMotherID(PhotonPairTrackID[0].first, false);
    m_TagLowEPi0PhotonMotherTrueID = PID_Truth.GetTrueMotherID(PhotonPairTrackID[0].second, false);
  }
  return StatusCode::SUCCESS;
}
