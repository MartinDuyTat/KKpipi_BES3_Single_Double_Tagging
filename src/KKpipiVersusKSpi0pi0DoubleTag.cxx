// Martin Duy Tat 3rd March 2021

// KKpipi
#include "KKpipi/KKpipiVersusKSpi0pi0DoubleTag.h"
#include "KKpipi/FindKKpipiTagInfo.h"
#include "KKpipi/FindKS.h"
#include "KKpipi/FindPi0Eta.h"
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
#include<utility>

KKpipiVersusKSpi0pi0DoubleTag::KKpipiVersusKSpi0pi0DoubleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKpipiVersusKSpi0pi0DoubleTag::~KKpipiVersusKSpi0pi0DoubleTag() {
}

StatusCode KKpipiVersusKSpi0pi0DoubleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KKpipi vs KSpi0pi0 Double Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KSpi0pi0DoubleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KSpi0pi0DoubleTag", CLID_ColumnWiseTuple, "Double tagged D->KKpipi vs D->KSpi0pi0 events");
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
      status = m_tuple->addItem("TagKSPiPluspx", m_TagKSPiPluspx);
      status = m_tuple->addItem("TagKSPiPluspy", m_TagKSPiPluspy);
      status = m_tuple->addItem("TagKSPiPluspz", m_TagKSPiPluspz);
      status = m_tuple->addItem("TagKSPiPlusenergy", m_TagKSPiPlusenergy);
      status = m_tuple->addItem("TagKSPiMinuspx", m_TagKSPiMinuspx);
      status = m_tuple->addItem("TagKSPiMinuspy", m_TagKSPiMinuspy);
      status = m_tuple->addItem("TagKSPiMinuspz", m_TagKSPiMinuspz);
      status = m_tuple->addItem("TagKSPiMinusenergy", m_TagKSPiMinusenergy);
      status = m_tuple->addItem("TagHighEPi0px1", m_TagHighEPi0px1);
      status = m_tuple->addItem("TagHighEPi0py1", m_TagHighEPi0py1);
      status = m_tuple->addItem("TagHighEPi0pz1", m_TagHighEPi0pz1);
      status = m_tuple->addItem("TagHighEPi0energy1", m_TagHighEPi0energy1);
      status = m_tuple->addItem("TagLowEPi0px1", m_TagLowEPi0px1);
      status = m_tuple->addItem("TagLowEPi0py1", m_TagLowEPi0py1);
      status = m_tuple->addItem("TagLowEPi0pz1", m_TagLowEPi0pz1);
      status = m_tuple->addItem("TagLowEPi0energy1", m_TagLowEPi0energy1);
      status = m_tuple->addItem("TagMgammagamma1", m_TagMgammagamma1);
      status = m_tuple->addItem("TagHighEPi0Constrainedpx1", m_TagHighEPi0Constrainedpx1);
      status = m_tuple->addItem("TagHighEPi0Constrainedpy1", m_TagHighEPi0Constrainedpy1);
      status = m_tuple->addItem("TagHighEPi0Constrainedpz1", m_TagHighEPi0Constrainedpz1);
      status = m_tuple->addItem("TagHighEPi0Constrainedenergy1", m_TagHighEPi0Constrainedenergy1);
      status = m_tuple->addItem("TagLowEPi0Constrainedpx1", m_TagLowEPi0Constrainedpx1);
      status = m_tuple->addItem("TagLowEPi0Constrainedpy1", m_TagLowEPi0Constrainedpy1);
      status = m_tuple->addItem("TagLowEPi0Constrainedpz1", m_TagLowEPi0Constrainedpz1);
      status = m_tuple->addItem("TagLowEPi0Constrainedenergy1", m_TagLowEPi0Constrainedenergy1);
      status = m_tuple->addItem("TagPi0Chi2Fit1", m_Pi0Chi2Fit1);
      status = m_tuple->addItem("TagHighEPi0px2", m_TagHighEPi0px2);
      status = m_tuple->addItem("TagHighEPi0py2", m_TagHighEPi0py2);
      status = m_tuple->addItem("TagHighEPi0pz2", m_TagHighEPi0pz2);
      status = m_tuple->addItem("TagHighEPi0energy2", m_TagHighEPi0energy2);
      status = m_tuple->addItem("TagLowEPi0px2", m_TagLowEPi0px2);
      status = m_tuple->addItem("TagLowEPi0py2", m_TagLowEPi0py2);
      status = m_tuple->addItem("TagLowEPi0pz2", m_TagLowEPi0pz2);
      status = m_tuple->addItem("TagLowEPi0energy2", m_TagLowEPi0energy2);
      status = m_tuple->addItem("TagMgammagamma2", m_TagMgammagamma2);
      status = m_tuple->addItem("TagHighEPi0Constrainedpx2", m_TagHighEPi0Constrainedpx2);
      status = m_tuple->addItem("TagHighEPi0Constrainedpy2", m_TagHighEPi0Constrainedpy2);
      status = m_tuple->addItem("TagHighEPi0Constrainedpz2", m_TagHighEPi0Constrainedpz2);
      status = m_tuple->addItem("TagHighEPi0Constrainedenergy2", m_TagHighEPi0Constrainedenergy2);
      status = m_tuple->addItem("TagLowEPi0Constrainedpx2", m_TagLowEPi0Constrainedpx2);
      status = m_tuple->addItem("TagLowEPi0Constrainedpy2", m_TagLowEPi0Constrainedpy2);
      status = m_tuple->addItem("TagLowEPi0Constrainedpz2", m_TagLowEPi0Constrainedpz2);
      status = m_tuple->addItem("TagLowEPi0Constrainedenergy2", m_TagLowEPi0Constrainedenergy2);
      status = m_tuple->addItem("TagPi0Chi2Fit2", m_Pi0Chi2Fit2);
      status = m_tuple->addItem("TagIsSameDMother", m_TagIsSameDMother);
      status = m_tuple->addItem("TagIsSameDMotherAll", m_TagIsSameDMotherAll);
      status = m_tuple->addItem("TagPIDTrue", m_TagPIDTrue);
      status = m_tuple->addItem("TagKSPiPlusTrueID", m_TagKSPiPlusTrueID);
      status = m_tuple->addItem("TagKSPiMinusTrueID", m_TagKSPiMinusTrueID);
      status = m_tuple->addItem("TagHighEPi0PhotonTrueID1", m_TagHighEPi0PhotonTrueID1);
      status = m_tuple->addItem("TagLowEPi0PhotonTrueID1", m_TagLowEPi0PhotonTrueID1);
      status = m_tuple->addItem("TagHighEPi0PhotonTrueID2", m_TagHighEPi0PhotonTrueID2);
      status = m_tuple->addItem("TagLowEPi0PhotonTrueID2", m_TagLowEPi0PhotonTrueID2);
      status = m_tuple->addItem("TagKSPiPlusMotherTrueID", m_TagKSPiPlusMotherTrueID);
      status = m_tuple->addItem("TagKSPiMinusMotherTrueID", m_TagKSPiMinusMotherTrueID);
      status = m_tuple->addItem("TagHighEPi0PhotonMotherTrueID1", m_TagHighEPi0PhotonMotherTrueID1);
      status = m_tuple->addItem("TagLowEPi0PhotonMotherTrueID1", m_TagLowEPi0PhotonMotherTrueID1);
      status = m_tuple->addItem("TagHighEPi0PhotonMotherTrueID2", m_TagHighEPi0PhotonMotherTrueID2);
      status = m_tuple->addItem("TagLowEPi0PhotonMotherTrueID2", m_TagLowEPi0PhotonMotherTrueID2);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KKpipi vs KSpi0pi0 Double Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KKpipiVersusKSpi0pi0DoubleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KKpipi vs KSpi0pi0 Double Tag Algorithm" << endreq;
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
  if(DTTool.findDTag(EvtRecDTag::kD0toKKPiPi, EvtRecDTag::kD0toKsPi0Pi0)) {
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
    m_RunNumber = eventHeader->runNumber();
    m_EventNumber = eventHeader->eventNumber();
    DTagToolIterator DTTool_Signal_iter = DTTool.dtag1();
    DTagToolIterator DTTool_Tag_iter = DTTool.dtag2();
    StatusCode FillTupleStatus = FillTuple(DTTool_Signal_iter, DTTool_Tag_iter, DTTool);
    if(FillTupleStatus != StatusCode::SUCCESS) {
      log << MSG::FATAL << "Assigning KSpi0pi0 tuple info failed" << endreq;
      return StatusCode::FAILURE;
    }
    m_tuple->write();
  }
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKSpi0pi0DoubleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KKpipi vs KSpi0pi0 Double Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKSpi0pi0DoubleTag::FillTuple(DTagToolIterator DTTool_Signal_iter, DTagToolIterator DTTool_Tag_iter, DTagTool &DTTool) {
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
  FindKS findKS(true);
  status = findKS.findKS(DTTool_Tag_iter, DTTool);
  if(status == StatusCode::SUCCESS) {
    m_TagKSFitSuccess = 1;
    m_TagDecayLengthFit = findKS.GetDecayLengthFit();
    m_TagDecayLengthErrorFit = findKS.GetDecayLengthErrorFit();
    m_TagChi2Fit = findKS.GetChi2Fit();
  } else {
    m_TagKSFitSuccess = 0;
  }
  m_TagDecayLengthVeeVertex = findKS.GetDecayLengthVeeVertex();
  m_TagChi2VeeVertex = findKS.GetChi2VeeVertex();
  m_TagKSMassVeeVertex = findKS.GetKSMassVeeVertex();
  m_TagKSPiPluspx = findKS.GetKSPiPlusP(0);
  m_TagKSPiPluspy = findKS.GetKSPiPlusP(1);
  m_TagKSPiPluspz = findKS.GetKSPiPlusP(2);
  m_TagKSPiPlusenergy = findKS.GetKSPiPlusP(3);
  m_TagKSPiMinuspx = findKS.GetKSPiMinusP(0);
  m_TagKSPiMinuspy = findKS.GetKSPiMinusP(1);
  m_TagKSPiMinuspz = findKS.GetKSPiMinusP(2);
  m_TagKSPiMinusenergy = findKS.GetKSPiMinusP(3);
  FindPi0Eta findPi0(2);
  findPi0.findPi0Eta(DTTool_Tag_iter, DTTool);
  m_TagHighEPi0px1 = findPi0.GetHighEPhotonP(0, 0);
  m_TagHighEPi0py1 = findPi0.GetHighEPhotonP(1, 0);
  m_TagHighEPi0pz1 = findPi0.GetHighEPhotonP(2, 0);
  m_TagHighEPi0energy1 = findPi0.GetHighEPhotonP(3, 0);
  m_TagLowEPi0px1 = findPi0.GetLowEPhotonP(0, 0);
  m_TagLowEPi0py1 = findPi0.GetLowEPhotonP(1, 0);
  m_TagLowEPi0pz1 = findPi0.GetLowEPhotonP(2, 0);
  m_TagLowEPi0energy1 = findPi0.GetLowEPhotonP(3, 0);
  m_TagMgammagamma1 = findPi0.GetMgammagamma(0);
  m_TagHighEPi0Constrainedpx1 = findPi0.GetHighEPhotonPConstrained(0, 0);
  m_TagHighEPi0Constrainedpy1 = findPi0.GetHighEPhotonPConstrained(1, 0);
  m_TagHighEPi0Constrainedpz1 = findPi0.GetHighEPhotonPConstrained(2, 0);
  m_TagHighEPi0Constrainedenergy1 = findPi0.GetHighEPhotonPConstrained(3, 0);
  m_TagLowEPi0Constrainedpx1 = findPi0.GetLowEPhotonPConstrained(0, 0);
  m_TagLowEPi0Constrainedpy1 = findPi0.GetLowEPhotonPConstrained(1, 0);
  m_TagLowEPi0Constrainedpz1 = findPi0.GetLowEPhotonPConstrained(2, 0);
  m_TagLowEPi0Constrainedenergy1 = findPi0.GetLowEPhotonPConstrained(3, 0);
  m_Pi0Chi2Fit1 = findPi0.GetChi2Fit(0);
  m_TagHighEPi0px2 = findPi0.GetHighEPhotonP(0, 1);
  m_TagHighEPi0py2 = findPi0.GetHighEPhotonP(1, 1);
  m_TagHighEPi0pz2 = findPi0.GetHighEPhotonP(2, 1);
  m_TagHighEPi0energy2 = findPi0.GetHighEPhotonP(3, 1);
  m_TagLowEPi0px2 = findPi0.GetLowEPhotonP(0, 1);
  m_TagLowEPi0py2 = findPi0.GetLowEPhotonP(1, 1);
  m_TagLowEPi0pz2 = findPi0.GetLowEPhotonP(2, 1);
  m_TagLowEPi0energy2 = findPi0.GetLowEPhotonP(3, 1);
  m_TagMgammagamma2 = findPi0.GetMgammagamma(1);
  m_TagHighEPi0Constrainedpx2 = findPi0.GetHighEPhotonPConstrained(0, 1);
  m_TagHighEPi0Constrainedpy2 = findPi0.GetHighEPhotonPConstrained(1, 1);
  m_TagHighEPi0Constrainedpz2 = findPi0.GetHighEPhotonPConstrained(2, 1);
  m_TagHighEPi0Constrainedenergy2 = findPi0.GetHighEPhotonPConstrained(3, 1);
  m_TagLowEPi0Constrainedpx2 = findPi0.GetLowEPhotonPConstrained(0, 1);
  m_TagLowEPi0Constrainedpy2 = findPi0.GetLowEPhotonPConstrained(1, 1);
  m_TagLowEPi0Constrainedpz2 = findPi0.GetLowEPhotonPConstrained(2, 1);
  m_TagLowEPi0Constrainedenergy2 = findPi0.GetLowEPhotonPConstrained(3, 1);
  m_Pi0Chi2Fit2 = findPi0.GetChi2Fit(1);
  if(m_RunNumber < 0) {
    std::vector<int> DaughterTrackIDs = findKS.GetDaughterTrackIDs();
    PIDTruth PID_Truth(DaughterTrackIDs, 2, this);
    m_TagIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    std::vector<std::pair<int, int> > PhotonPairTrackID;
    PhotonPairTrackID.push_back(std::make_pair(findPi0.GetHighEPhotonTrackID(0), findPi0.GetLowEPhotonTrackID(0)));
    PhotonPairTrackID.push_back(std::make_pair(findPi0.GetHighEPhotonTrackID(1), findPi0.GetLowEPhotonTrackID(1)));
    PID_Truth = PIDTruth(DaughterTrackIDs, 2, this, PhotonPairTrackID);
    m_TagIsSameDMotherAll = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[6] = {211, -211, 0, 0, 0, 0};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 6);
    m_TagPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_TagKSPiPlusTrueID = ReconstructedPID[0];
    m_TagKSPiMinusTrueID = ReconstructedPID[1];
    m_TagHighEPi0PhotonTrueID1 = ReconstructedPID[2];
    m_TagLowEPi0PhotonTrueID1 = ReconstructedPID[3];
    m_TagHighEPi0PhotonTrueID2 = ReconstructedPID[4];
    m_TagLowEPi0PhotonTrueID2 = ReconstructedPID[5];
    m_TagKSPiPlusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[0], true);
    m_TagKSPiMinusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[1], true);
    m_TagHighEPi0PhotonMotherTrueID1 = PID_Truth.GetTrueMotherID(DaughterTrackIDs[2], false);
    m_TagLowEPi0PhotonMotherTrueID1 = PID_Truth.GetTrueMotherID(DaughterTrackIDs[3], false);
    m_TagHighEPi0PhotonMotherTrueID2 = PID_Truth.GetTrueMotherID(DaughterTrackIDs[4], false);
    m_TagLowEPi0PhotonMotherTrueID2 = PID_Truth.GetTrueMotherID(DaughterTrackIDs[5], false);
  }
  return StatusCode::SUCCESS;
}
