// Martin Duy Tat 12th March 2021

// KKpipi
#include "KKpipi/KSetaSingleTag.h"
#include "KKpipi/FindKS.h"
#include "KKpipi/FindPi0Eta.h"
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

KSetaSingleTag::KSetaSingleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KSetaSingleTag::~KSetaSingleTag() {
}

StatusCode KSetaSingleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KSeta Single Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KSetaSingleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KSetaSingleTag", CLID_ColumnWiseTuple, "Single tagged D->KSeta events");
    if(m_tuple) {
      status = m_tuple->addItem("Run", m_RunNumber);
      status = m_tuple->addItem("Event", m_EventNumber);
      status = m_tuple->addItem("NumberOfParticles", m_NumberParticles, 0, 100);
      status = m_tuple->addIndexedItem("ParticleIDs", m_NumberParticles, m_pdgID);
      status = m_tuple->addIndexedItem("MotherIndex", m_NumberParticles, m_MotherIndex);
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
      status = m_tuple->addItem("KSPiPluspx", m_KSPiPluspx);
      status = m_tuple->addItem("KSPiPluspy", m_KSPiPluspy);
      status = m_tuple->addItem("KSPiPluspz", m_KSPiPluspz);
      status = m_tuple->addItem("KSPiPlusenergy", m_KSPiPlusenergy);
      status = m_tuple->addItem("KSPiMinuspx", m_KSPiMinuspx);
      status = m_tuple->addItem("KSPiMinuspy", m_KSPiMinuspy);
      status = m_tuple->addItem("KSPiMinuspz", m_KSPiMinuspz);
      status = m_tuple->addItem("KSPiMinusenergy", m_KSPiMinusenergy);
      status = m_tuple->addItem("KSPiPluspxFit", m_KSPiPluspxFit);
      status = m_tuple->addItem("KSPiPluspyFit", m_KSPiPluspyFit);
      status = m_tuple->addItem("KSPiPluspzFit", m_KSPiPluspzFit);
      status = m_tuple->addItem("KSPiPlusenergyFit", m_KSPiPlusenergyFit);
      status = m_tuple->addItem("KSPiMinuspxFit", m_KSPiMinuspxFit);
      status = m_tuple->addItem("KSPiMinuspyFit", m_KSPiMinuspyFit);
      status = m_tuple->addItem("KSPiMinuspzFit", m_KSPiMinuspzFit);
      status = m_tuple->addItem("KSPiMinusenergyFit", m_KSPiMinusenergyFit);
      status = m_tuple->addItem("HighEEtapx", m_HighEEtapx);
      status = m_tuple->addItem("HighEEtapy", m_HighEEtapy);
      status = m_tuple->addItem("HighEEtapz", m_HighEEtapz);
      status = m_tuple->addItem("HighEEtaenergy", m_HighEEtaenergy);
      status = m_tuple->addItem("LowEEtapx", m_LowEEtapx);
      status = m_tuple->addItem("LowEEtapy", m_LowEEtapy);
      status = m_tuple->addItem("LowEEtapz", m_LowEEtapz);
      status = m_tuple->addItem("LowEEtaenergy", m_LowEEtaenergy);
      status = m_tuple->addItem("HighEEtaConstrainedpx", m_HighEEtaConstrainedpx);
      status = m_tuple->addItem("HighEEtaConstrainedpy", m_HighEEtaConstrainedpy);
      status = m_tuple->addItem("HighEEtaConstrainedpz", m_HighEEtaConstrainedpz);
      status = m_tuple->addItem("HighEEtaConstrainedenergy", m_HighEEtaConstrainedenergy);
      status = m_tuple->addItem("LowEEtaConstrainedpx", m_LowEEtaConstrainedpx);
      status = m_tuple->addItem("LowEEtaConstrainedpy", m_LowEEtaConstrainedpy);
      status = m_tuple->addItem("LowEEtaConstrainedpz", m_LowEEtaConstrainedpz);
      status = m_tuple->addItem("LowEEtaConstrainedenergy", m_LowEEtaConstrainedenergy);
      status = m_tuple->addItem("EtaChi2Fit", m_EtaChi2Fit);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KSeta Single Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KSetaSingleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KSeta Single Tag Algorithm" << endreq;
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
  if(DTTool.findSTag(EvtRecDTag::kD0toKsEta)) {
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

StatusCode KSetaSingleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KSeta Single Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KSetaSingleTag::FillTuple(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
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
  }
  m_DMass = (*DTTool_iter)->mass();
  m_MBC = (*DTTool_iter)->mBC();
  m_DeltaE = (*DTTool_iter)->deltaE();
  m_BeamE = (*DTTool_iter)->beamE();
  m_Dpx = (*DTTool_iter)->p4().x();
  m_Dpy = (*DTTool_iter)->p4().y();
  m_Dpz = (*DTTool_iter)->p4().z();
  m_Denergy = (*DTTool_iter)->p4().t();
  FindKS findKS(true);
  StatusCode status = findKS.findKS(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_DecayLengthVeeVertex = findKS.GetDecayLengthVeeVertex();
  m_Chi2VeeVertex = findKS.GetChi2VeeVertex();
  m_KSMassVeeVertex = findKS.GetKSMassVeeVertex();
  m_DecayLengthFit = findKS.GetDecayLengthFit();
  m_DecayLengthErrorFit = findKS.GetDecayLengthErrorFit();
  m_Chi2Fit = findKS.GetChi2Fit();
  m_KSMassFit = findKS.GetKSMassFit();
  m_KSPiPluspx = findKS.GetKSPiPlusP(0);
  m_KSPiPluspy = findKS.GetKSPiPlusP(1);
  m_KSPiPluspz = findKS.GetKSPiPlusP(2);
  m_KSPiPlusenergy = findKS.GetKSPiPlusP(3);
  m_KSPiMinuspx = findKS.GetKSPiMinusP(0);
  m_KSPiMinuspy = findKS.GetKSPiMinusP(1);
  m_KSPiMinuspz = findKS.GetKSPiMinusP(2);
  m_KSPiMinusenergy = findKS.GetKSPiMinusP(3);
  m_KSPiPluspxFit = findKS.GetKSPiPlusPFit(0);
  m_KSPiPluspyFit = findKS.GetKSPiPlusPFit(1);
  m_KSPiPluspzFit = findKS.GetKSPiPlusPFit(2);
  m_KSPiPlusenergyFit = findKS.GetKSPiPlusPFit(3);
  m_KSPiMinuspxFit = findKS.GetKSPiMinusPFit(0);
  m_KSPiMinuspyFit = findKS.GetKSPiMinusPFit(1);
  m_KSPiMinuspzFit = findKS.GetKSPiMinusPFit(2);
  m_KSPiMinusenergyFit = findKS.GetKSPiMinusPFit(3);
  FindPi0Eta findEta(1, "eta");
  findEta.findPi0Eta(DTTool_iter, DTTool);
  m_HighEEtapx = findEta.GetHighEPhotonP(0);
  m_HighEEtapy = findEta.GetHighEPhotonP(1);
  m_HighEEtapz = findEta.GetHighEPhotonP(2);
  m_HighEEtaenergy = findEta.GetHighEPhotonP(3);
  m_LowEEtapx = findEta.GetLowEPhotonP(0);
  m_LowEEtapy = findEta.GetLowEPhotonP(1);
  m_LowEEtapz = findEta.GetLowEPhotonP(2);
  m_LowEEtaenergy = findEta.GetLowEPhotonP(3);
  m_HighEEtaConstrainedpx = findEta.GetHighEPhotonPConstrained(0);
  m_HighEEtaConstrainedpy = findEta.GetHighEPhotonPConstrained(1);
  m_HighEEtaConstrainedpz = findEta.GetHighEPhotonPConstrained(2);
  m_HighEEtaConstrainedenergy = findEta.GetHighEPhotonPConstrained(3);
  m_LowEEtaConstrainedpx = findEta.GetLowEPhotonPConstrained(0);
  m_LowEEtaConstrainedpy = findEta.GetLowEPhotonPConstrained(1);
  m_LowEEtaConstrainedpz = findEta.GetLowEPhotonPConstrained(2);
  m_LowEEtaConstrainedenergy = findEta.GetLowEPhotonPConstrained(3);
  m_EtaChi2Fit = findEta.GetChi2Fit();
  return StatusCode::SUCCESS;
}
