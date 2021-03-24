// Martin Duy Tat 12th March 2021

// KKpipi
#include "KKpipi/KSpi0pi0SingleTag.h"
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

KSpi0pi0SingleTag::KSpi0pi0SingleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KSpi0pi0SingleTag::~KSpi0pi0SingleTag() {
}

StatusCode KSpi0pi0SingleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KSpi0pi0 Single Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KSpi0pi0SingleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KSpi0pi0SingleTag", CLID_ColumnWiseTuple, "Single tagged D->KSpi0pi0 events");
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
      status = m_tuple->addItem("HighEPi0px1", m_HighEPi0px1);
      status = m_tuple->addItem("HighEPi0py1", m_HighEPi0py1);
      status = m_tuple->addItem("HighEPi0pz1", m_HighEPi0pz1);
      status = m_tuple->addItem("HighEPi0energy1", m_HighEPi0energy1);
      status = m_tuple->addItem("LowEPi0px1", m_LowEPi0px1);
      status = m_tuple->addItem("LowEPi0py1", m_LowEPi0py1);
      status = m_tuple->addItem("LowEPi0pz1", m_LowEPi0pz1);
      status = m_tuple->addItem("LowEPi0energy1", m_LowEPi0energy1);
      status = m_tuple->addItem("HighEPi0Constrainedpx1", m_HighEPi0Constrainedpx1);
      status = m_tuple->addItem("HighEPi0Constrainedpy1", m_HighEPi0Constrainedpy1);
      status = m_tuple->addItem("HighEPi0Constrainedpz1", m_HighEPi0Constrainedpz1);
      status = m_tuple->addItem("HighEPi0Constrainedenergy1", m_HighEPi0Constrainedenergy1);
      status = m_tuple->addItem("LowEPi0Constrainedpx1", m_LowEPi0Constrainedpx1);
      status = m_tuple->addItem("LowEPi0Constrainedpy1", m_LowEPi0Constrainedpy1);
      status = m_tuple->addItem("LowEPi0Constrainedpz1", m_LowEPi0Constrainedpz1);
      status = m_tuple->addItem("LowEPi0Constrainedenergy1", m_LowEPi0Constrainedenergy1);
      status = m_tuple->addItem("Pi0Chi2Fit1", m_Pi0Chi2Fit1);
      status = m_tuple->addItem("HighEPi0px2", m_HighEPi0px2);
      status = m_tuple->addItem("HighEPi0py2", m_HighEPi0py2);
      status = m_tuple->addItem("HighEPi0pz2", m_HighEPi0pz2);
      status = m_tuple->addItem("HighEPi0energy2", m_HighEPi0energy2);
      status = m_tuple->addItem("LowEPi0px2", m_LowEPi0px2);
      status = m_tuple->addItem("LowEPi0py2", m_LowEPi0py2);
      status = m_tuple->addItem("LowEPi0pz2", m_LowEPi0pz2);
      status = m_tuple->addItem("LowEPi0energy2", m_LowEPi0energy2);
      status = m_tuple->addItem("HighEPi0Constrainedpx2", m_HighEPi0Constrainedpx2);
      status = m_tuple->addItem("HighEPi0Constrainedpy2", m_HighEPi0Constrainedpy2);
      status = m_tuple->addItem("HighEPi0Constrainedpz2", m_HighEPi0Constrainedpz2);
      status = m_tuple->addItem("HighEPi0Constrainedenergy2", m_HighEPi0Constrainedenergy2);
      status = m_tuple->addItem("LowEPi0Constrainedpx2", m_LowEPi0Constrainedpx2);
      status = m_tuple->addItem("LowEPi0Constrainedpy2", m_LowEPi0Constrainedpy2);
      status = m_tuple->addItem("LowEPi0Constrainedpz2", m_LowEPi0Constrainedpz2);
      status = m_tuple->addItem("LowEPi0Constrainedenergy2", m_LowEPi0Constrainedenergy2);
      status = m_tuple->addItem("Pi0Chi2Fit2", m_Pi0Chi2Fit2);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KSpi0pi0 Single Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KSpi0pi0SingleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KSpi0pi0 Single Tag Algorithm" << endreq;
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
  if(DTTool.findSTag(EvtRecDTag::kD0toKsPi0Pi0)) {
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

StatusCode KSpi0pi0SingleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KSpi0pi0 Single Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KSpi0pi0SingleTag::FillTuple(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
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
  FindPi0Eta findPi0(2);
  findPi0.findPi0Eta(DTTool_iter, DTTool);
  m_HighEPi0px1 = findPi0.GetHighEPhotonP(0, 0);
  m_HighEPi0py1 = findPi0.GetHighEPhotonP(1, 0);
  m_HighEPi0pz1 = findPi0.GetHighEPhotonP(2, 0);
  m_HighEPi0energy1 = findPi0.GetHighEPhotonP(3, 0);
  m_LowEPi0px1 = findPi0.GetLowEPhotonP(0, 0);
  m_LowEPi0py1 = findPi0.GetLowEPhotonP(1, 0);
  m_LowEPi0pz1 = findPi0.GetLowEPhotonP(2, 0);
  m_LowEPi0energy1 = findPi0.GetLowEPhotonP(3, 0);
  m_HighEPi0Constrainedpx1 = findPi0.GetHighEPhotonPConstrained(0, 0);
  m_HighEPi0Constrainedpy1 = findPi0.GetHighEPhotonPConstrained(1, 0);
  m_HighEPi0Constrainedpz1 = findPi0.GetHighEPhotonPConstrained(2, 0);
  m_HighEPi0Constrainedenergy1 = findPi0.GetHighEPhotonPConstrained(3, 0);
  m_LowEPi0Constrainedpx1 = findPi0.GetLowEPhotonPConstrained(0, 0);
  m_LowEPi0Constrainedpy1 = findPi0.GetLowEPhotonPConstrained(1, 0);
  m_LowEPi0Constrainedpz1 = findPi0.GetLowEPhotonPConstrained(2, 0);
  m_LowEPi0Constrainedenergy1 = findPi0.GetLowEPhotonPConstrained(3, 0);
  m_Pi0Chi2Fit1 = findPi0.GetChi2Fit(0);
  m_HighEPi0px2 = findPi0.GetHighEPhotonP(0, 1);
  m_HighEPi0py2 = findPi0.GetHighEPhotonP(1, 1);
  m_HighEPi0pz2 = findPi0.GetHighEPhotonP(2, 1);
  m_HighEPi0energy2 = findPi0.GetHighEPhotonP(3, 1);
  m_LowEPi0px2 = findPi0.GetLowEPhotonP(0, 1);
  m_LowEPi0py2 = findPi0.GetLowEPhotonP(1, 1);
  m_LowEPi0pz2 = findPi0.GetLowEPhotonP(2, 1);
  m_LowEPi0energy2 = findPi0.GetLowEPhotonP(3, 1);
  m_HighEPi0Constrainedpx2 = findPi0.GetHighEPhotonPConstrained(0, 1);
  m_HighEPi0Constrainedpy2 = findPi0.GetHighEPhotonPConstrained(1, 1);
  m_HighEPi0Constrainedpz2 = findPi0.GetHighEPhotonPConstrained(2, 1);
  m_HighEPi0Constrainedenergy2 = findPi0.GetHighEPhotonPConstrained(3, 1);
  m_LowEPi0Constrainedpx2 = findPi0.GetLowEPhotonPConstrained(0, 1);
  m_LowEPi0Constrainedpy2 = findPi0.GetLowEPhotonPConstrained(1, 1);
  m_LowEPi0Constrainedpz2 = findPi0.GetLowEPhotonPConstrained(2, 1);
  m_LowEPi0Constrainedenergy2 = findPi0.GetLowEPhotonPConstrained(3, 1);
  m_Pi0Chi2Fit2 = findPi0.GetChi2Fit(1);
  return StatusCode::SUCCESS;
}
