// Martin Duy Tat 12th March 2021

// KKpipi
#include "KKpipi/pipipi0SingleTag.h"
#include "KKpipi/FindhhTagInfo.h"
#include "KKpipi/FindKS.h"
#include "KKpipi/FindMCInfo.h"
#include "KKpipi/FindPi0Eta.h"
#include "KKpipi/ParticleMasses.h"
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
// ROOT
#include "TMath.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McTruth/McParticle.h"
// STL
#include<vector>
#include<string>

pipipi0SingleTag::pipipi0SingleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

pipipi0SingleTag::~pipipi0SingleTag() {
}

StatusCode pipipi0SingleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing pipipi0 Single Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/pipipi0SingleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/pipipi0SingleTag", CLID_ColumnWiseTuple, "Single tagged D->pipipi0 events");
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
      status = m_tuple->addItem("PiPluspx", m_PiPluspx);
      status = m_tuple->addItem("PiPluspy", m_PiPluspy);
      status = m_tuple->addItem("PiPluspz", m_PiPluspz);
      status = m_tuple->addItem("PiPlusenergy", m_PiPlusenergy);
      status = m_tuple->addItem("PiMinuspx", m_PiMinuspx);
      status = m_tuple->addItem("PiMinuspy", m_PiMinuspy);
      status = m_tuple->addItem("PiMinuspz", m_PiMinuspz);
      status = m_tuple->addItem("PiMinusenergy", m_PiMinusenergy);
      status = m_tuple->addItem("HighEPi0px", m_HighEPi0px);
      status = m_tuple->addItem("HighEPi0py", m_HighEPi0py);
      status = m_tuple->addItem("HighEPi0pz", m_HighEPi0pz);
      status = m_tuple->addItem("HighEPi0energy", m_HighEPi0energy);
      status = m_tuple->addItem("LowEPi0px", m_LowEPi0px);
      status = m_tuple->addItem("LowEPi0py", m_LowEPi0py);
      status = m_tuple->addItem("LowEPi0pz", m_LowEPi0pz);
      status = m_tuple->addItem("LowEPi0energy", m_LowEPi0energy);
      status = m_tuple->addItem("HighEPi0Constrainedpx", m_HighEPi0Constrainedpx);
      status = m_tuple->addItem("HighEPi0Constrainedpy", m_HighEPi0Constrainedpy);
      status = m_tuple->addItem("HighEPi0Constrainedpz", m_HighEPi0Constrainedpz);
      status = m_tuple->addItem("HighEPi0Constrainedenergy", m_HighEPi0Constrainedenergy);
      status = m_tuple->addItem("LowEPi0Constrainedpx", m_LowEPi0Constrainedpx);
      status = m_tuple->addItem("LowEPi0Constrainedpy", m_LowEPi0Constrainedpy);
      status = m_tuple->addItem("LowEPi0Constrainedpz", m_LowEPi0Constrainedpz);
      status = m_tuple->addItem("LowEPi0Constrainedenergy", m_LowEPi0Constrainedenergy);
      status = m_tuple->addItem("Pi0Chi2Fit", m_Pi0Chi2Fit);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for pipipi0 Single Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode pipipi0SingleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing pipipi0 Single Tag Algorithm" << endreq;
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
  if(DTTool.findSTag(EvtRecDTag::kD0toPiPiPi0)) {
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

StatusCode pipipi0SingleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing pipipi0 Single Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode pipipi0SingleTag::FillTuple(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
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
  FindhhTagInfo findpipiTagInfo("pipi");
  StatusCode status = findpipiTagInfo.CalculateTagInfo(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_PiPluspx = findpipiTagInfo.GethPlusP(0);
  m_PiPluspy = findpipiTagInfo.GethPlusP(1);
  m_PiPluspz = findpipiTagInfo.GethPlusP(2);
  m_PiPlusenergy = findpipiTagInfo.GethPlusP(3);
  m_PiMinuspx = findpipiTagInfo.GethMinusP(0);
  m_PiMinuspy = findpipiTagInfo.GethMinusP(1);
  m_PiMinuspz = findpipiTagInfo.GethMinusP(2);
  m_PiMinusenergy = findpipiTagInfo.GethMinusP(3);
  SmartRefVector<EvtRecTrack> Tracks = (*DTTool_iter)->tracks();
  double Mpipi = TMath::Sqrt(TMath::Power(m_PiPlusenergy + m_PiMinusenergy, 2) - TMath::Power(m_PiPluspx + m_PiMinuspx, 2) - TMath::Power(m_PiPluspy + m_PiMinuspy, 2) - TMath::Power(m_PiPluspz + m_PiMinuspz, 2));
  m_KSFitSuccess = 0;
  if(TMath::Abs(Mpipi - MASS::KS_MASS) < 0.020) {
    FindKS findKS(false);
    std::vector<int> PionTrackIDs;
    for(SmartRefVector<EvtRecTrack>::iterator Track_iter = Tracks.begin(); Track_iter != Tracks.end(); Track_iter++) {
      if(DTTool.isPion(*Track_iter)) {
	PionTrackIDs.push_back((*Track_iter)->trackId());
      }
    }
    StatusCode statuscode = findKS.findKS(DTTool_iter, DTTool, PionTrackIDs);
    m_KSFitSuccess = 0;
    if(statuscode == StatusCode::SUCCESS) {
      m_KSFitSuccess = 1;
      m_DecayLengthVeeVertex = findKS.GetDecayLengthVeeVertex();
      m_Chi2VeeVertex = findKS.GetChi2VeeVertex();
      m_KSMassVeeVertex = findKS.GetKSMassVeeVertex();
      m_DecayLengthFit = findKS.GetDecayLengthFit();
      m_DecayLengthErrorFit = findKS.GetDecayLengthErrorFit();
      m_Chi2Fit = findKS.GetChi2Fit();
      m_KSMassFit = findKS.GetKSMassFit();
    }
  }
  FindPi0Eta findPi0;
  findPi0.findPi0Eta(DTTool_iter, DTTool);
  m_HighEPi0px = findPi0.GetHighEPhotonP(0);
  m_HighEPi0py = findPi0.GetHighEPhotonP(1);
  m_HighEPi0pz = findPi0.GetHighEPhotonP(2);
  m_HighEPi0energy = findPi0.GetHighEPhotonP(3);
  m_LowEPi0px = findPi0.GetLowEPhotonP(0);
  m_LowEPi0py = findPi0.GetLowEPhotonP(1);
  m_LowEPi0pz = findPi0.GetLowEPhotonP(2);
  m_LowEPi0energy = findPi0.GetLowEPhotonP(3);
  m_HighEPi0Constrainedpx = findPi0.GetHighEPhotonPConstrained(0);
  m_HighEPi0Constrainedpy = findPi0.GetHighEPhotonPConstrained(1);
  m_HighEPi0Constrainedpz = findPi0.GetHighEPhotonPConstrained(2);
  m_HighEPi0Constrainedenergy = findPi0.GetHighEPhotonPConstrained(3);
  m_LowEPi0Constrainedpx = findPi0.GetLowEPhotonPConstrained(0);
  m_LowEPi0Constrainedpy = findPi0.GetLowEPhotonPConstrained(1);
  m_LowEPi0Constrainedpz = findPi0.GetLowEPhotonPConstrained(2);
  m_LowEPi0Constrainedenergy = findPi0.GetLowEPhotonPConstrained(3);
  m_Pi0Chi2Fit = findPi0.GetChi2Fit();
  return StatusCode::SUCCESS;
}
