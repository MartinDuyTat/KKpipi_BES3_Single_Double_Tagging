// Martin Duy Tat 12th March 2021

// KKpipi
#include "KKpipi/KSpi0SingleTag.h"
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

KSpi0SingleTag::KSpi0SingleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KSpi0SingleTag::~KSpi0SingleTag() {
}

StatusCode KSpi0SingleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KSpi0 Single Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KSpi0SingleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KSpi0SingleTag", CLID_ColumnWiseTuple, "Single tagged D->KSpi0 events");
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
      status = m_tuple->addItem("KSFitSuccess", m_KSFitSuccess);
      status = m_tuple->addItem("KSDecayLengthVeeVertex", m_DecayLengthVeeVertex);
      status = m_tuple->addItem("KSChi2VeeVertex", m_Chi2VeeVertex);
      status = m_tuple->addItem("KSMassVeeVertex", m_KSMassVeeVertex);
      status = m_tuple->addItem("KSDecayLengthFit", m_DecayLengthFit);
      status = m_tuple->addItem("KSDecayLengthErrorFit", m_DecayLengthErrorFit);
      status = m_tuple->addItem("KSChi2Fit", m_Chi2Fit);
      status = m_tuple->addItem("KSPiPluspx", m_KSPiPluspx);
      status = m_tuple->addItem("KSPiPluspy", m_KSPiPluspy);
      status = m_tuple->addItem("KSPiPluspz", m_KSPiPluspz);
      status = m_tuple->addItem("KSPiPlusenergy", m_KSPiPlusenergy);
      status = m_tuple->addItem("KSPiMinuspx", m_KSPiMinuspx);
      status = m_tuple->addItem("KSPiMinuspy", m_KSPiMinuspy);
      status = m_tuple->addItem("KSPiMinuspz", m_KSPiMinuspz);
      status = m_tuple->addItem("KSPiMinusenergy", m_KSPiMinusenergy);
      status = m_tuple->addItem("HighEPi0px", m_HighEPi0px);
      status = m_tuple->addItem("HighEPi0py", m_HighEPi0py);
      status = m_tuple->addItem("HighEPi0pz", m_HighEPi0pz);
      status = m_tuple->addItem("HighEPi0energy", m_HighEPi0energy);
      status = m_tuple->addItem("LowEPi0px", m_LowEPi0px);
      status = m_tuple->addItem("LowEPi0py", m_LowEPi0py);
      status = m_tuple->addItem("LowEPi0pz", m_LowEPi0pz);
      status = m_tuple->addItem("LowEPi0energy", m_LowEPi0energy);
      status = m_tuple->addItem("Mgammagamma", m_Mgammagamma);
      status = m_tuple->addItem("HighEPi0Constrainedpx", m_HighEPi0Constrainedpx);
      status = m_tuple->addItem("HighEPi0Constrainedpy", m_HighEPi0Constrainedpy);
      status = m_tuple->addItem("HighEPi0Constrainedpz", m_HighEPi0Constrainedpz);
      status = m_tuple->addItem("HighEPi0Constrainedenergy", m_HighEPi0Constrainedenergy);
      status = m_tuple->addItem("LowEPi0Constrainedpx", m_LowEPi0Constrainedpx);
      status = m_tuple->addItem("LowEPi0Constrainedpy", m_LowEPi0Constrainedpy);
      status = m_tuple->addItem("LowEPi0Constrainedpz", m_LowEPi0Constrainedpz);
      status = m_tuple->addItem("LowEPi0Constrainedenergy", m_LowEPi0Constrainedenergy);
      status = m_tuple->addItem("Pi0Chi2Fit", m_Pi0Chi2Fit);
      status = m_tuple->addItem("IsSameDMother", m_IsSameDMother);
      status = m_tuple->addItem("IsSameDMotherAll", m_IsSameDMotherAll);
      status = m_tuple->addItem("PIDTrue", m_PIDTrue);
      status = m_tuple->addItem("KSPiPlusTrueID", m_KSPiPlusTrueID);
      status = m_tuple->addItem("KSPiMinusTrueID", m_KSPiMinusTrueID);
      status = m_tuple->addItem("HighEPi0PhotonTrueID", m_HighEPi0PhotonTrueID);
      status = m_tuple->addItem("LowEPi0PhotonTrueID", m_LowEPi0PhotonTrueID);
      status = m_tuple->addItem("KSPiPlusMotherTrueID", m_KSPiPlusMotherTrueID);
      status = m_tuple->addItem("KSPiMinusMotherTrueID", m_KSPiMinusMotherTrueID);
      status = m_tuple->addItem("HighEPi0PhotonMotherTrueID", m_HighEPi0PhotonMotherTrueID);
      status = m_tuple->addItem("LowEPi0PhotonMotherTrueID", m_LowEPi0PhotonMotherTrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KSpi0 Single Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KSpi0SingleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KSpi0 Single Tag Algorithm" << endreq;
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

StatusCode KSpi0SingleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KSpi0 Single Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KSpi0SingleTag::FillTuple(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
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
  if(status == StatusCode::SUCCESS) {
    m_KSFitSuccess = 1;
    m_DecayLengthFit = findKS.GetDecayLengthFit();
    m_DecayLengthErrorFit = findKS.GetDecayLengthErrorFit();
    m_Chi2Fit = findKS.GetChi2Fit();
  } else {
    m_KSFitSuccess = 0;
  }
  m_DecayLengthVeeVertex = findKS.GetDecayLengthVeeVertex();
  m_Chi2VeeVertex = findKS.GetChi2VeeVertex();
  m_KSMassVeeVertex = findKS.GetKSMassVeeVertex();
  m_KSPiPluspx = findKS.GetKSPiPlusP(0);
  m_KSPiPluspy = findKS.GetKSPiPlusP(1);
  m_KSPiPluspz = findKS.GetKSPiPlusP(2);
  m_KSPiPlusenergy = findKS.GetKSPiPlusP(3);
  m_KSPiMinuspx = findKS.GetKSPiMinusP(0);
  m_KSPiMinuspy = findKS.GetKSPiMinusP(1);
  m_KSPiMinuspz = findKS.GetKSPiMinusP(2);
  m_KSPiMinusenergy = findKS.GetKSPiMinusP(3);
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
  m_Mgammagamma = findPi0.GetMgammagamma();
  m_HighEPi0Constrainedpx = findPi0.GetHighEPhotonPConstrained(0);
  m_HighEPi0Constrainedpy = findPi0.GetHighEPhotonPConstrained(1);
  m_HighEPi0Constrainedpz = findPi0.GetHighEPhotonPConstrained(2);
  m_HighEPi0Constrainedenergy = findPi0.GetHighEPhotonPConstrained(3);
  m_LowEPi0Constrainedpx = findPi0.GetLowEPhotonPConstrained(0);
  m_LowEPi0Constrainedpy = findPi0.GetLowEPhotonPConstrained(1);
  m_LowEPi0Constrainedpz = findPi0.GetLowEPhotonPConstrained(2);
  m_LowEPi0Constrainedenergy = findPi0.GetLowEPhotonPConstrained(3);
  m_Pi0Chi2Fit = findPi0.GetChi2Fit();
  if(m_RunNumber < 0) {
    std::vector<int> DaughterTrackIDs = findKS.GetDaughterTrackIDs();
    PIDTruth PID_Truth(DaughterTrackIDs, 2, this);
    m_IsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    std::vector<std::pair<int, int>> PhotonPairTrackID;
    PhotonPairTrackID.push_back(std::make_pair(findPi0.GetHighEPhotonTrackID(), findPi0.GetLowEPhotonTrackID()));
    PID_Truth = PIDTruth(DaughterTrackIDs, 2, this, PhotonPairTrackID);
    m_IsSameDMotherAll = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {211, -211, 0, 0};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_PIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_KSPiPlusTrueID = ReconstructedPID[0];
    m_KSPiMinusTrueID = ReconstructedPID[1];
    m_HighEPi0PhotonTrueID = ReconstructedPID[2];
    m_LowEPi0PhotonTrueID = ReconstructedPID[3];
    m_KSPiPlusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[0], true);
    m_KSPiMinusMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[1], true);
    m_HighEPi0PhotonMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[2], false);
    m_LowEPi0PhotonMotherTrueID = PID_Truth.GetTrueMotherID(DaughterTrackIDs[3], false);
  }
  return StatusCode::SUCCESS;
}
