// Martin Duy Tat 25th March 2021

// KKpipi
#include "KKpipi/KpipipiSingleTag.h"
#include "KKpipi/FindKpipipiTagInfo.h"
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

KpipipiSingleTag::KpipipiSingleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KpipipiSingleTag::~KpipipiSingleTag() {
}

StatusCode KpipipiSingleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing Kpipipi Single Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KpipipiSingleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KpipipiSingleTag", CLID_ColumnWiseTuple, "Single tagged D->Kpipipi events");
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
      status = m_tuple->addItem("Kpx", m_Kpx);
      status = m_tuple->addItem("Kpy", m_Kpy);
      status = m_tuple->addItem("Kpz", m_Kpz);
      status = m_tuple->addItem("Kenergy", m_Kenergy);
      status = m_tuple->addItem("Pi1px", m_Pi1px);
      status = m_tuple->addItem("Pi1py", m_Pi1py);
      status = m_tuple->addItem("Pi1pz", m_Pi1pz);
      status = m_tuple->addItem("Pi1energy", m_Pi1energy);
      status = m_tuple->addItem("Pi2px", m_Pi2px);
      status = m_tuple->addItem("Pi2py", m_Pi2py);
      status = m_tuple->addItem("Pi2pz", m_Pi2pz);
      status = m_tuple->addItem("Pi2energy", m_Pi2energy);
      status = m_tuple->addItem("Pi3px", m_Pi3px);
      status = m_tuple->addItem("Pi3py", m_Pi3py);
      status = m_tuple->addItem("Pi3pz", m_Pi3pz);
      status = m_tuple->addItem("Pi3energy", m_Pi3energy);
      status = m_tuple->addItem("KCharge", m_KCharge);
      status = m_tuple->addItem("Pi1Charge", m_Pi1Charge);
      status = m_tuple->addItem("Pi2Charge", m_Pi2Charge);
      status = m_tuple->addItem("Pi3Charge", m_Pi3Charge);
      status = m_tuple->addItem("12KSFitSuccess", m_12KSFitSuccess);
      status = m_tuple->addItem("12KSDecayLengthVeeVertex", m_12DecayLengthVeeVertex);
      status = m_tuple->addItem("12KSChi2VeeVertex", m_12Chi2VeeVertex);
      status = m_tuple->addItem("12KSMassVeeVertex", m_12KSMassVeeVertex);
      status = m_tuple->addItem("12KSDecayLengthFit", m_12DecayLengthFit);
      status = m_tuple->addItem("12KSDecayLengthErrorFit", m_12DecayLengthErrorFit);
      status = m_tuple->addItem("12KSChi2Fit", m_12Chi2Fit);
      status = m_tuple->addItem("12KSMassFit", m_12KSMassFit);
      status = m_tuple->addItem("13KSFitSuccess", m_13KSFitSuccess);
      status = m_tuple->addItem("13KSDecayLengthVeeVertex", m_13DecayLengthVeeVertex);
      status = m_tuple->addItem("13KSChi2VeeVertex", m_13Chi2VeeVertex);
      status = m_tuple->addItem("13KSMassVeeVertex", m_13KSMassVeeVertex);
      status = m_tuple->addItem("13KSDecayLengthFit", m_13DecayLengthFit);
      status = m_tuple->addItem("13KSDecayLengthErrorFit", m_13DecayLengthErrorFit);
      status = m_tuple->addItem("13KSChi2Fit", m_13Chi2Fit);
      status = m_tuple->addItem("13KSMassFit", m_13KSMassFit);
      status = m_tuple->addItem("IsSameDMother", m_IsSameDMother);
      status = m_tuple->addItem("PIDTrue", m_PIDTrue);
      status = m_tuple->addItem("KTrueID", m_KTrueID);
      status = m_tuple->addItem("Pi1TrueID", m_Pi1TrueID);
      status = m_tuple->addItem("Pi2TrueID", m_Pi2TrueID);
      status = m_tuple->addItem("Pi3TrueID", m_Pi3TrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for Kpipipi Single Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KpipipiSingleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing Kpipipi Single Tag Algorithm" << endreq;
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
  if(DTTool.findSTag(EvtRecDTag::kD0toKPiPiPi)) {
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

StatusCode KpipipiSingleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing Kpipipi Single Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KpipipiSingleTag::FillTuple(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
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
  FindKpipipiTagInfo findKpipipiTagInfo;
  StatusCode status = findKpipipiTagInfo.CalculateTagInfo(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_Kpx = findKpipipiTagInfo.GetKP(0);
  m_Kpy = findKpipipiTagInfo.GetKP(1);
  m_Kpz = findKpipipiTagInfo.GetKP(2);
  m_Kenergy = findKpipipiTagInfo.GetKP(3);
  m_Pi1px = findKpipipiTagInfo.GetPi1P(0);
  m_Pi1py = findKpipipiTagInfo.GetPi1P(1);
  m_Pi1pz = findKpipipiTagInfo.GetPi1P(2);
  m_Pi1energy = findKpipipiTagInfo.GetPi1P(3);
  m_Pi2px = findKpipipiTagInfo.GetPi2P(0);
  m_Pi2py = findKpipipiTagInfo.GetPi2P(1);
  m_Pi2pz = findKpipipiTagInfo.GetPi2P(2);
  m_Pi2energy = findKpipipiTagInfo.GetPi2P(3);
  m_Pi3px = findKpipipiTagInfo.GetPi3P(0);
  m_Pi3py = findKpipipiTagInfo.GetPi3P(1);
  m_Pi3pz = findKpipipiTagInfo.GetPi3P(2);
  m_Pi3energy = findKpipipiTagInfo.GetPi3P(3);
  m_KCharge = findKpipipiTagInfo.GetKCharge();
  m_Pi1Charge = findKpipipiTagInfo.GetPi1Charge();
  m_Pi2Charge = findKpipipiTagInfo.GetPi2Charge();
  m_Pi3Charge = findKpipipiTagInfo.GetPi3Charge();
  m_12KSFitSuccess = findKpipipiTagInfo.GetKSFitSuccess12();
  m_12DecayLengthVeeVertex = findKpipipiTagInfo.GetDecayLengthVeeVertex12();
  m_12Chi2VeeVertex = findKpipipiTagInfo.GetChi2VeeVertex12();
  m_12KSMassVeeVertex = findKpipipiTagInfo.GetKSMassVeeVertex12();
  m_12DecayLengthFit = findKpipipiTagInfo.GetDecayLengthFit12();
  m_12DecayLengthErrorFit = findKpipipiTagInfo.GetDecayLengthErrorFit12();
  m_12Chi2Fit = findKpipipiTagInfo.GetChi2Fit12();
  m_12KSMassFit = findKpipipiTagInfo.GetKSMassFit12();
  m_13KSFitSuccess = findKpipipiTagInfo.GetKSFitSuccess13();
  m_13DecayLengthVeeVertex = findKpipipiTagInfo.GetDecayLengthVeeVertex13();
  m_13Chi2VeeVertex = findKpipipiTagInfo.GetChi2VeeVertex13();
  m_13KSMassVeeVertex = findKpipipiTagInfo.GetKSMassVeeVertex13();
  m_13DecayLengthFit = findKpipipiTagInfo.GetDecayLengthFit13();
  m_13DecayLengthErrorFit = findKpipipiTagInfo.GetDecayLengthErrorFit13();
  m_13Chi2Fit = findKpipipiTagInfo.GetChi2Fit13();
  m_13KSMassFit = findKpipipiTagInfo.GetKSMassFit13();
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKpipipiTagInfo.GetDaughterTrackID(), this);
    m_IsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {321*m_KCharge, 211*m_Pi1Charge, 211*m_Pi2Charge, 211*m_Pi3Charge};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_PIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_KTrueID = ReconstructedPID[0];
    m_Pi1TrueID = ReconstructedPID[1];
    m_Pi2TrueID = ReconstructedPID[2];
    m_Pi3TrueID = ReconstructedPID[3];
  }
  return StatusCode::SUCCESS;
}
