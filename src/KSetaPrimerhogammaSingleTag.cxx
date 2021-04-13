// Martin Duy Tat 12th March 2021

// KKpipi
#include "KKpipi/KSetaPrimerhogammaSingleTag.h"
#include "KKpipi/FindKS.h"
#include "KKpipi/FindhhTagInfo.h"
#include "KKpipi/FindMCInfo.h"
#include "KKpipi/PIDTruth.h"
#include "KKpipi/ParticleMasses.h"
#include "KKpipi/KKpipiUtilities.h"
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

KSetaPrimerhogammaSingleTag::KSetaPrimerhogammaSingleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KSetaPrimerhogammaSingleTag::~KSetaPrimerhogammaSingleTag() {
}

StatusCode KSetaPrimerhogammaSingleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KSetaPrime(rhogamma) Single Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KSetaPrimerhogammaSingleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KSetaPrimerhogammaSingleTag", CLID_ColumnWiseTuple, "Single tagged D->KSetaPrime(rhogamma) events");
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
      status = m_tuple->addItem("pipiKSFitSuccess", m_pipiKSFitSuccess);
      status = m_tuple->addItem("pipiKSDecayLengthVeeVertex", m_pipiDecayLengthVeeVertex);
      status = m_tuple->addItem("pipiKSChi2VeeVertex", m_pipiChi2VeeVertex);
      status = m_tuple->addItem("pipiKSMassVeeVertex", m_pipiKSMassVeeVertex);
      status = m_tuple->addItem("pipiKSDecayLengthFit", m_pipiDecayLengthFit);
      status = m_tuple->addItem("pipiKSDecayLengthErrorFit", m_pipiDecayLengthErrorFit);
      status = m_tuple->addItem("pipiKSChi2Fit", m_pipiChi2Fit);
      status = m_tuple->addItem("pipiKSMassFit", m_pipiKSMassFit);
      status = m_tuple->addItem("PiPluspx", m_PiPluspx);
      status = m_tuple->addItem("PiPluspy", m_PiPluspy);
      status = m_tuple->addItem("PiPluspz", m_PiPluspz);
      status = m_tuple->addItem("PiPlusenergy", m_PiPlusenergy);
      status = m_tuple->addItem("PiMinuspx", m_PiMinuspx);
      status = m_tuple->addItem("PiMinuspy", m_PiMinuspy);
      status = m_tuple->addItem("PiMinuspz", m_PiMinuspz);
      status = m_tuple->addItem("PiMinusenergy", m_PiMinusenergy);
      status = m_tuple->addItem("Gammapx", m_Gammapx);
      status = m_tuple->addItem("Gammapy", m_Gammapy);
      status = m_tuple->addItem("Gammapz", m_Gammapz);
      status = m_tuple->addItem("Gammaenergy", m_Gammaenergy);
      status = m_tuple->addItem("Mpipigamma", m_Mpipigamma);
      status = m_tuple->addItem("PhotonAngleSeparation", m_PhotonAngleSeparation);
      status = m_tuple->addItem("PhotonThetaSeparation", m_PhotonThetaSeparation);
      status = m_tuple->addItem("PhotonPhiSeparation", m_PhotonPhiSeparation);
      status = m_tuple->addItem("NumberShowers", m_NumberShowers);
      status = m_tuple->addItem("IsSameDMother", m_IsSameDMother);
      status = m_tuple->addItem("PIDTrue", m_PIDTrue);
      status = m_tuple->addItem("KSPiPlusTrueID", m_KSPiPlusTrueID);
      status = m_tuple->addItem("KSPiMinusTrueID", m_KSPiMinusTrueID);
      status = m_tuple->addItem("EtaPPiPlusTrueID", m_EtaPPiPlusTrueID);
      status = m_tuple->addItem("EtaPPiMinusTrueID", m_EtaPPiMinusTrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KSetaPrime(rhogamma) Single Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KSetaPrimerhogammaSingleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KSetaPrime(rhogamma) Single Tag Algorithm" << endreq;
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
  if(DTTool.findSTag(EvtRecDTag::kD0toKsEPRhoGam)) {
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(), "/Event/EventHeader");
    m_RunNumber = eventHeader->runNumber();
    m_EventNumber = eventHeader->eventNumber();
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

StatusCode KSetaPrimerhogammaSingleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KSetaPrime(rhogamma) Single Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KSetaPrimerhogammaSingleTag::FillTuple(DTagToolIterator DTTool_iter, DTagTool &DTTool) {
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
  FindhhTagInfo findpipiInfo("pipi", findKS.GetDaughterTrackIDs());
  status = findpipiInfo.CalculateTagInfo(DTTool_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_PiPluspx = findpipiInfo.GethPlusP(0);
  m_PiPluspy = findpipiInfo.GethPlusP(1);
  m_PiPluspz = findpipiInfo.GethPlusP(2);
  m_PiPlusenergy = findpipiInfo.GethPlusP(3);
  m_PiMinuspx = findpipiInfo.GethMinusP(0);
  m_PiMinuspy = findpipiInfo.GethMinusP(1);
  m_PiMinuspz = findpipiInfo.GethMinusP(2);
  m_PiMinusenergy = findpipiInfo.GethMinusP(3);
  double Mpipi = TMath::Sqrt(TMath::Power(m_PiPlusenergy + m_PiMinusenergy, 2) - TMath::Power(m_PiPluspx + m_PiMinuspx, 2) - TMath::Power(m_PiPluspy + m_PiMinuspy, 2) - TMath::Power(m_PiPluspz + m_PiMinuspz, 2));
  m_pipiKSFitSuccess = 0;
  if(Mpipi - MASS::KS_MASS < 0.050 && Mpipi - MASS::KS_MASS > -0.060) {
    FindKS findKSFromPiPi(false);
    std::vector<int> PionTrackIDs;
    PionTrackIDs.push_back(findpipiInfo.GetPiPlusTrackID());
    PionTrackIDs.push_back(findpipiInfo.GetPiMinusTrackID());
    StatusCode statuscode = findKSFromPiPi.findKS(DTTool_iter, DTTool, PionTrackIDs);
    if(statuscode == StatusCode::SUCCESS) {
      m_pipiKSFitSuccess = 1;
      m_pipiDecayLengthVeeVertex = findKSFromPiPi.GetDecayLengthVeeVertex();
      m_pipiChi2VeeVertex = findKSFromPiPi.GetChi2VeeVertex();
      m_pipiKSMassVeeVertex = findKSFromPiPi.GetKSMassVeeVertex();
      m_pipiDecayLengthFit = findKSFromPiPi.GetDecayLengthFit();
      m_pipiDecayLengthErrorFit = findKSFromPiPi.GetDecayLengthErrorFit();
      m_pipiChi2Fit = findKSFromPiPi.GetChi2Fit();
      m_pipiKSMassFit = findKSFromPiPi.GetKSMassFit();
    }
  }
  // Get photon momentum
  SmartRefVector<EvtRecTrack> Showers = (*DTTool_iter)->showers();
  m_NumberShowers = Showers.size();
  RecEmcShower *PhotonShower = Showers[0]->emcShower();
  // Get EMC position of shower
  CLHEP::Hep3Vector EMCPosition(PhotonShower->x(), PhotonShower->y(), PhotonShower->z());
  // Find separation to nearest charged track
  double Angle, Theta, Phi;
  KKpipiUtilities::GetPhotonAngularSeparation(EMCPosition, Angle, Theta, Phi);
  m_PhotonAngleSeparation = Angle;
  m_PhotonThetaSeparation = Theta;
  m_PhotonPhiSeparation = Phi;
  CLHEP::HepLorentzVector PhotonP = KKpipiUtilities::GetPhoton4Vector(PhotonShower->energy(), PhotonShower->theta(), PhotonShower->phi());
  m_Gammapx = PhotonP[0];
  m_Gammapy = PhotonP[1];
  m_Gammapz = PhotonP[2];
  m_Gammaenergy = PhotonP[3];
  m_Mpipigamma = TMath::Sqrt(TMath::Power(m_PiPlusenergy + m_PiMinusenergy + m_Gammaenergy, 2)
                           - TMath::Power(m_PiPluspx + m_PiMinuspx + m_Gammapx, 2)
                           - TMath::Power(m_PiPluspy + m_PiMinuspy + m_Gammapy, 2)
                           - TMath::Power(m_PiPluspz + m_PiMinuspz + m_Gammapz, 2));
  if(m_RunNumber < 0) {
    std::vector<int> DaughterTrackIDs = findKS.GetDaughterTrackIDs();
    std::vector<int> EtaPDaughterTrackIDs = findpipiInfo.GetDaughterTrackID();
    DaughterTrackIDs.insert(DaughterTrackIDs.end(), EtaPDaughterTrackIDs.begin(), EtaPDaughterTrackIDs.end());
    PIDTruth PID_Truth(DaughterTrackIDs, this);
    m_IsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {211, -211, 211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_PIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_KSPiPlusTrueID = ReconstructedPID[0];
    m_KSPiMinusTrueID = ReconstructedPID[1];
    m_EtaPPiPlusTrueID = ReconstructedPID[2];
    m_EtaPPiMinusTrueID = ReconstructedPID[3];
  }
  return StatusCode::SUCCESS;
}
