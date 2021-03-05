// Martin Duy Tat 4th March 2021

// KKpipi
#include "KKpipi/KKpipiVersusKSetaPrimerhogammaDoubleTag.h"
#include "KKpipi/FindKKpipiTagInfo.h"
#include "KKpipi/FindhhTagInfo.h"
#include "KKpipi/FindKS.h"
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
// ROOT
#include "TMath.h"
// STL
#include<vector>
#include<string>

KKpipiVersusKSetaPrimerhogammaDoubleTag::KKpipiVersusKSetaPrimerhogammaDoubleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKpipiVersusKSetaPrimerhogammaDoubleTag::~KKpipiVersusKSetaPrimerhogammaDoubleTag() {
}

StatusCode KKpipiVersusKSetaPrimerhogammaDoubleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KKpipi vs KSetaPrime(rhogamma) Double Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KSetaPrimerhogammaDoubleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KSetaPrimerhogammaDoubleTag", CLID_ColumnWiseTuple, "Double tagged D->KKpipi vs D->KSetaPrime(rhogamma) events");
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
      status = m_tuple->addItem("SignalKSFitSuccess", m_SignalKSFitSuccess);
      status = m_tuple->addItem("SignalKSDecayLengthVeeVertex", m_SignalDecayLengthVeeVertex);
      status = m_tuple->addItem("SignalKSChi2VeeVertex", m_SignalChi2VeeVertex);
      status = m_tuple->addItem("SignalKSMassVeeVertex", m_SignalKSMassVeeVertex);
      status = m_tuple->addItem("SignalKSDecayLengthFit", m_SignalDecayLengthFit);
      status = m_tuple->addItem("SignalKSDecayLengthErrorFit", m_SignalDecayLengthErrorFit);
      status = m_tuple->addItem("SignalKSChi2Fit", m_SignalChi2Fit);
      status = m_tuple->addItem("SignalKSMassFit", m_SignalKSMassFit);
      status = m_tuple->addItem("TagKSDecayLengthVeeVertex", m_TagDecayLengthVeeVertex);
      status = m_tuple->addItem("TagKSChi2VeeVertex", m_TagChi2VeeVertex);
      status = m_tuple->addItem("TagKSMassVeeVertex", m_TagKSMassVeeVertex);
      status = m_tuple->addItem("TagKSDecayLengthFit", m_TagDecayLengthFit);
      status = m_tuple->addItem("TagKSDecayLengthErrorFit", m_TagDecayLengthErrorFit);
      status = m_tuple->addItem("TagKSChi2Fit", m_TagChi2Fit);
      status = m_tuple->addItem("TagKSMassFit", m_TagKSMassFit);
      status = m_tuple->addItem("TagKSPiPluspx", m_TagKSPiPluspx);
      status = m_tuple->addItem("TagKSPiPluspy", m_TagKSPiPluspy);
      status = m_tuple->addItem("TagKSPiPluspz", m_TagKSPiPluspz);
      status = m_tuple->addItem("TagKSPiPlusenergy", m_TagKSPiPlusenergy);
      status = m_tuple->addItem("TagKSPiMinuspx", m_TagKSPiMinuspx);
      status = m_tuple->addItem("TagKSPiMinuspy", m_TagKSPiMinuspy);
      status = m_tuple->addItem("TagKSPiMinuspz", m_TagKSPiMinuspz);
      status = m_tuple->addItem("TagKSPiMinusenergy", m_TagKSPiMinusenergy);
      status = m_tuple->addItem("TagKSPiPluspxFit", m_TagKSPiPluspxFit);
      status = m_tuple->addItem("TagKSPiPluspyFit", m_TagKSPiPluspyFit);
      status = m_tuple->addItem("TagKSPiPluspzFit", m_TagKSPiPluspzFit);
      status = m_tuple->addItem("TagKSPiPlusenergyFit", m_TagKSPiPlusenergyFit);
      status = m_tuple->addItem("TagKSPiMinuspxFit", m_TagKSPiMinuspxFit);
      status = m_tuple->addItem("TagKSPiMinuspyFit", m_TagKSPiMinuspyFit);
      status = m_tuple->addItem("TagKSPiMinuspzFit", m_TagKSPiMinuspzFit);
      status = m_tuple->addItem("TagKSPiMinusenergyFit", m_TagKSPiMinusenergyFit);
      status = m_tuple->addItem("TagpipiKSFitSuccess", m_TagpipiKSFitSuccess);
      status = m_tuple->addItem("TagpipiKSDecayLengthVeeVertex", m_TagpipiDecayLengthVeeVertex);
      status = m_tuple->addItem("TagpipiKSChi2VeeVertex", m_TagpipiChi2VeeVertex);
      status = m_tuple->addItem("TagpipiKSMassVeeVertex", m_TagpipiKSMassVeeVertex);
      status = m_tuple->addItem("TagpipiKSDecayLengthFit", m_TagpipiDecayLengthFit);
      status = m_tuple->addItem("TagpipiKSDecayLengthErrorFit", m_TagpipiDecayLengthErrorFit);
      status = m_tuple->addItem("TagpipiKSChi2Fit", m_TagpipiChi2Fit);
      status = m_tuple->addItem("TagpipiKSMassFit", m_TagpipiKSMassFit);
      status = m_tuple->addItem("TagPiPluspx", m_TagPiPluspx);
      status = m_tuple->addItem("TagPiPluspy", m_TagPiPluspy);
      status = m_tuple->addItem("TagPiPluspz", m_TagPiPluspz);
      status = m_tuple->addItem("TagPiPlusenergy", m_TagPiPlusenergy);
      status = m_tuple->addItem("TagPiMinuspx", m_TagPiMinuspx);
      status = m_tuple->addItem("TagPiMinuspy", m_TagPiMinuspy);
      status = m_tuple->addItem("TagPiMinuspz", m_TagPiMinuspz);
      status = m_tuple->addItem("TagPiMinusenergy", m_TagPiMinusenergy);
      status = m_tuple->addItem("TagGammapx", m_TagGammapx);
      status = m_tuple->addItem("TagGammapy", m_TagGammapy);
      status = m_tuple->addItem("TagGammapz", m_TagGammapz);
      status = m_tuple->addItem("TagGammaenergy", m_TagGammaenergy);
      status = m_tuple->addItem("TagNumberShowers", m_NumberShowers);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KKpipi vs KSetaPrime(rhogamma) Double Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KKpipiVersusKSetaPrimerhogammaDoubleTag::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KKpipi vs KSetaPrime(rhogamma) Double Tag Algorithm" << endreq;
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
  if(DTTool.findDTag(EvtRecDTag::kD0toKKPiPi, EvtRecDTag::kD0toKsEPRhoGam)) {
    DTagToolIterator DTTool_Signal_iter = DTTool.dtag1();
    DTagToolIterator DTTool_Tag_iter = DTTool.dtag2();
    StatusCode FillTupleStatus = FillTuple(DTTool_Signal_iter, DTTool_Tag_iter, DTTool);
    if(FillTupleStatus != StatusCode::SUCCESS) {
        log << MSG::FATAL << "Assigning tuple info failed" << endreq;
	return StatusCode::FAILURE;
    }
    m_tuple->write();
  }
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKSetaPrimerhogammaDoubleTag::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KKpipi vs KSetaPrime(rhogamma) Double Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKSetaPrimerhogammaDoubleTag::FillTuple(DTagToolIterator DTTool_Signal_iter, DTagToolIterator DTTool_Tag_iter, DTagTool &DTTool) {
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
  m_SignalKSFitSuccess = findKKpipiTagInfo.GetKSFitSuccess();
  m_SignalDecayLengthVeeVertex = findKKpipiTagInfo.GetDecayLengthVeeVertex();
  m_SignalChi2VeeVertex = findKKpipiTagInfo.GetChi2VeeVertex();
  m_SignalKSMassVeeVertex = findKKpipiTagInfo.GetKSMassVeeVertex();
  m_SignalDecayLengthFit = findKKpipiTagInfo.GetDecayLengthFit();
  m_SignalDecayLengthErrorFit = findKKpipiTagInfo.GetDecayLengthErrorFit();
  m_SignalChi2Fit = findKKpipiTagInfo.GetChi2Fit();
  m_SignalKSMassFit = findKKpipiTagInfo.GetKSMassFit();
  FindKS findKS(true);
  status = findKS.findKS(DTTool_Tag_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_TagDecayLengthVeeVertex = findKS.GetDecayLengthVeeVertex();
  m_TagChi2VeeVertex = findKS.GetChi2VeeVertex();
  m_TagKSMassVeeVertex = findKS.GetKSMassVeeVertex();
  m_TagDecayLengthFit = findKS.GetDecayLengthFit();
  m_TagDecayLengthErrorFit = findKS.GetDecayLengthErrorFit();
  m_TagChi2Fit = findKS.GetChi2Fit();
  m_TagKSMassFit = findKS.GetKSMassFit();
  m_TagKSPiPluspx = findKS.GetKSPiPlusP(0);
  m_TagKSPiPluspy = findKS.GetKSPiPlusP(1);
  m_TagKSPiPluspz = findKS.GetKSPiPlusP(2);
  m_TagKSPiPlusenergy = findKS.GetKSPiPlusP(3);
  m_TagKSPiMinuspx = findKS.GetKSPiMinusP(0);
  m_TagKSPiMinuspy = findKS.GetKSPiMinusP(1);
  m_TagKSPiMinuspz = findKS.GetKSPiMinusP(2);
  m_TagKSPiMinusenergy = findKS.GetKSPiMinusP(3);
  m_TagKSPiPluspxFit = findKS.GetKSPiPlusPFit(0);
  m_TagKSPiPluspyFit = findKS.GetKSPiPlusPFit(1);
  m_TagKSPiPluspzFit = findKS.GetKSPiPlusPFit(2);
  m_TagKSPiPlusenergyFit = findKS.GetKSPiPlusPFit(3);
  m_TagKSPiMinuspxFit = findKS.GetKSPiMinusPFit(0);
  m_TagKSPiMinuspyFit = findKS.GetKSPiMinusPFit(1);
  m_TagKSPiMinuspzFit = findKS.GetKSPiMinusPFit(2);
  m_TagKSPiMinusenergyFit = findKS.GetKSPiMinusPFit(3);
  FindhhTagInfo findpipiTagInfo("pipi", findKS.GetDaughterTrackIDs());
  status = findpipiTagInfo.CalculateTagInfo(DTTool_Tag_iter, DTTool);
  if(status != StatusCode::SUCCESS) {
    return status;
  }
  m_TagPiPluspx = findpipiTagInfo.GethPlusP(0);
  m_TagPiPluspy = findpipiTagInfo.GethPlusP(1);
  m_TagPiPluspz = findpipiTagInfo.GethPlusP(2);
  m_TagPiPlusenergy = findpipiTagInfo.GethPlusP(3);
  m_TagPiMinuspx = findpipiTagInfo.GethMinusP(0);
  m_TagPiMinuspy = findpipiTagInfo.GethMinusP(1);
  m_TagPiMinuspz = findpipiTagInfo.GethMinusP(2);
  m_TagPiMinusenergy = findpipiTagInfo.GethMinusP(3);
  std::vector<int> VetoKSIDs;
  VetoKSIDs.push_back(DTTool.ksId(DTTool_Tag_iter)[0]);
  FindKS findKSFromPiPi(false, VetoKSIDs);
  StatusCode statuscode = findKSFromPiPi.findKS(DTTool_Tag_iter, DTTool);
  m_TagpipiKSFitSuccess = 0;
  if(statuscode == StatusCode::SUCCESS) {
    m_TagpipiKSFitSuccess = 1;
    m_TagpipiDecayLengthVeeVertex = findKSFromPiPi.GetDecayLengthVeeVertex();
    m_TagpipiChi2VeeVertex = findKSFromPiPi.GetChi2VeeVertex();
    m_TagpipiKSMassVeeVertex = findKSFromPiPi.GetKSMassVeeVertex();
    m_TagpipiDecayLengthFit = findKSFromPiPi.GetDecayLengthFit();
    m_TagpipiDecayLengthErrorFit = findKSFromPiPi.GetDecayLengthErrorFit();
    m_TagpipiChi2Fit = findKSFromPiPi.GetChi2Fit();
    m_TagpipiKSMassFit = findKSFromPiPi.GetKSMassFit();
  }
  // Get photon momentum
  SmartRefVector<EvtRecTrack> Showers = (*DTTool_Tag_iter)->showers();
  m_NumberShowers = Showers.size();
  RecEmcShower *PhotonShower = Showers[0]->emcShower();
  double GammaEnergy = PhotonShower->energy();
  double GammaTheta = PhotonShower->theta();
  double GammaPhi = PhotonShower->phi();
  m_TagGammapx = GammaEnergy*TMath::Sin(GammaTheta)*TMath::Cos(GammaPhi);
  m_TagGammapy = GammaEnergy*TMath::Sin(GammaTheta)*TMath::Sin(GammaPhi);
  m_TagGammapz = GammaEnergy*TMath::Cos(GammaTheta);
  m_TagGammaenergy = GammaEnergy;
  return StatusCode::SUCCESS;
}