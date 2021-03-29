// Martin Duy Tat 26th March 2021

// KKpipi
#include "KKpipi/KKpipiVersusKLDoubleTags.h"
#include "KKpipi/FindKKpipiTagInfo.h"
#include "KKpipi/FindMCInfo.h"
#include "KKpipi/FindKL.h"
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
// ROOT
#include "TMath.h"
// CLHEP
#include "CLHEP/Vector/LorentzVector.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "McDecayModeSvc/McDecayModeSvc.h"
#include "McTruth/McParticle.h"
// STL
#include<vector>
#include<string>

KKpipiVersusKLDoubleTags::KKpipiVersusKLDoubleTags(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKpipiVersusKLDoubleTags::~KKpipiVersusKLDoubleTags() {
}

StatusCode KKpipiVersusKLDoubleTags::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KKpipi vs KLX Double Tagging" << endreq;
  StatusCode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/KLDoubleTags");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/KLDoubleTags", CLID_ColumnWiseTuple, "Double tagged D->KKpipi vs D->KLX events");
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
      status = m_tuple->addItem("SignalDpx", m_SignalDpx);
      status = m_tuple->addItem("SignalDpy", m_SignalDpy);
      status = m_tuple->addItem("SignalDpz", m_SignalDpz);
      status = m_tuple->addItem("SignalDenergy", m_SignalDenergy);
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
      status = m_tuple->addItem("SignalIsSameDMother", m_SignalIsSameDMother);
      status = m_tuple->addItem("SignalPIDTrue", m_SignalPIDTrue);
      status = m_tuple->addItem("SignalKPlusTrueID", m_SignalKPlusTrueID);
      status = m_tuple->addItem("SignalKMinusTrueID", m_SignalKMinusTrueID);
      status = m_tuple->addItem("SignalPiPlusTrueID", m_SignalPiPlusTrueID);
      status = m_tuple->addItem("SignalPiMinusTrueID", m_SignalPiMinusTrueID);
      status = m_tuple->addItem("TagFoundPionPair", m_TagFoundPionPair);
      status = m_tuple->addItem("TagPiPluspx", m_TagPiPluspx);
      status = m_tuple->addItem("TagPiPluspy", m_TagPiPluspy);
      status = m_tuple->addItem("TagPiPluspz", m_TagPiPluspz);
      status = m_tuple->addItem("TagPiPlusenergy", m_TagPiPlusenergy);
      status = m_tuple->addItem("TagPiMinuspx", m_TagPiMinuspx);
      status = m_tuple->addItem("TagPiMinuspy", m_TagPiMinuspy);
      status = m_tuple->addItem("TagPiMinuspz", m_TagPiMinuspz);
      status = m_tuple->addItem("TagPiMinusenergy", m_TagPiMinusenergy);
      status = m_tuple->addItem("TagNumberPi0", m_TagNumberPi0, 0, 100);
      status = m_tuple->addIndexedItem("TagPi0HighEPhotonpx", m_TagNumberPi0, m_TagPi0HighEPhotonpx);
      status = m_tuple->addIndexedItem("TagPi0HighEPhotonpy", m_TagNumberPi0, m_TagPi0HighEPhotonpy);
      status = m_tuple->addIndexedItem("TagPi0HighEPhotonpz", m_TagNumberPi0, m_TagPi0HighEPhotonpz);
      status = m_tuple->addIndexedItem("TagPi0HighEPhotonenergy", m_TagNumberPi0, m_TagPi0HighEPhotonenergy);
      status = m_tuple->addIndexedItem("TagPi0LowEPhotonpx", m_TagNumberPi0, m_TagPi0LowEPhotonpx);
      status = m_tuple->addIndexedItem("TagPi0LowEPhotonpy", m_TagNumberPi0, m_TagPi0LowEPhotonpy);
      status = m_tuple->addIndexedItem("TagPi0LowEPhotonpz", m_TagNumberPi0, m_TagPi0LowEPhotonpz);
      status = m_tuple->addIndexedItem("TagPi0LowEPhotonenergy", m_TagNumberPi0, m_TagPi0LowEPhotonenergy);
      status = m_tuple->addIndexedItem("TagPi0HighEPhotonpxConstrained", m_TagNumberPi0, m_TagPi0HighEPhotonpxConstrained);
      status = m_tuple->addIndexedItem("TagPi0HighEPhotonpyConstrained", m_TagNumberPi0, m_TagPi0HighEPhotonpyConstrained);
      status = m_tuple->addIndexedItem("TagPi0HighEPhotonpzConstrained", m_TagNumberPi0, m_TagPi0HighEPhotonpzConstrained);
      status = m_tuple->addIndexedItem("TagPi0HighEPhotonenergyConstrained", m_TagNumberPi0, m_TagPi0HighEPhotonenergyConstrained);
      status = m_tuple->addIndexedItem("TagPi0LowEPhotonpxConstrained", m_TagNumberPi0, m_TagPi0LowEPhotonpxConstrained);
      status = m_tuple->addIndexedItem("TagPi0LowEPhotonpyConstrained", m_TagNumberPi0, m_TagPi0LowEPhotonpyConstrained);
      status = m_tuple->addIndexedItem("TagPi0LowEPhotonpzConstrained", m_TagNumberPi0, m_TagPi0LowEPhotonpzConstrained);
      status = m_tuple->addIndexedItem("TagPi0LowEPhotonenergyConstrained", m_TagNumberPi0, m_TagPi0LowEPhotonenergyConstrained);
      status = m_tuple->addIndexedItem("TagPi0Chi2Fit", m_TagNumberPi0, m_TagPi0Chi2Fit);
      status = m_tuple->addIndexedItem("TagPi0HighEPhotonTrackID", m_TagNumberPi0, m_TagPi0HighEPhotonTrackID);
      status = m_tuple->addIndexedItem("TagPi0LowEPhotonTrackID", m_TagNumberPi0, m_TagPi0LowEPhotonTrackID);
      status = m_tuple->addItem("TagNumberEta", m_TagNumberEta, 0, 100);
      status = m_tuple->addIndexedItem("TagEtaHighEPhotonpx", m_TagNumberEta, m_TagEtaHighEPhotonpx);
      status = m_tuple->addIndexedItem("TagEtaHighEPhotonpy", m_TagNumberEta, m_TagEtaHighEPhotonpy);
      status = m_tuple->addIndexedItem("TagEtaHighEPhotonpz", m_TagNumberEta, m_TagEtaHighEPhotonpz);
      status = m_tuple->addIndexedItem("TagEtaHighEPhotonenergy", m_TagNumberEta, m_TagEtaHighEPhotonenergy);
      status = m_tuple->addIndexedItem("TagEtaLowEPhotonpx", m_TagNumberEta, m_TagEtaLowEPhotonpx);
      status = m_tuple->addIndexedItem("TagEtaLowEPhotonpy", m_TagNumberEta, m_TagEtaLowEPhotonpy);
      status = m_tuple->addIndexedItem("TagEtaLowEPhotonpz", m_TagNumberEta, m_TagEtaLowEPhotonpz);
      status = m_tuple->addIndexedItem("TagEtaLowEPhotonenergy", m_TagNumberEta, m_TagEtaLowEPhotonenergy);
      status = m_tuple->addIndexedItem("TagEtaHighEPhotonpxConstrained", m_TagNumberEta, m_TagEtaHighEPhotonpxConstrained);
      status = m_tuple->addIndexedItem("TagEtaHighEPhotonpyConstrained", m_TagNumberEta, m_TagEtaHighEPhotonpyConstrained);
      status = m_tuple->addIndexedItem("TagEtaHighEPhotonpzConstrained", m_TagNumberEta, m_TagEtaHighEPhotonpzConstrained);
      status = m_tuple->addIndexedItem("TagEtaHighEPhotonenergyConstrained", m_TagNumberEta, m_TagEtaHighEPhotonenergyConstrained);
      status = m_tuple->addIndexedItem("TagEtaLowEPhotonpxConstrained", m_TagNumberEta, m_TagEtaLowEPhotonpxConstrained);
      status = m_tuple->addIndexedItem("TagEtaLowEPhotonpyConstrained", m_TagNumberEta, m_TagEtaLowEPhotonpyConstrained);
      status = m_tuple->addIndexedItem("TagEtaLowEPhotonpzConstrained", m_TagNumberEta, m_TagEtaLowEPhotonpzConstrained);
      status = m_tuple->addIndexedItem("TagEtaLowEPhotonenergyConstrained", m_TagNumberEta, m_TagEtaLowEPhotonenergyConstrained);
      status = m_tuple->addIndexedItem("TagEtaChi2Fit", m_TagNumberEta, m_TagEtaChi2Fit);
      status = m_tuple->addIndexedItem("TagEtaHighEPhotonTrackID", m_TagNumberEta, m_TagEtaHighEPhotonTrackID);
      status = m_tuple->addIndexedItem("TagEtaLowEPhotonTrackID", m_TagNumberEta, m_TagEtaLowEPhotonTrackID);
      status = m_tuple->addItem("TagNumberGamma", m_TagNumberGamma, 0, 100);
      status = m_tuple->addIndexedItem("TagPhotonEnergy", m_TagNumberGamma, m_TagPhotonEnergy);
      status = m_tuple->addIndexedItem("TagPhotonPx", m_TagNumberGamma, m_TagPhotonPx);
      status = m_tuple->addIndexedItem("TagPhotonPy", m_TagNumberGamma, m_TagPhotonPy);
      status = m_tuple->addIndexedItem("TagPhotonPz", m_TagNumberGamma, m_TagPhotonPz);
      status = m_tuple->addIndexedItem("TagPhotonAngleSeparation", m_TagNumberGamma, m_TagPhotonAngleSeparation);
      status = m_tuple->addIndexedItem("TagPhotonThetaSeparation", m_TagNumberGamma, m_TagPhotonThetaSeparation);
      status = m_tuple->addIndexedItem("TagPhotonPhiSeparation", m_TagNumberGamma, m_TagPhotonPhiSeparation);
      status = m_tuple->addIndexedItem("TagPhotonTrackID", m_TagNumberPi0, m_TagPhotonTrackID);
      status = m_tuple->addItem("TagMissingEnergy", m_TagMissingEnergy);
      status = m_tuple->addItem("TagMissingMass2", m_TagMissingMass2);
      status = m_tuple->addItem("TagIsSameDMother", m_TagIsSameDMother);
      status = m_tuple->addItem("TagPIDTrue", m_TagPIDTrue);
      status = m_tuple->addItem("TagPiPlusTrueID", m_TagPiPlusTrueID);
      status = m_tuple->addItem("TagPiMinusTrueID", m_TagPiMinusTrueID);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KKpipi vs KLX Double Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}

StatusCode KKpipiVersusKLDoubleTags::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Executing KKpipi vs KLX Double Tag Algorithm" << endreq;
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
  if(DTTool.findSTag(EvtRecDTag::kD0toKKPiPi)) {
    DTagToolIterator DTTool_Signal_iter = DTTool.stag();
    StatusCode FillTupleStatus = FillTuple(DTTool_Signal_iter, DTTool);
    if(FillTupleStatus == StatusCode::RECOVERABLE) {
      return StatusCode::SUCCESS;
    } else if(FillTupleStatus == StatusCode::FAILURE) {
      log << MSG::FATAL << "Assigning KL tuple info failed" << endreq;
      return StatusCode::FAILURE;
    } else {
      m_tuple->write();
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKLDoubleTags::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Finalizing KKpipi vs KL Double Tagging" << endreq;
  return StatusCode::SUCCESS;
}

StatusCode KKpipiVersusKLDoubleTags::FillTuple(DTagToolIterator DTTool_Signal_iter, DTagTool &DTTool) {
  // First check if there are any KL candidates, otherwise no point in saving all the other stuff
  FindKL findKL;
  StatusCode FoundKL = findKL.findKL(DTTool_Signal_iter, DTTool);
  if(FoundKL == StatusCode::FAILURE) {
    return StatusCode::RECOVERABLE;
  }
  // Save MC truth information
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
  // Save information about signal side
  m_SignalDMass = (*DTTool_Signal_iter)->mass();
  m_SignalMBC = (*DTTool_Signal_iter)->mBC();
  m_SignalDeltaE = (*DTTool_Signal_iter)->deltaE();
  m_SignalBeamE = (*DTTool_Signal_iter)->beamE();
  m_SignalDpx = (*DTTool_Signal_iter)->p4().x();
  m_SignalDpy = (*DTTool_Signal_iter)->p4().y();
  m_SignalDpz = (*DTTool_Signal_iter)->p4().z();
  m_SignalDenergy = (*DTTool_Signal_iter)->p4().t();
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
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKKpipiTagInfo.GetDaughterTrackID(), this);
    m_SignalIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[4] = {321, -321, 211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 4);
    m_SignalPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_SignalKPlusTrueID = ReconstructedPID[0];
    m_SignalKMinusTrueID = ReconstructedPID[1];
    m_SignalPiPlusTrueID = ReconstructedPID[2];
    m_SignalPiMinusTrueID = ReconstructedPID[3];
  }
  // Save information about tag side
  m_TagFoundPionPair = findKL.GetFoundPionPair();
  if(m_TagFoundPionPair == 1) {
    m_TagPiPluspx = findKL.GetPiPlusP(0);
    m_TagPiPluspy = findKL.GetPiPlusP(1);
    m_TagPiPluspz = findKL.GetPiPlusP(2);
    m_TagPiPlusenergy = findKL.GetPiPlusP(3);
    m_TagPiMinuspx = findKL.GetPiMinusP(0);
    m_TagPiMinuspy = findKL.GetPiMinusP(1);
    m_TagPiMinuspz = findKL.GetPiMinusP(2);
    m_TagPiMinusenergy = findKL.GetPiMinusP(3);
  }
  m_TagNumberPi0 = findKL.GetNumberPi0();
  for(int j = 0; j < m_TagNumberPi0; j++) {
    m_TagPi0HighEPhotonpx[j] = findKL.GetPi0HighEPhotonP(0, j);
    m_TagPi0HighEPhotonpy[j] = findKL.GetPi0HighEPhotonP(1, j);
    m_TagPi0HighEPhotonpz[j] = findKL.GetPi0HighEPhotonP(2, j);
    m_TagPi0HighEPhotonenergy[j] = findKL.GetPi0HighEPhotonP(3, j);
    m_TagPi0LowEPhotonpx[j] = findKL.GetPi0LowEPhotonP(0, j);
    m_TagPi0LowEPhotonpy[j] = findKL.GetPi0LowEPhotonP(1, j);
    m_TagPi0LowEPhotonpz[j] = findKL.GetPi0LowEPhotonP(2, j);
    m_TagPi0LowEPhotonenergy[j] = findKL.GetPi0LowEPhotonP(3, j);
    m_TagPi0HighEPhotonpxConstrained[j] = findKL.GetPi0HighEPhotonPConstrained(0, j);
    m_TagPi0HighEPhotonpyConstrained[j] = findKL.GetPi0HighEPhotonPConstrained(1, j);
    m_TagPi0HighEPhotonpzConstrained[j] = findKL.GetPi0HighEPhotonPConstrained(2, j);
    m_TagPi0HighEPhotonenergyConstrained[j] = findKL.GetPi0HighEPhotonPConstrained(3, j);
    m_TagPi0LowEPhotonpxConstrained[j] = findKL.GetPi0LowEPhotonPConstrained(0, j);
    m_TagPi0LowEPhotonpyConstrained[j] = findKL.GetPi0LowEPhotonPConstrained(1, j);
    m_TagPi0LowEPhotonpzConstrained[j] = findKL.GetPi0LowEPhotonPConstrained(2, j);
    m_TagPi0LowEPhotonenergyConstrained[j] = findKL.GetPi0LowEPhotonPConstrained(3, j);
    m_TagPi0Chi2Fit[j] = findKL.GetPi0Chi2Fit(j);
    m_TagPi0HighEPhotonTrackID[j] = findKL.GetPi0HighEPhotonTrackID(j);
    m_TagPi0LowEPhotonTrackID[j] = findKL.GetPi0LowEPhotonTrackID(j);
  }
  m_TagNumberEta = findKL.GetNumberEta();
  for(int j = 0; j < m_TagNumberEta; j++) {
    m_TagEtaHighEPhotonpx[j] = findKL.GetEtaHighEPhotonP(0, j);
    m_TagEtaHighEPhotonpy[j] = findKL.GetEtaHighEPhotonP(1, j);
    m_TagEtaHighEPhotonpz[j] = findKL.GetEtaHighEPhotonP(2, j);
    m_TagEtaHighEPhotonenergy[j] = findKL.GetEtaHighEPhotonP(3, j);
    m_TagEtaLowEPhotonpx[j] = findKL.GetEtaLowEPhotonP(0, j);
    m_TagEtaLowEPhotonpy[j] = findKL.GetEtaLowEPhotonP(1, j);
    m_TagEtaLowEPhotonpz[j] = findKL.GetEtaLowEPhotonP(2, j);
    m_TagEtaLowEPhotonenergy[j] = findKL.GetEtaLowEPhotonP(3, j);
    m_TagEtaHighEPhotonpxConstrained[j] = findKL.GetEtaHighEPhotonPConstrained(0, j);
    m_TagEtaHighEPhotonpyConstrained[j] = findKL.GetEtaHighEPhotonPConstrained(1, j);
    m_TagEtaHighEPhotonpzConstrained[j] = findKL.GetEtaHighEPhotonPConstrained(2, j);
    m_TagEtaHighEPhotonenergyConstrained[j] = findKL.GetEtaHighEPhotonPConstrained(3, j);
    m_TagEtaLowEPhotonpxConstrained[j] = findKL.GetEtaLowEPhotonPConstrained(0, j);
    m_TagEtaLowEPhotonpyConstrained[j] = findKL.GetEtaLowEPhotonPConstrained(1, j);
    m_TagEtaLowEPhotonpzConstrained[j] = findKL.GetEtaLowEPhotonPConstrained(2, j);
    m_TagEtaLowEPhotonenergyConstrained[j] = findKL.GetEtaLowEPhotonPConstrained(3, j);
    m_TagEtaChi2Fit[j] = findKL.GetEtaChi2Fit(j);
    m_TagEtaHighEPhotonTrackID[j] = findKL.GetEtaHighEPhotonTrackID(j);
    m_TagEtaLowEPhotonTrackID[j] = findKL.GetEtaLowEPhotonTrackID(j);
  }
  m_TagNumberGamma = findKL.GetNumberGamma();
  for(int j = 0; j < m_TagNumberGamma; j++) {
    m_TagPhotonPx[j] = findKL.GetPhotonEnergy(0, j);
    m_TagPhotonPy[j] = findKL.GetPhotonEnergy(1, j);
    m_TagPhotonPz[j] = findKL.GetPhotonEnergy(2, j);
    m_TagPhotonEnergy[j] = findKL.GetPhotonEnergy(3, j);
    m_TagPhotonAngleSeparation[j] = findKL.GetPhotonAngleSeparation(j);
    m_TagPhotonThetaSeparation[j] = findKL.GetPhotonThetaSeparation(j);
    m_TagPhotonPhiSeparation[j] = findKL.GetPhotonPhiSeparation(j);
    m_TagPhotonTrackID[j] = findKL.GetPhotonTrackID(j);
  }
  FillMissingMassEnergy();
  if(m_RunNumber < 0) {
    PIDTruth PID_Truth(findKL.GetDaughterTrackID(), this);
    m_TagIsSameDMother = PID_Truth.SameDMother() ? 1 : 0;
    int SomeArray[2] = {211, -211};
    std::vector<int> ReconstructedPID(SomeArray, SomeArray + 2);
    m_TagPIDTrue = PID_Truth.FindTrueID(ReconstructedPID) ? 1 : 0;
    m_TagPiPlusTrueID = ReconstructedPID[0];
    m_TagPiMinusTrueID = ReconstructedPID[1];
  }
  return StatusCode::SUCCESS;
}

void KKpipiVersusKLDoubleTags::FillMissingMassEnergy() {
  double Pi0Energy = 0.0, Pi0Px = 0.0, Pi0Py = 0.0, Pi0Pz = 0.0;
  for(int j = 0; j < m_TagNumberPi0; j++) {
    Pi0Energy += m_TagPi0HighEPhotonenergy[j] + m_TagPi0LowEPhotonenergy[j];
    Pi0Px += m_TagPi0HighEPhotonpx[j] + m_TagPi0LowEPhotonpx[j];
    Pi0Py += m_TagPi0HighEPhotonpy[j] + m_TagPi0LowEPhotonpy[j];
    Pi0Pz += m_TagPi0HighEPhotonpz[j] + m_TagPi0LowEPhotonpz[j];
  }
  double EtaEnergy = 0.0, EtaPx = 0.0, EtaPy = 0.0, EtaPz = 0.0;
  for(int j = 0; j < m_TagNumberEta; j++) {
    EtaEnergy += m_TagEtaHighEPhotonenergy[j] + m_TagEtaLowEPhotonenergy[j];
    EtaPx += m_TagEtaHighEPhotonpx[j] + m_TagEtaLowEPhotonpx[j];
    EtaPy += m_TagEtaHighEPhotonpy[j] + m_TagEtaLowEPhotonpy[j];
    EtaPz += m_TagEtaHighEPhotonpz[j] + m_TagEtaLowEPhotonpz[j];
  }
  double EMiss = m_SignalBeamE - Pi0Energy - EtaEnergy;
  double PxMiss = -m_SignalDpx - Pi0Px - EtaPx;
  double PyMiss = -m_SignalDpy - Pi0Py - EtaPy;
  double PzMiss = -m_SignalDpz - Pi0Pz - EtaPz;
  if(m_TagFoundPionPair == 1) {
    EMiss -= m_TagPiPlusenergy + m_TagPiMinusenergy;
    PxMiss -= m_TagPiPluspx + m_TagPiMinuspx;
    PyMiss -= m_TagPiPluspy + m_TagPiMinuspy;
    PzMiss -= m_TagPiPluspz + m_TagPiMinuspz;
  }
  m_TagMissingEnergy = EMiss;
  m_TagMissingMass2 = EMiss*EMiss - PxMiss*PxMiss - PyMiss*PyMiss - PzMiss*PzMiss;
}
