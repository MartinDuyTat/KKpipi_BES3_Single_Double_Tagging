// Martin Duy Tat 12th February 2021

// KKpipi
#include "KKpipi/KKpipi.h"
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/StatusCode.h"
// STL
#include<string>

KKpipi::KKpipi(const std::string& name, ISvcLocator* pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("KKpipiVersusKpiDoubleTag", m_recKpiTag = true);
  declareProperty("KKpipiVersusKKDoubleTag", m_recKKTag = true);
  declareProperty("KKpipiVersuspipiDoubleTag", m_recpipiTag = true);
  declareProperty("KKpipiVersusKpipi0DoubleTag", m_recKpipi0Tag = true);
  declareProperty("KKpipiVersuspipipi0DoubleTag", m_recpipipi0Tag = true);
  declareProperty("KKpipiVersusKSpi0DoubleTag", m_recKSpi0Tag = true);
  declareProperty("KKpipiVersusKSpi0pi0DoubleTag", m_recKSpi0pi0Tag = true);
  declareProperty("KKpipiVersusKSetaDoubleTag", m_recKSetaTag = true);
  declareProperty("KKpipiVersusKSpipipi0DoubleTag", m_recKSpipipi0Tag = true);
  declareProperty("KKpipiVersusKSKKDoubleTag", m_recKSKKTag = true);
  declareProperty("KKpipiVersusKSetaPrimepipietaDoubleTag", m_recKSetaPrimepipietaTag = true);
  declareProperty("KKpipiVersusKSetaPrimerhogammaDoubleTag", m_recKSetaPrimerhogammaTag = true);
  declareProperty("KKpipiVersusKKpipiDoubleTag", m_recKKpipiTag = true);
  declareProperty("KKpipiVersusKSpipiDoubleTag", m_recKSpipiTag = true);
}

KKpipi::~KKpipi() {
}

StatusCode KKpipi::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Creating KKpipi Double Tag Algorithm" << endreq;
  StatusCode sc;
  if(m_recKKTag) {
    sc = createSubAlgorithm("KKpipiVersusKKDoubleTag", "KKpipiVersusKKDoubleTag", m_KKTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKKDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("KKSingleTag", "KKSingleTag", m_KKSingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKSingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recpipiTag) {
    sc = createSubAlgorithm("KKpipiVersuspipiDoubleTag", "KKpipiVersuspipiDoubleTag", m_pipiTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersuspipiDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("pipiSingleTag", "pipiSingleTag", m_pipiSingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating pipiSingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKpiTag) {
    sc = createSubAlgorithm("KKpipiVersusKpiDoubleTag", "KKpipiVersusKpiDoubleTag", m_KpiTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKpiDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("KpiSingleTag", "KpiSingleTag", m_KpiSingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KpiSingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKpipi0Tag) {
    sc = createSubAlgorithm("KKpipiVersusKpipi0DoubleTag", "KKpipiVersusKpipi0DoubleTag", m_Kpipi0Tag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKpipi0DoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("Kpipi0SingleTag", "Kpipi0SingleTag", m_Kpipi0SingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating Kpipi0SingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recpipipi0Tag) {
    sc = createSubAlgorithm("KKpipiVersuspipipi0DoubleTag", "KKpipiVersuspipipi0DoubleTag", m_pipipi0Tag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersuspipipi0DoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("pipipi0SingleTag", "pipipi0SingleTag", m_pipipi0SingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating pipipi0SingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSpi0Tag) {
    sc = createSubAlgorithm("KKpipiVersusKSpi0DoubleTag", "KKpipiVersusKSpi0DoubleTag", m_KSpi0Tag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKSpi0DoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("KSpi0SingleTag", "KSpi0SingleTag", m_KSpi0SingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSpi0SingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSpi0pi0Tag) {
    sc = createSubAlgorithm("KKpipiVersusKSpi0pi0DoubleTag", "KKpipiVersusKSpi0pi0DoubleTag", m_KSpi0pi0Tag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKSpi0pi0DoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("KSpi0pi0SingleTag", "KSpi0pi0SingleTag", m_KSpi0pi0SingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSpi0pi0SingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSetaTag) {
    sc = createSubAlgorithm("KKpipiVersusKSetaDoubleTag", "KKpipiVersusKSetaDoubleTag", m_KSetaTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKSetaDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("KSetaSingleTag", "KSetaSingleTag", m_KSetaSingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSetaSingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSetaPrimepipietaTag) {
    sc = createSubAlgorithm("KKpipiVersusKSetaPrimepipietaDoubleTag", "KKpipiVersusKSetaPrimepipietaDoubleTag", m_KSetaPrimepipietaTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKSetaPrimepipietaDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("KSetaPrimepipietaSingleTag", "KSetaPrimepipietaSingleTag", m_KSetaPrimepipietaSingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSetaPrimepipietaSingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSetaPrimerhogammaTag) {
    sc = createSubAlgorithm("KKpipiVersusKSetaPrimerhogammaDoubleTag", "KKpipiVersusKSetaPrimerhogammaDoubleTag", m_KSetaPrimerhogammaTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKSetaPrimerhogammaDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("KSetaPrimerhogammaSingleTag", "KSetaPrimerhogammaSingleTag", m_KSetaPrimerhogammaSingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSetaPrimerhogammaSingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSpipipi0Tag) {
    sc = createSubAlgorithm("KKpipiVersusKSpipipi0DoubleTag", "KKpipiVersusKSpipipi0DoubleTag", m_KSpipipi0Tag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKSpipipi0DoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("KSpipipi0SingleTag", "KSpipipi0SingleTag", m_KSpipipi0SingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSpipipi0SingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSKKTag) {
    sc = createSubAlgorithm("KKpipiVersusKSKKDoubleTag", "KKpipiVersusKSKKDoubleTag", m_KSKKTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKSKKDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("KSKKSingleTag", "KSKKSingleTag", m_KSKKSingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSKKSingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSpipiTag) {
    sc = createSubAlgorithm("KKpipiVersusKSpipiDoubleTag", "KKpipiVersusKSpipiDoubleTag", m_KSpipiTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKSpipiDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("KSpipiDoubleTag", "KSpipiDoubleTag", m_KSpipiSingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KSpipiSingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKKpipiTag) {
    sc = createSubAlgorithm("KKpipiVersusKKpipiDoubleTag", "KKpipiVersusKKpipiDoubleTag", m_KKpipiTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKKpipiDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
    sc = createSubAlgorithm("KKpipiSingleTag", "KKpipiSingleTag", m_KKpipiSingleTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiSingleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
}

StatusCode KKpipi::execute() {
  for(std::vector<Algorithm*>::const_iterator it = subAlgorithms()->begin(); it != subAlgorithms()->end(); it++) {
    StatusCode sc = (*it)->execute();
    if(sc.isFailure()) {
      return StatusCode::FAILURE;
    }
  }
}

StatusCode KKpipi::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "KKpipi Double Tag Algorithm finalized" << endmsg;
  return StatusCode::SUCCESS;
}
