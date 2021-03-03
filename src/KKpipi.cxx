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
}

KKpipi::~KKpipi() {
}

StatusCode KKpipi::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Creating KKpipi Double Tag Algorithm" << endreq;
  StatusCode sc;
  if(m_recKpiTag) {
    sc = createSubAlgorithm("KKpipiVersusKpiDoubleTag", "KKpipiVersusKpiDoubleTag", m_KpiTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKpiDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKKTag) {
    sc = createSubAlgorithm("KKpipiVersusKKDoubleTag", "KKpipiVersusKKDoubleTag", m_KKTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKKDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recpipiTag) {
    sc = createSubAlgorithm("KKpipiVersuspipiDoubleTag", "KKpipiVersuspipiDoubleTag", m_pipiTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersuspipiDoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKpipi0Tag) {
    sc = createSubAlgorithm("KKpipiVersusKpipi0DoubleTag", "KKpipiVersusKpipi0DoubleTag", m_Kpipi0Tag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKpipi0DoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recpipipi0Tag) {
    sc = createSubAlgorithm("KKpipiVersuspipipi0DoubleTag", "KKpipiVersuspipipi0DoubleTag", m_pipipi0Tag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersuspipipi0DoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSpi0Tag) {
    sc = createSubAlgorithm("KKpipiVersusKSpi0DoubleTag", "KKpipiVersusKSpi0DoubleTag", m_KSpi0Tag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKSpi0DoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSpi0pi0Tag) {
    sc = createSubAlgorithm("KKpipiVersusKSpi0pi0DoubleTag", "KKpipiVersusKSpi0pi0DoubleTag", m_KSpi0pi0Tag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKSpi0pi0DoubleTag" << endreq;
      return StatusCode::FAILURE;
    }
  }
  if(m_recKSetaTag) {
    sc = createSubAlgorithm("KKpipiVersusKSetaDoubleTag", "KKpipiVersusKSetaDoubleTag", m_KSetaTag);
    if(sc.isFailure()) {
      log << MSG::ERROR << "Error while creating KKpipiVersusKSetaDoubleTag" << endreq;
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
