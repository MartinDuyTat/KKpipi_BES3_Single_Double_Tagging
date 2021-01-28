// Martin Duy Tat 28th January 2021

// Header file
#include "KKpipiSingleTag.h"
// Gaudi
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
// Event information
#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecDTag.h"
// Boss
#include "DTagTool/DTagTool.h"
#include "SimplePIDSvc/ISimplePIDSvc.h"
// STL
#include<vector>
#include<string>

KKpipiSingleTag::KKpipiSingleTag(const std::string &name, ISvcLocator *pSvcLocator): Algorithm(name, pSvcLocator) {
  declareProperty("dummy", m_dummy = 0);
}

KKpipiSingleTag::~KKpipiSingleTag() {
}

KKpipiSingleTag::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "Initializing KKpipi Single Tagging" << endreq;
  Statuscode status;
  NTuplePtr ntp(ntupleSvc(), "KKPIPI/SingleTag");
  if(ntp) {
    m_tuple = ntp;
  } else {
    m_tuple = ntupleSvc()->book("KKPIPI/SingleTag", CLID_ColumnWiseTuple, "Single tagged D->KKpipi events");
    if(m_tuple) {
      status = m_tuple->addItem("Run", m_RunNumber);
      status = m_tuple->addItem("Event", m_EventNumber);
      status = m_tuple->addItem("NumberOfParticles", m_NumberParticles);
      status = m_tuple->addIndexedItem("ParticleIDs", m_pdgID);
      status = m_tuple->addIndexedItem("MotherIndex", m_MotherIndex);
      status = m_tuple->addItem("GeneratorNumberOfParticles", m_GeneratorNumberParticles);
      status = m_tuple->addItem("GeneratorParticleIDs", m_GeneratorPDGID);
      status = m_tuple->addItem("GeneratorMotherID", m_MotherID);
      status = m_tuple->addItem("True_P", m_TrueMomentum);
      status = m_tuple->addItem("True_PT", m_TruePT);
      status = m_tuple->addItem("True_phi", m_TruePhi);
      status = m_tuple->addItem("True_theta", m_Theta);
      status = m_tuple->addItem("Charm", m_Charm);
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
      status = m_tuple->addItem("KPluspx", m_KPluspx);
      status = m_tuple->addItem("KPluspy", m_KPluspy);
      status = m_tuple->addItem("KPluspz", m_KPluspz);
      status = m_tuple->addItem("KPlusenergy", m_KPlusenergy);
      status = m_tuple->addItem("KMinuspx", m_KMinuspx);
      status = m_tuple->addItem("KMinuspy", m_KMinuspy);
      status = m_tuple->addItem("KMinuspz", m_KMinuspz);
      status = m_tuple->addItem("KMinusenergy", m_KMinusenergy);
    } else {
      log << MSG::ERROR << "Cannot book NTuple for KKpipi Single Tags" << endmsg;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;
  }
}
