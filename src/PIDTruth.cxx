// Martin Duy Tat 23rd March 2021

// KKpipi
#include "KKpipi/PIDTruth.h"
// Gaudi
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/SmartDataPtr.h"
// BOSS
#include "EventNavigator/EventNavigator.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "McTruth/McParticle.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include 
// STL
#include <vector>
#include <algorithm>

PIDTruth::PIDTruth(const std::vector<int> &TrackID, const Algorithm *algorithm): m_TrackID(Track_ID), m_algorithm(algorithm) {
}

int PIDTruth::MCTKPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG) const {
  SmartDataPtr < EvtRecEvent > evtRecEvent(m_algorithm->eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(m_algorithm->eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(m_algorithm->eventSvc(), "/Event/Navigator");
  int ismatched = 0;
  for (int i = 0; i < evtRecEvent -> totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (( * itTrk) -> trackId() == tkID) {
      if (!( * itTrk) -> isMdcKalTrackValid()) {
        continue;
      }
      RecMdcKalTrack * mdcTrk = ( * itTrk) -> mdcKalTrack();
      McParticleVector particles = navigator -> getMcParticles(mdcTrk);
      int temp_hits = -2;
      int thebestmatchedPDG = -2;
      int thebestParent = -2;
      int thebestGParent = -2;
      int thebestGGParent = -2;
      int thebestGGGParent = -2;
      int thebestGGGGParent = -2;
      int thebestGGGGGParent = -2;
      int thebestGGGGGGParent = -2;
      int thebestGGGGGGGParent = -2;
      for (unsigned int i = 0; i < particles.size(); i++) {
        int relevance = navigator -> getMcParticleRelevance(mdcTrk, particles[i]);
        if (relevance > temp_hits) {
          temp_hits = relevance;
          thebestmatchedPDG = particles[i] -> particleProperty();
          thebestParent = particles[i] -> mother().particleProperty();
          thebestGParent = particles[i] -> mother().mother().particleProperty();
          thebestGGParent = particles[i] -> mother().mother().mother().particleProperty();
          thebestGGGParent = particles[i] -> mother().mother().mother().mother().particleProperty();
          thebestGGGGParent = particles[i] -> mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().particleProperty();
          thebestGGGGGGGParent = particles[i] -> mother().mother().mother().mother().mother().mother().mother().mother().particleProperty();
        }
      }
      if (mcParPDG == -1 && (GParPDG == 0 || GParPDG == thebestGParent)) {
        return thebestParent;
      }
      if (mcParPDG == -1 && GParPDG == 0) {
        return thebestParent;
      }
      if (mcPDG == -1) {
        return thebestmatchedPDG;
      }
      if (mcPDG == 0 && mcParPDG == 0) {
        if (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG) {
          return 1;
        }
      }
      if (GParPDG == 0 && mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG) {
          return 1;
        }
      }
      if (mcPDG == 0 && GParPDG == 0) {
        if (thebestParent == mcParPDG) {
          return 1;
        }
      }
      if (mcParPDG == 0) {
        if (thebestmatchedPDG == mcPDG &&
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          return 1;
        }
      } else {
        if (thebestmatchedPDG == mcPDG && thebestParent == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          ismatched = 1;
        }
      }
    }
  }
  return ismatched;
}

bool PIDTruth::SameDMother() const {
  std::vector<int> IsD0Mother, IsD0barMother;
  for(auto TrackID : m_TrackID) {
    IsD0Mother.push_back(MCTKPIDCHG(TrackID, 0, 0, 421));
    IsD0barMother.push_back(MCTKPIDCHG(TrackID, 0, 0, -421));
  }
  if(std::all_of(IsD0Mother.begin(), IsD0Mother.end(), [](int i) {return i == 1})) {
    return true;
  } else if (std::all_of(IsD0barMother.begin(), IsD0barMother.end(), [](int i) {return i == 1})) {
    return true;
  } else {
    return false;
  }
}
