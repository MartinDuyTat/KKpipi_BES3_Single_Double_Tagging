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
// STL
#include <vector>
#include <utility>
#include <algorithm>

PIDTruth::PIDTruth(const std::vector<int> &TrackID, int NumberCharged, const Algorithm *algorithm, const std::vector<std::pair<int, int>> &PhotonPairTrackID): m_TrackID(TrackID), m_NumberCharged(NumberCharged), m_algorithm(algorithm), m_PhotonPairTrackID(PhotonPairTrackID) {
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

int PIDTruth::MCSHPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG) const {
  SmartDataPtr < EvtRecEvent > evtRecEvent(m_algorithm->eventSvc(), "/Event/EvtRec/EvtRecEvent");
  SmartDataPtr < EvtRecTrackCol > evtRecTrkCol(m_algorithm->eventSvc(), "/Event/EvtRec/EvtRecTrackCol");
  SmartDataPtr < EventNavigator > navigator(m_algorithm->eventSvc(), "/Event/Navigator");
  int ismatched = -1;
  for (int i = evtRecEvent -> totalCharged(); i < evtRecEvent -> totalTracks(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol -> begin() + i;
    if (!( * itTrk) -> isEmcShowerValid()) {
      continue;
    }
    RecEmcShower * emcTrk = ( * itTrk) -> emcShower();
    if (( * itTrk) -> trackId() == tkID) {
      McParticleVector particles = navigator -> getMcParticles(emcTrk);
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
        int relevance = navigator -> getMcParticleRelevance(emcTrk, particles[i]);
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

      if (mcParPDG == -1) {
        return thebestParent;
      }
      if (mcPDG == -1) {
        return thebestmatchedPDG;
      }

      if (mcPDG == 0 && mcParPDG == 0) {
        if (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
          thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG) {
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
          (thebestParent == GParPDG || thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG ||
            thebestGGGGParent == GParPDG || thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          return 1;
        }
      } else {
        if (thebestmatchedPDG == mcPDG && thebestParent == mcParPDG &&
          (thebestGParent == GParPDG || thebestGGParent == GParPDG || thebestGGGParent == GParPDG || thebestGGGGParent == GParPDG ||
            thebestGGGGGParent == GParPDG || thebestGGGGGGParent == GParPDG || thebestGGGGGGGParent == GParPDG)) {
          ismatched = 1;
        }
      }
    }
  }
  return ismatched;
}

int PIDTruth::FindDOrigin(int TrackID, bool Charged) const {
  if(Charged) {
    if(MCTKPIDCHG(TrackID, 0, 0, 421) == 1) {
      return 421;
    } else if(MCTKPIDCHG(TrackID, 0, 0, -421) == 1) {
      return -421;
    } else {
      return 0;
    }
  } else {
    if(MCSHPIDCHG(TrackID, 0, 0, 421) == 1) {
      return 421;
    } else if(MCSHPIDCHG(TrackID, 0, 0, -421) == 1) {
      return -421;
    } else {
      return 0;
    }
  }
}

bool PIDTruth::SameDMother() const {
  // Check D meson ID of first track ID
  int DMotherID = FindDOrigin(m_TrackID[0], m_NumberCharged != 0);
  // Loop over all other track IDs and compare their D meson ID to the first one
  for(std::vector<int>::const_iterator iter = m_TrackID.begin() + 1; iter != m_TrackID.end(); iter++) {
    int ThisDMotherID = iter - m_TrackID.begin() < m_NumberCharged ? FindDOrigin(*iter, true) : FindDOrigin(*iter, false);
    // If their D meson IDs are not equal, these tracks don't originate from the same D mother
    if(DMotherID != ThisDMotherID) {
      return false;
    }
  }
  // Repeat, but for pairs of photons only require one of them to pass the truth matching
  for(std::vector<std::pair<int, int>>::const_iterator iter = m_PhotonPairTrackID.begin(); iter != m_PhotonPairTrackID.end(); iter++) {
    int ThisMotherIDFirst = FindDOrigin(iter->first, false);
    int ThisMotherIDSecond = FindDOrigin(iter->second, false);
    if(DMotherID != ThisMotherIDFirst && DMotherID != ThisMotherIDSecond) {
      return false;
    }
  }
  // If we make it this far, all tracks and photons originate from the same D meson
  return true;
}

bool PIDTruth::FindTrueID(std::vector<int> &ParticleID) const {
  std::vector<int> TrueID;
  bool PIDMatch = true;
  // Loop over all track IDs
  for(unsigned int i = 0; i < m_TrackID.size(); i++) {
    if(i < m_NumberCharged) {
      // Charged tracks
      TrueID.push_back(MCTKPIDCHG(m_TrackID[i], -1, 0, 0));
    } else {
      // Neutral tracks
      TrueID.push_back(MCSHPIDCHG(m_TrackID[i], -1, 0, 0));
    }
    // Check if PID matches
    if(PIDMatch && ParticleID[i] != TrueID[i] && ParticleID[i] != 0) {
      PIDMatch = false;
    }
  }
  // Repeat for the photon pairs to obtain their true PID
  for(unsigned int i = 0; i < m_PhotonPairTrackID.size(); i++) {
    TrueID.push_back(MCSHPIDCHG(m_PhotonPairTrackID[i].first, -1, 0, 0));
    TrueID.push_back(MCSHPIDCHG(m_PhotonPairTrackID[i].second, -1, 0, 0));
    // No PID truth matching for photon pairs
  }
  std::swap(ParticleID, TrueID);
  return PIDMatch;
}

int PIDTruth::GetTrueMotherID(int TrackID, bool Charged) const {
  if(Charged) {
    return MCTKPIDCHG(TrackID, 0, -1, 0);
  } else {
    return MCSHPIDCHG(TrackID, 0, -1, 0);
  }
}
