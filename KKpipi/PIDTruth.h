// Martin Duy Tat 23rd March 2021
/**
 * PIDTruth is a class used to check the true PID values of reconstructed tracks
 */

// Gaudi
#include "GaudiKernel/Algorithm.h"
// STL
#include <vector>

#ifndef PIDTRUTH
#define PIDTRUTH

class PIDTruth {
  public:
    /**
     * Constructor that takes in a vector of track IDs of reconstructed tracks and an Algorithm object and saves these
     * @param TrackID Vector of track IDs
     * @param algorithm Pointer to the Algorithm object that created this class, need this to access the MC info!
     */
  PIDTruth(const std::vector<int> &TrackID, const Algorithm *algorithm);
    /**
     * Function from Alex Gilman that matches reconstructed particles with generator particles by comparing the hits of the reconstructed object and the truth level trajectory
     * @param tkID Track ID of particle of interest
     * @param mcPDG Set to -1 to find PID of the particle of interest, set to the reconstructed PID to check if the PID has been correctly assigned, set to 0 if we're not interested in this
     * @param mcParPDG Set to -1 to find PID of the parent, set to the reconstructed PID to check if the PID has been correctly assigned to the parent, set to 0 if we're not interested in this
     * @param GParPDG Set to the PID of any parent or grand parent if we want to check if the particle of interest originated from this particle through a (potentially long) decay chain, set to 0 if we're not interested in this
     * @return Returns 0 for false, 1 for true or the PDG PID number where appropriate
     */
    int MCTKPIDCHG(int tkID, int mcPDG, int mcParPDG, int GParPDG) const;
    /**
     * Function that checks if all the particles originate from the same \f$D\f$ meson
     */
    bool SameDMother() const;
  private:
    /**
     * Vector with track IDs of reconstructed daughters
     */
    std::vector<int> m_TrackID;
    /**
     * Pointer to the Algorithm object that created this class, need this to access the MC info!
     */
    Algorithm *m_algorithm;
}

#endif
