// Martin Duy Tat 12th February 2021
/**
 * KKpipi is a BOSS algorithm
 * It runs a double tag analysis on the signal mode \f$D\to K^+K^-\pi^+\pi^-\f$
 */

#ifndef KKPIPI
#define KKPIPI

// Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
// STL
#include<string>

class KKpipi: public Algorithm {
  public:
    /**
     * Default constructor for an algorithm where all necessary properties are declared
     */
    KKpipi(const std::string& name, ISvcLocator* pSvcLocator);
    /**
     * Trivial destructor
     */
    ~KKpipi();
    /**
     * This function runs when algorithm is initialized
     */
    StatusCode initialize();
    /**
     * Execution of the algorithm
     */
    StatusCode execute();
    /**
     * This function runs when algorithm is finalized
     */
    StatusCode finalize();
  private:
    /**
     * \f$K\pi\f$ tag mode
     */
    Algorithm *m_KpiTag;
    /**
     * \f$KK\f$ tag mode
     */
    Algorithm *m_KKTag;
    /**
     * \f$KK\f$ tag mode
     */
    Algorithm *m_pipiTag;
    /**
     * \f$K\pi\pi^0\f$ tag mode
     */
    Algorithm *m_Kpipi0Tag;
    /**
     * \f$\pi\pi\pi^0\f$ tag mode
     */
    Algorithm *m_pipipi0Tag;
    /**
     * \f$K_S\pi^0\f$ tag mode
     */
    Algorithm *m_KSpi0Tag;
    /**
     * Turn on \f$K\pi\f$ tag mode
     */
    bool m_recKpiTag;
    /**
     * Turn on \f$KK\f$ tag mode
     */
    bool m_recKKTag;
    /**
     * Turn on \f$KK\f$ tag mode
     */
    bool m_recpipiTag;
    /**
     * Turn on \f$K\pipi^0\f$ tag mode
     */
    bool m_recKpipi0Tag;
    /**
     * Turn on \f$\pi\pipi^0\f$ tag mode
     */
    bool m_recpipipi0Tag;
    /**
     * Turn on \f$K_S\pi^0\f$ tag mode
     */
    bool m_recKSpi0Tag;
};

#endif