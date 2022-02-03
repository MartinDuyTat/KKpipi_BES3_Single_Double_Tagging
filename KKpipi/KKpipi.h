// Martin Duy Tat 12th February 2021
/**
 * KKpipi is a BOSS algorithm
 * It runs a double tag analysis on the signal mode \f$D\to K^+K^-\pi^+\pi^-\f$
 * It runs a single tag analysis on the modes selected for double tag analysis
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
     * \f$KK\f$ tag mode
     */
    Algorithm *m_KKTag;
    /**
     * \f$KK\f$ tag mode, single tag
     */
    Algorithm *m_KKSingleTag;
    /**
     * \f$\pi\pi\f$ tag mode
     */
    Algorithm *m_pipiTag;
    /**
     * \f$\pi\pi\f$ tag mode, single tag
     */
    Algorithm *m_pipiSingleTag;
    /**
     * \f$K\pi\f$ tag mode
     */
    Algorithm *m_KpiTag;
    /**
     * \f$K\pi\f$ tag mode, single tag
     */
    Algorithm *m_KpiSingleTag;
    /**
     * \f$K\pi\pi^0\f$ tag mode
     */
    Algorithm *m_Kpipi0Tag;
    /**
     * \f$K\pi\pi^0\f$ tag mode, single tag
     */
    Algorithm *m_Kpipi0SingleTag;
    /**
     * \f$\pi\pi\pi^0\f$ tag mode
     */
    Algorithm *m_pipipi0Tag;
    /**
     * \f$\pi\pi\pi^0\f$ tag mode, single tag
     */
    Algorithm *m_pipipi0SingleTag;
    /**
     * \f$K_S\pi^0\f$ tag mode
     */
    Algorithm *m_KSpi0Tag;
    /**
     * \f$K_S\pi^0\f$ tag mode, single tag
     */
    Algorithm *m_KSpi0SingleTag;
    /**
     * \f$K_S\pi^0\pi^0\f$ tag mode
     */
    Algorithm *m_KSpi0pi0Tag;
    /**
     * \f$K_S\pi^0\pi^0\f$ tag mode, single tag
     */
    Algorithm *m_KSpi0pi0SingleTag;
    /**
     * \f$K_S\eta\f$ tag mode
     */
    Algorithm *m_KSetaTag;
    /**
     * \f$K_S\eta\f$ tag mode, single tag
     */
    Algorithm *m_KSetaSingleTag;
    /**
     * \f$K_S\eta'(\pi\pi\eta)\f$ tag mode
     */
    Algorithm *m_KSetaPrimepipietaTag;
    /**
     * \f$K_S\eta'(\pi\pi\eta)\f$ tag mode, single tag
     */
    Algorithm *m_KSetaPrimepipietaSingleTag;
    /**
     * \f$K_S\eta'(\rho\gamma)\f$ tag mode
     */
    Algorithm *m_KSetaPrimerhogammaTag;
    /**
     * \f$K_S\eta'(\rho\gamma)\f$ tag mode, single tag
     */
    Algorithm *m_KSetaPrimerhogammaSingleTag;
    /**
     * \f$K_S\pi\pi\pi^0\f$ tag mode
     */
    Algorithm *m_KSpipipi0Tag;
    /**
     * \f$K_S\pi\pi\pi^0\f$ tag mode, single tag
     */
    Algorithm *m_KSpipipi0SingleTag;
    /**
     * \f$K_SKK\f$ tag mode
     */
    Algorithm *m_KSKKTag;
    /**
     * \f$K_SKK\f$ tag mode, single tag
     */
    Algorithm *m_KSKKSingleTag;
    /**
     * \f$K_S\pi\pi\f$ tag mode
     */
    Algorithm *m_KSpipiTag;
    /**
     * \f$K_S\pi\pi\f$ tag mode, partially reconstructed
     */
    Algorithm *m_KSpipiPartRecoTag;
    /**
     * \f$K_S\pi\pi\f$ tag mode, single tag
     */
    Algorithm *m_KSpipiSingleTag;
    /**
     * \f$KK\pi\pi\f$ tag mode
     */
    Algorithm *m_KKpipiTag;
    /**
     * \f$KK\pi\pi\f$ tag mode, single tag
     */
    Algorithm *m_KKpipiSingleTag;
    /**
     * \f$K\pi\pi\pi\f$ tag mode
     */
    Algorithm *m_KpipipiTag;
    /**
     * \f$K\pi\pi\pi\f$ tag mode, single tag
     */
    Algorithm *m_KpipipiSingleTag;
    /**
     * \f$K_L\pi^0\f$ tag mode
     */
    Algorithm *m_KLpi0Tag;
    /**
     * \f$K_L\pi^0\pi^0\f$ tag mode
     */
    Algorithm *m_KLpi0pi0Tag;
    /**
     * \f$K_L\pi\pi\pi^0\f$ tag mode
     */
    Algorithm *m_KLpipipi0Tag;
    /**
     * \f$K_L\pi\pi\f$ tag mode
     */
    Algorithm *m_KLpipiTag;
    /**
     * \f$K_LKK\f$ tag mode
     */
    Algorithm *m_KLKKTag;
    /**
     * \f$Ke\nu\f$ tag mode
     */
    Algorithm *m_KeNuTag;
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
    /**
     * Turn on \f$K_S\pi^0\pi^0\f$ tag mode
     */
    bool m_recKSpi0pi0Tag;
    /**
     * Turn on \f$K_S\eta\f$ tag mode
     */
    bool m_recKSetaTag;
    /**
     * Turn on \f$K_S\pi\pi\pi^0\f$ tag mode
     */
    bool m_recKSpipipi0Tag;
    /**
     * Turn on \f$K_SKK\f$ tag mode
     */
    bool m_recKSKKTag;
    /**
     * Turn on \f$K_S\eta'(\pi\pi\eta)\f$ tag mode
     */
    bool m_recKSetaPrimepipietaTag;
    /**
     * Turn on \f$K_S\eta'(\rho\gamma)\f$ tag mode
     */
    bool m_recKSetaPrimerhogammaTag;
    /**
     * Turn on \f$KK\pi\pi\f$ tag mode
     */
    bool m_recKKpipiTag;
    /**
     * Turn on \f$K_S\pi\pi\f$ tag mode
     */
    bool m_recKSpipiTag;
    /**
     * Turn on \f$K\pi\pi\pi\f$ tag mode
     */
    bool m_recKpipipiTag;
    /**
     * Turn on \f$K_L\pi^0\f$ tag mode
     */
    bool m_recKLpi0Tag;
    /**
     * Turn on \f$K_L\pi^0\pi^0\f$ tag mode
     */
    bool m_recKLpi0pi0Tag;
    /**
     * Turn on \f$K_L\pi\pi\pi^0\f$ tag mode
     */
    bool m_recKLpipipi0Tag;
    /**
     * Turn on \f$K_L\pi\pi\f$ tag mode
     */
    bool m_recKLpipiTag;
    /**
     * Turn on \f$K_LKK\f$ tag mode
     */
    bool m_recKLKKTag;
    /**
     * Turn on \f$Ke\nu\f$ tag mode
     */
    bool m_recKeNuTag;
};

#endif
