#include "GaudiKernel/DeclareFactoryEntries.h"
#include "KKpipi/KKpipi.h"
#include "KKpipi/KKpipiSingleTag.h"
#include "KKpipi/KKpipiVersusKpiDoubleTag.h"
#include "KKpipi/KKpipiVersusKpipi0DoubleTag.h"
#include "KKpipi/KKpipiVersusKKDoubleTag.h"
#include "KKpipi/KKpipiVersuspipiDoubleTag.h"
#include "KKpipi/KKpipiVersuspipipi0DoubleTag.h"
#include "KKpipi/KKpipiVersusKSpi0DoubleTag.h"
#include "KKpipi/KKpipiVersusKSpi0pi0DoubleTag.h"
#include "KKpipi/KKpipiVersusKSetaDoubleTag.h"
#include "KKpipi/KKpipiVersusKSpipipi0DoubleTag.h"
#include "KKpipi/KKpipiVersusKSKKDoubleTag.h"
#include "KKpipi/KKpipiVersusKSetaPrimepipietaDoubleTag.h"
#include "KKpipi/KKpipiVersusKSetaPrimerhogammaDoubleTag.h"
#include "KKpipi/KKpipiVersusKKpipiDoubleTag.h"
#include "KKpipi/KKpipiVersusKSpipiDoubleTag.h"
#include "KKpipi/KKSingleTag.h"
#include "KKpipi/pipiSingleTag.h"
#include "KKpipi/KpiSingleTag.h"
#include "KKpipi/Kpipi0SingleTag.h"
#include "KKpipi/pipipi0SingleTag.h"
#include "KKpipi/KSpi0SingleTag.h"
#include "KKpipi/KSpi0pi0SingleTag.h"
#include "KKpipi/KSetaSingleTag.h"
#include "KKpipi/KSetaPrimepipietaSingleTag.h"
#include "KKpipi/KSetaPrimerhogammaSingleTag.h"
#include "KKpipi/KSpipipi0SingleTag.h"
#include "KKpipi/KSKKSingleTag.h"
#include "KKpipi/KSpipiSingleTag.h"

DECLARE_ALGORITHM_FACTORY(KKpipi)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKKDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersuspipiDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKpiDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKpipi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersuspipipi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSpi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSpi0pi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSetaDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSetaPrimepipietaDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSetaPrimerhogammaDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSpipipi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSKKDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSpipiDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKKpipiDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKSingleTag)
DECLARE_ALGORITHM_FACTORY(pipiSingleTag)
DECLARE_ALGORITHM_FACTORY(KpiSingleTag)
DECLARE_ALGORITHM_FACTORY(Kpipi0SingleTag)
DECLARE_ALGORITHM_FACTORY(pipipi0SingleTag)
DECLARE_ALGORITHM_FACTORY(KSpi0SingleTag)
DECLARE_ALGORITHM_FACTORY(KSpi0pi0SingleTag)
DECLARE_ALGORITHM_FACTORY(KSetaSingleTag)
DECLARE_ALGORITHM_FACTORY(KSetaPrimepipietaSingleTag)
DECLARE_ALGORITHM_FACTORY(KSetaPrimerhogammaSingleTag)
DECLARE_ALGORITHM_FACTORY(KSpipipi0SingleTag)
DECLARE_ALGORITHM_FACTORY(KSKKSingleTag)
DECLARE_ALGORITHM_FACTORY(KSpipiSingleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiSingleTag)

DECLARE_FACTORY_ENTRIES(KKpipi) {
  DECLARE_ALGORITHM(KKpipi)
}


