#include "GaudiKernel/DeclareFactoryEntries.h"
#include "KKpipi/KKpipiSingleTag.h"
#include "KKpipi/KKpipiVersusKpiDoubleTag.h"
#include "KKpipi/KKpipiVersusKpipi0DoubleTag.h"
#include "KKpipi/KKpipiVersusKKDoubleTag.h"
#include "KKpipi/KKpipiVersusKSpi0DoubleTag.h"

DECLARE_ALGORITHM_FACTORY(KKpipiSingleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKpiDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKpipi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKKDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSpi0DoubleTag)

DECLARE_FACTORY_ENTRIES(KKpipi) {
  DECLARE_ALGORITHM(KKpipi)
}


