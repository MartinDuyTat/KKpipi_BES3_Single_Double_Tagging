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

DECLARE_ALGORITHM_FACTORY(KKpipi)
DECLARE_ALGORITHM_FACTORY(KKpipiSingleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKpiDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKpipi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKKDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersuspipiDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersuspipipi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSpi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSpi0pi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSetaDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSpipipi0DoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSKKDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSetaPrimepipietaDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSetaPrimerhogammaDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKKpipiDoubleTag)
DECLARE_ALGORITHM_FACTORY(KKpipiVersusKSpipiDoubleTag)

DECLARE_FACTORY_ENTRIES(KKpipi) {
  DECLARE_ALGORITHM(KKpipi)
}


