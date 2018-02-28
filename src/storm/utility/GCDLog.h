#pragma once

#include "storm/adapters/RationalFunctionAdapter.h"

#ifdef STORM_HAVE_CARL

#include <carl/core/RationalFunction.h>

#ifdef CARL_RATIONAL_FUNCTION_HAVE_GCD_LOG
#define STORM_HAVE_GCD_LOG 1

namespace storm {
namespace utility {

    class GCDLog {
    public:
        static void enable();
        static void disable();
    };

}
}

#endif
#endif

