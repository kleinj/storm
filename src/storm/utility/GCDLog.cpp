#include "storm/utility/GCDLog.h"

#ifdef STORM_HAVE_GCD_LOG

#include <carl/core/GcdLog.h>
#include "storm/adapters/RationalFunctionAdapter.h"

namespace storm {
namespace utility {

    static const bool debug = false;
    static std::chrono::nanoseconds gcd_time = std::chrono::nanoseconds::zero();
    static std::chrono::nanoseconds gcd_time_min = std::chrono::nanoseconds::zero();
    static std::chrono::nanoseconds gcd_time_max = std::chrono::nanoseconds::zero();
    static std::size_t gcd_calls = 0;

    static void hook(void *rat_fun, bool first, std::chrono::nanoseconds duration) {
        const storm::RationalFunction* f = static_cast<const storm::RationalFunction*>(rat_fun);

        if (debug) {
            std::cout << (first? "Before GCD: " : "After GCD: " ) << *f << std::endl;
            if (!first) {
                std::cout << duration.count() << "ns" << std::endl;
            }
        }
	if (!first) {
	  //	  std::cout << "GCD took: " << duration.count() << "ns<< std::endl;
	  gcd_time += duration;

	  if (gcd_calls == 0) {
	    gcd_time_min = duration;
	    gcd_time_max = duration;
	  } else {
	    if (gcd_time_min > duration)
	      gcd_time_min = duration;
	    if (gcd_time_max < duration)
	      gcd_time_max = duration;
	  }
	  
	  gcd_calls++;
	}
    }

    void GCDLog::enable() {
        gcd_time = std::chrono::nanoseconds::zero();
	gcd_calls = 0;
        carl::GcdLog::register_hook(hook);
    }

    void GCDLog::disable() {
        carl::GcdLog::reset();
	std::cout << "GCD calls: " << gcd_calls
		  << ", GCD time: " << std::chrono::duration<double, std::milli>(gcd_time).count() << "ms (combined), "
		  << std::chrono::duration<double, std::milli>(gcd_time_min).count() << "ms (min), "
		  << std::chrono::duration<double, std::milli>(gcd_time_max).count() << "ms (max), "
		  << std::endl;
    }

}
}

#endif
