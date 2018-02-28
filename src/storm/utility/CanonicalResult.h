#pragma once

#include <string>
#include <sstream>
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/adapters/RationalNumberAdapter.h"

namespace storm {
namespace utility {

class CanonicalResult {
public:
    static std::string getCanonicalString(double value);

#ifdef STORM_HAVE_CARL
    static std::string getCanonicalString(storm::RationalNumber value);
    static std::string getCanonicalString(storm::RationalFunction value);
#endif

    static std::string getHash(const std::string& s);

    template <typename T>
    static std::string getHashFromString(const T& v) {
        std::stringstream s;
        s << v;
        return getHash(s.str());
    }
};


}
}
