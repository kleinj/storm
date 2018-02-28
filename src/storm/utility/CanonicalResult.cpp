#include "storm/utility/CanonicalResult.h"
#include "boost/uuid/sha1.hpp"

namespace storm {
namespace utility {

std::string CanonicalResult::getCanonicalString(double value) {
    return std::to_string(value);
}

#ifdef STORM_HAVE_CARL

std::string CanonicalResult::getCanonicalString(storm::RationalNumber value) {
    std::stringstream out;
    out << value;
    return out.str();
}

std::string CanonicalResult::getCanonicalString(storm::RationalFunction value) {
    typedef carl::RationalFunction<RawPolynomial, true> NonFactorizedRationalFunction;

    auto num = carl::computePolynomial(value.nominator());
    auto den = carl::computePolynomial(value.denominator());

    NonFactorizedRationalFunction expanded(num, den);

    num = expanded.nominator();
    den = expanded.denominator();
    num.makeOrdered();
    den.makeOrdered();

    std::stringstream out;
    out << num << "/" << den;
    return out.str();
}

std::string CanonicalResult::getHash(const std::string& value) {
    boost::uuids::detail::sha1 sha1;
    sha1.process_bytes(value.c_str(), value.length());

    unsigned int digest[5];
    sha1.get_digest(digest);

    std::stringstream out;
    out << std::hex; // hex output
    for (unsigned int i=0;i<5;i++) {
        out << digest[i];
    }
    return out.str();
}

#endif

}
}
