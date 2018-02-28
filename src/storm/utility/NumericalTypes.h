#pragma once

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"

namespace storm {
namespace utility {

namespace NumericalTypes {

    template <typename ValueType>
    inline std::string getTypeName();

    template <typename ValueType>
    inline std::size_t getStringComplexity(const ValueType& v) {
        std::stringstream ss;
        ss << v;
        return ss.str().length();
    }

    template <typename ValueType>
    inline std::string getStats(const ValueType& v) {
        return std::to_string(getStringComplexity(v));
    }

    template <typename ValueType>
    inline std::size_t getComplexity(const ValueType& v) {
        return getStringComplexity(v);
    }

    template <typename ValueType>
    inline std::string getStats(const std::vector<ValueType>& values) {
        std::stringstream ss;
        for (std::size_t i = 0; i < values.size(); i++) {
            ss << i << ": " << getStats(values[i]) << "\n";
        }

        return ss.str();
    }

    // specializations

    template <>
    inline std::string getTypeName<double>() {
        return "double";
    }

#ifdef STORM_HAVE_CARL

    template <>
    inline std::string getTypeName<storm::RationalFunctionCoefficient>() {
        return "RationalFunctionCoefficient";
    }

    template <>
    inline std::string getTypeName<storm::RationalFunction>() {
        return "RationalFunction(factorized,simplify)";
    }

    template <>
    inline std::string getTypeName<PlainRationalFunction>() {
        return "RationalFunction(plain,simplify)";
    }

    template <>
    inline std::string getTypeName<PlainRationalFunctionNoSimplify>() {
        return "RationalFunction(plain,nosimplify)";
    }

    template <>
    inline std::string getTypeName<FactorizedRationalFunctionNoSimplify>() {
        return "RationalFunction(factorized,nosimplify)";
    }

    template <>
    inline std::string getTypeName<storm::RawPolynomial>() {
        return "MultivariatePolynomial(with rational numbers)";
    }

    template <>
    inline std::string getTypeName<storm::Polynomial>() {
        return "Product of MultivariatePolynomial(with rational numbers)";
    }

    template <typename P, bool AS>
    inline std::string getStats(const typename carl::RationalFunction<P,AS>& f) {
        std::stringstream ss;
        ss << "N[" << getStats(f.nominator());
        ss << "]/D[" << getStats(f.denominator()) << "]";

        if (!AS) {
            carl::RationalFunction<P,AS> simplified(f);
            simplified.simplify();
            ss << " sN[" << getStats(simplified.nominator());
            ss << "]/sD[" << getStats(simplified.denominator()) << "]";
        }

        return ss.str();
    }

    template <>
    inline std::string getStats(const storm::RawPolynomial& p) {
        std::stringstream ss;

        ss << "L=" << getStringComplexity(p);
        ss << ",M=" << (p.isConstant() ? 1 : p.nrTerms());
        ss << ",D=" << (p.isZero() ? 0 : p.totalDegree());

        return ss.str();
    }

    template <>
    inline std::string getStats(const storm::Polynomial& p) {
        std::stringstream ss;

        ss << "L=" << getStringComplexity(p);
        ss << ",M=" << (p.isConstant() ? 1 : p.nrTerms());
        ss << ",D=" << (p.isZero() ? 0 : p.totalDegree());

        return ss.str();
    }

#endif


}

}
}
