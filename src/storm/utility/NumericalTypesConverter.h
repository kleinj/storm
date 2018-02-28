#pragma once

#include "storm/utility/macros.h"
#include "storm/exceptions/IllegalArgumentException.h"

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"

namespace storm {
namespace utility {

template <typename From, typename To>
struct NumericalTypesConverter {
    To convert(const From& from) {
        // default: use constructor from To
        return To(from);
    }
};

// ----------- specializations ---------------------------------

// identity (just pass const reference)

template <typename T>
struct NumericalTypesConverter<T,T> {
    const T& convert(const T& from) {
        return from;
    }
};

template <>
struct NumericalTypesConverter<storm::RationalFunction, storm::PlainRationalFunctionNoSimplify> {

    PlainRationalFunctionNoSimplify convert(const storm::RationalFunction& from) {
        storm::PlainRationalFunctionNoSimplify result(carl::computePolynomial(from.nominator()));
        result /= carl::computePolynomial(from.denominator());
        return result;
    }
};

template <>
struct NumericalTypesConverter<storm::FactorizedRationalFunctionNoSimplify, storm::RationalFunction> {
    storm::RationalFunction convert(const storm::FactorizedRationalFunctionNoSimplify& from) {
        storm::RationalFunction result(from.nominator(), from.denominator());
        return result;
    }
};

template <>
struct NumericalTypesConverter<storm::RationalFunction, storm::PlainRationalFunction> {
    storm::PlainRationalFunction convert(const storm::RationalFunction& from) {
        storm::PlainRationalFunction result(carl::computePolynomial(from.nominator()));
        result /= carl::computePolynomial(from.denominator());
        return result;
    }
};

template <>
struct NumericalTypesConverter<storm::RationalFunction, storm::RawPolynomial> {
    storm::RawPolynomial convert(const storm::RationalFunction& from) {
        if (from.isConstant()) {
            return storm::RawPolynomial(from.constantPart());
        }
        storm::RawPolynomial result(carl::computePolynomial(from.nominator()));
        STORM_LOG_THROW(from.denominator().isConstant(), storm::exceptions::IllegalArgumentException, "Can not convert rational function with non-constant denominator to polynomial: " << from);
        if (!from.denominator().isOne()) {
            result /= from.denominator().constantPart();
        }
        return result;
    }
};

template <>
struct NumericalTypesConverter<storm::PlainRationalFunctionNoSimplify, storm::RationalFunctionCoefficient> {
    storm::RationalFunctionCoefficient convert(const storm::PlainRationalFunctionNoSimplify& from) {
        STORM_LOG_THROW(from.isConstant(), storm::exceptions::IllegalArgumentException, "Can not convert non-constant rational function to rational number: " << from);
        return storm::RationalFunctionCoefficient(from.constantPart());
    }
};

template <>
struct NumericalTypesConverter<storm::PlainRationalFunction, storm::RationalFunctionCoefficient> {
    storm::RationalFunctionCoefficient convert(const storm::PlainRationalFunction& from) {
        STORM_LOG_THROW(from.isConstant(), storm::exceptions::IllegalArgumentException, "Can not convert non-constant rational function to rational number: " << from);
        return storm::RationalFunctionCoefficient(from.constantPart());
    }
};

template <>
struct NumericalTypesConverter<storm::RationalFunction, storm::RationalFunctionCoefficient> {
    storm::RationalFunctionCoefficient convert(const storm::RationalFunction& from) {
        STORM_LOG_THROW(from.isConstant(), storm::exceptions::IllegalArgumentException, "Can not convert non-constant rational function to rational number: " << from);
        return storm::RationalFunctionCoefficient(from.constantPart());
    }
};

template <>
struct NumericalTypesConverter<storm::FactorizedRationalFunctionNoSimplify, storm::RationalFunctionCoefficient> {
    storm::RationalFunctionCoefficient convert(const storm::FactorizedRationalFunctionNoSimplify& from) {
        STORM_LOG_THROW(from.isConstant(), storm::exceptions::IllegalArgumentException, "Can not convert non-constant rational function to rational number: " << from);
        return storm::RationalFunctionCoefficient(from.constantPart());
    }
};

// RationalFunction -> FactorizedRationalFunctionNoSimplify

template <>
struct NumericalTypesConverter<storm::RationalFunction, storm::FactorizedRationalFunctionNoSimplify> {
    storm::FactorizedRationalFunctionNoSimplify convert(const storm::RationalFunction& from) {
        return storm::FactorizedRationalFunctionNoSimplify(from.nominator(), from.denominator());
    }
};

// RawPolynomial -> RationalFunction

template <>
struct NumericalTypesConverter<storm::RawPolynomial, storm::RationalFunction> {
    NumericalTypesConverter() : cache(new carl::Cache<carl::PolynomialFactorizationPair<RawPolynomial>>()) {
    }

    storm::RationalFunction convert(const storm::RawPolynomial& from) {
        return storm::RationalFunction(storm::Polynomial(from, cache));
    }

private:
    std::shared_ptr<carl::Cache<carl::PolynomialFactorizationPair<RawPolynomial>>> cache;
};


// PlainRationalFunctionNoSimplify -> RationalFunction

template <>
struct NumericalTypesConverter<storm::PlainRationalFunctionNoSimplify, storm::RationalFunction> {
    NumericalTypesConverter() : cache(new carl::Cache<carl::PolynomialFactorizationPair<RawPolynomial>>()) {
    }

    storm::RationalFunction convert(const storm::PlainRationalFunctionNoSimplify& from) {
        return storm::RationalFunction(storm::Polynomial(from.nominator(), cache), storm::Polynomial(from.denominator(), cache));
    }

private:
    std::shared_ptr<carl::Cache<carl::PolynomialFactorizationPair<RawPolynomial>>> cache;
};


// PlainRationalFunction -> RationalFunction

template <>
struct NumericalTypesConverter<storm::PlainRationalFunction, storm::RationalFunction> {
    NumericalTypesConverter() : cache(new carl::Cache<carl::PolynomialFactorizationPair<RawPolynomial>>()) {
    }

    storm::RationalFunction convert(const storm::PlainRationalFunction& from) {
        return storm::RationalFunction(storm::Polynomial(from.nominator(), cache), storm::Polynomial(from.denominator(), cache));
    }

private:
    std::shared_ptr<carl::Cache<carl::PolynomialFactorizationPair<RawPolynomial>>> cache;
};


}
}
