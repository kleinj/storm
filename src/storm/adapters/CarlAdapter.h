#pragma once

#include "storm/adapters/NumberAdapter.h"

#include <carl/core/MultivariatePolynomial.h>
#include <carl/core/RationalFunction.h>
#include <carl/core/VariablePool.h>
#include <carl/core/FactorizedPolynomial.h>
#include <carl/core/Relation.h>
#include <carl/core/SimpleConstraint.h>
#include <carl/util/stringparser.h>

namespace carl {
    // Define hash values for all polynomials and rational function.
    template<typename C, typename O, typename P>
    inline size_t hash_value(carl::MultivariatePolynomial<C,O,P> const& p) {
        std::hash<carl::MultivariatePolynomial<C,O,P>> h;
        return h(p);
    }
    
    template<typename Pol>
    inline size_t hash_value(carl::FactorizedPolynomial<Pol> const& p) {
        std::hash<FactorizedPolynomial<Pol>> h;
        return h(p);
    }

    template<typename Pol, bool AutoSimplify>
    inline size_t hash_value(carl::RationalFunction<Pol, AutoSimplify> const& f)  {
        std::hash<Pol> h;
        return h(f.nominator()) ^ h(f.denominator());
    }
    
    template<typename Number>
    inline size_t hash_value(carl::Interval<Number> const& i) {
        std::hash<Interval<Number>> h;
        return h(i);
    }

}

namespace storm {
    typedef carl::Variable RationalFunctionVariable;

#if defined(STORM_HAVE_CLN) && defined(STORM_USE_CLN_RF)
    typedef cln::cl_RA RationalFunctionCoefficient;
#elif defined(STORM_HAVE_GMP) && !defined(STORM_USE_CLN_RF)
    typedef mpq_class RationalFunctionCoefficient;
#elif defined(STORM_USE_CLN_RF)
#error CLN is to be used, but is not available.
#else
#error GMP is to be used, but is not available.
#endif
    
    typedef carl::MultivariatePolynomial<RationalFunctionCoefficient> RawPolynomial;
    typedef carl::FactorizedPolynomial<RawPolynomial> Polynomial;
	typedef carl::Cache<carl::PolynomialFactorizationPair<RawPolynomial>> RawPolynomialCache;
    typedef carl::Relation CompareRelation;
    
    typedef carl::RationalFunction<Polynomial, true> RationalFunction;
    typedef carl::Interval<double> Interval;
    template<typename T> using ArithConstraint = carl::SimpleConstraint<T>;
}
