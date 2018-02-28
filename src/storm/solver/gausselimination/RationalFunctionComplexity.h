#pragma once

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/solver/gausselimination/PolynomialComplexity.h"
#include <string>

namespace storm {
namespace solver {
namespace gausselimination {


class RationalFunctionComplexity {
public:

    template <typename P, bool AS>
    RationalFunctionComplexity(const typename carl::RationalFunction<P,AS>& f)
     :  numeratorComplexity(f.nominator()),
        denominatorComplexity(f.denominator()) {
    }

    const PolynomialComplexity& numerator() const {
        return numeratorComplexity;
    }

    const PolynomialComplexity& denominator() const {
        return denominatorComplexity;
    }

    std::size_t maxReductionDegree(const RationalFunctionComplexity& reducedFunction) const {
        std::size_t reduction_num = numerator().getTotalDegree() - reducedFunction.numerator().getTotalDegree();
        std::size_t reduction_denom = denominator().getTotalDegree() - reducedFunction.denominator().getTotalDegree();

        if (reduction_num > reduction_denom)
            return reduction_num;
        return reduction_denom;
    }

    std::size_t maxReductionTerms(const RationalFunctionComplexity& reducedFunction) const {
        std::size_t reduction_num = numerator().getNumberOfTerms() - reducedFunction.numerator().getNumberOfTerms();
        std::size_t reduction_denom = denominator().getNumberOfTerms() - reducedFunction.denominator().getNumberOfTerms();

        if (reduction_num > reduction_denom)
            return reduction_num;
        return reduction_denom;
    }

    template <typename VT1, typename VT2>
    static std::string reportMaxReduction(const std::vector<VT1>& a, const std::vector<VT2>& b) {
        assert(a.size() == b.size());
        std::size_t maxRedDegree = 0;
        std::size_t maxRedTerms = 0;
        for (std::size_t i = 0; i < a.size(); i++) {
            RationalFunctionComplexity after(a.at(i));
            RationalFunctionComplexity before(b.at(i));

            std::size_t redDegree = before.maxReductionDegree(after);
            std::size_t redTerms = before.maxReductionTerms(after);

            if (redDegree > maxRedDegree) {
                maxRedDegree = redDegree;
            }
            if (redTerms > maxRedTerms) {
                maxRedTerms = redTerms;
            }
        }

        return "Max reduction: Degree = " + std::to_string(maxRedDegree) + ", Terms = " + std::to_string(maxRedTerms);

    }

private:
    PolynomialComplexity numeratorComplexity;
    PolynomialComplexity denominatorComplexity;
};

}
}
}
