#pragma once


#pragma once

namespace storm {
namespace solver {
namespace gausselimination {

class PolynomialComplexity {
public:

    PolynomialComplexity(const storm::RawPolynomial& p) {
        totalDegree = p.isZero() ? 0 : p.totalDegree();
        numberOfTerms = p.nrTerms();
    }

    PolynomialComplexity(const storm::Polynomial& p) {
        totalDegree = p.isZero() ? 0 : p.totalDegree();
        numberOfTerms = p.isConstant() ? 1 : p.nrTerms();
    }

    std::size_t getTotalDegree() const {
        return totalDegree;
    }

    std::size_t getNumberOfTerms() const {
        return numberOfTerms;
    }

private:
    std::size_t totalDegree;
    std::size_t numberOfTerms;
};

}
}
}
