#pragma once

#include "storm/storage/SparseMatrix.h"
#include <iostream>
#include <sstream>

namespace storm {
    class RationalFunctionStatistics {
    public:
        template <typename ValueType>
        static void analyzeSLE(const storm::storage::SparseMatrix<ValueType>& matrix, const std::vector<ValueType>& b) {
            (void) matrix; // IGNORE
            (void) b; // IGNORE
            return;
        }

    private:
        static void countMonomials(const RawPolynomial& poly, std::map<carl::Monomial::Arg, std::size_t>& monomials) {
            for (auto& term : poly) {
                const carl::Monomial::Arg monomial = term.monomial();
                if (!monomial)
                    continue;
                monomials[monomial]++;
            }
        }
    };

    template <>
    inline void RationalFunctionStatistics::analyzeSLE<RationalFunction>(const storm::storage::SparseMatrix<RationalFunction>& matrix, const std::vector<RationalFunction>& b) {
        std::size_t numberOfParametricStates = 0;
        std::size_t numberOfParametricTransitions = 0;
        std::size_t numberOfNnzTransitions = 0;
        std::size_t numberOfParametricBEntries = 0;

        std::set<RationalFunctionVariable> variables;
        std::map<carl::Monomial::Arg, std::size_t> monomials;

        std::stringstream out;

        for (std::size_t s = 0; s < matrix.getRowCount(); s++) {
            bool stateIsParametric = false;
            for (auto& entries : matrix.getRow(s)) {
                const RationalFunction& p = entries.getValue();
                if (!p.isConstant()) {
                    stateIsParametric = true;
                    numberOfParametricTransitions++;

                    p.gatherVariables(variables);
                    if (!p.nominatorAsPolynomial().isConstant()) {
                        countMonomials(p.nominatorAsPolynomial().polynomial(), monomials);
                    }
                    if (!p.denominatorAsPolynomial().isConstant()) {
                        countMonomials(p.denominatorAsPolynomial().polynomial(), monomials);
                    }
                }
                numberOfNnzTransitions++;
            }

            if (stateIsParametric) {
                numberOfParametricStates++;
            }
        }

        for (const RationalFunction& p : b) {
            if (!p.isConstant()) {
                numberOfParametricBEntries++;

                p.gatherVariables(variables);
                if (!p.nominatorAsPolynomial().isConstant()) {
                    countMonomials(p.nominatorAsPolynomial().polynomial(), monomials);
                }
                if (!p.denominatorAsPolynomial().isConstant()) {
                    countMonomials(p.denominatorAsPolynomial().polynomial(), monomials);
                }
            }
        }

        out << "\nMatrix statistics:\n";
        out << " Number of parametric states:               " << numberOfParametricStates << " of " << matrix.getRowCount() << "\n";
        out << " Number of parametric transitions:          " << numberOfParametricTransitions << " of " << numberOfNnzTransitions << "\n";
        out << " Number of parametric entries in b-vector:  " << numberOfParametricBEntries << " of " << matrix.getRowCount() << "\n";
        out << " Variables: " << variables << "\n";

        out << " Monomial count:\n";
        for (auto& entry : monomials) {
            out << "   " << *(entry.first) << ": " << entry.second << "\n";
        }

        STORM_LOG_INFO(out.str());
    }

}
