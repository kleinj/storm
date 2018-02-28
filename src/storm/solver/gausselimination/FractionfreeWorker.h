#pragma once

#include <vector>

#include "storm/utility/NumericalTypes.h"
#include "storm/utility/NumericalTypesConverter.h"
#include "storm/settings/modules/GaussEliminationSettings.h"

namespace storm {
namespace solver {
namespace gausselimination {

template <typename SparseValueTypeT, typename IntermediateRationalTypeT, typename NoFractionTypeT>
class FractionfreeWorker {
public:
    typedef SparseValueTypeT SparseValueType;
    typedef NoFractionTypeT DenseValueType;
    typedef IntermediateRationalTypeT IntermediateValueType;

    template <class Matrix>
    void prepare(Matrix& matrix, boost::optional<const storm::storage::SCCInfo&> sccInfo = boost::none) {
        (void)matrix; // IGNORE
        sccInfoOpt = sccInfo;

        handleSCCsSeparately = (bool)sccInfoOpt;
        handleSCCsSeparately &= storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().isIndividualSCCsSet();
        if (handleSCCsSeparately) {
            useGlobalDenominator = !storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().isSCCDenominatorsSet();
        }

        solnNumerator.resize(matrix.getRowCount());
    }

    template <class Matrix>
    void startTriangulation(Matrix& matrix) {
        (void)matrix; // IGNORE
    }

    template <class Matrix>
    void endTriangulation(Matrix& matrix) {
        (void)matrix; // IGNORE
    }

    template <class Matrix>
    void startTriangulation(Matrix& matrix, std::size_t k) {
        (void)matrix; // IGNORE
        (void)k; // IGNORE
    }

    template <class Matrix>
    void endTriangulation(Matrix& matrix, std::size_t k, const DenseValueType& a_kk) {
        (void)matrix; // IGNORE
        (void)k; // IGNORE
        cancelValueForRow = a_kk;
    }

    template <class Matrix>
    void startTriangulation(Matrix& matrix, std::size_t k, std::size_t i,
                            const DenseValueType& a_kk,
                            const DenseValueType& a_ik) {
        (void)matrix; // IGNORE
        (void)k; // IGNORE
        (void)i; // IGNORE
        (void)a_kk; // IGNORE
        (void)a_ik; // IGNORE
    }

    template <class Matrix>
    DenseValueType doTriangulation(Matrix& matrix, std::size_t k, std::size_t i, std::size_t j,
                                   const DenseValueType& a_kk,
                                   const DenseValueType& a_ij,
                                   const DenseValueType& a_ik,
                                   const DenseValueType& a_kj) {
        if (debug) {
            std::cout << "(" << k << "," << i << "," << j << "): before triangulate: " << a_ij << "\n";
            std::cout << "(" << k << "," << i << "," << j << "): mult with: " << a_kk << "\n";
        }

        auto v = a_ij * a_kk; // a'_i,j = a_i,j * a_k,k

        if (debug) std::cout << "(" << k << "," << i << "," << j << "): subtract: " << (a_kj * a_ik) << "\n";
        v -= a_kj * a_ik;  // a'_i,j -= a_k,j * a_i,k

        if (debug) std::cout << "(" << k << "," << i << "," << j << "): after triangulate: " << v << "\n";

        bool cancel = (k >= 1);
        if (handleSCCsSeparately) {
            // if we handle each SCC separately, we can only cancel if the pivot row
            // is not the first one (index = 0) *in this SCC*
            cancel = (sccInfoOpt.get().getIndexInSCC(k) >= 1);
        }

        // cancel
        if (cancel) {
            if (debug) std::cout << "(" << k << "," << i << "," << j << "): before cancel: " << v << " with " << cancelValueForRow.get() << "\n";

            // TODO /=
            v = v / cancelValueForRow.get();

            if (debug) std::cout << "(" << k << "," << i << "," << j << "): after cancel with " << cancelValueForRow.get() << ": " << v << "\n";
        }

        return v;
    }

    template <class Matrix>
    void startSolve(const Matrix& triangularMatrix, std::vector<IntermediateValueType>& vector) {
        (void)vector;  // IGNORE

        typedef std::shared_ptr<DenseValueType> Value_ptr;

        // prepare denominators
        if (!useGlobalDenominator) {
            solnDenominator.resize(triangularMatrix.getRowCount());
            assert(handleSCCsSeparately);
            std::size_t last_state = triangularMatrix.getRowCount()-1;
            solnDenominator.at(last_state) = Value_ptr(new DenseValueType(triangularMatrix(last_state, last_state)));
            if (debug) {
                std::cout << "denom[" << last_state << "] = " << *solnDenominator.at(last_state) << "\n";
            }

            auto& sccInfo = sccInfoOpt.get();
            for (std::size_t i = last_state; i > 0; i--) {
                std::size_t state = i - 1;
                // i runs from last_state ... 1
                // => state from last_state-1 to ... 0
                if (sccInfo.getSCC(state) == sccInfo.getSCC(state + 1)) {
                    // same SCC, copy denominator
                    solnDenominator.at(state) = solnDenominator.at(state+1);
                } else {
                    // state is first state in new SCC, extend denominator
                    const DenseValueType& a = triangularMatrix(state, state);
                    solnDenominator.at(state) =
                            Value_ptr(new DenseValueType(std::move(a * *solnDenominator.at(state+1))));
                }
                if (debug) {
                    std::cout << "denom[" << state << "] = " << *solnDenominator.at(state) << "\n";
                }
            }

            extendForSCC.resize(sccInfo.getNumberOfSCCs(), DenseValueType(1));
        } else {
            if (handleSCCsSeparately) {
                // global denominator is the product of the
                // diagonal entries for all the *last* states in each SCC
                auto& sccInfo = sccInfoOpt.get();
                globalDenominator = DenseValueType(1);
                for (std::size_t scc = 0; scc < sccInfo.getNumberOfSCCs(); scc++) {
                    std::size_t last_in_scc = sccInfo.getLastStateInSCC(scc);
                    globalDenominator *= triangularMatrix(last_in_scc, last_in_scc);
                }
            } else {
                // global denominator is the last diagonal entry in the matrix
                globalDenominator = triangularMatrix(triangularMatrix.getRowCount()-1, triangularMatrix.getRowCount()-1);
            }

            if (debug) {
                std::cout << "Global denominator = " << globalDenominator;
            }
        }
    }

    template <class Matrix>
    void startSolve_sparse(const Matrix& triangularMatrix, std::vector<IntermediateValueType>& vector) {
        (void)vector;  // IGNORE

        typedef std::shared_ptr<DenseValueType> Value_ptr;

        // prepare denominators
        if (!useGlobalDenominator) {
            solnDenominator.resize(triangularMatrix.getRowCount());
            assert(handleSCCsSeparately);
            std::size_t last_state = triangularMatrix.getRowCount()-1;
            solnDenominator.at(last_state) = Value_ptr(new DenseValueType(triangularMatrix.getValueLeft(last_state, last_state)));
            if (debug) {
                std::cout << "denom[" << last_state << "] = " << *solnDenominator.at(last_state) << "\n";
            }

            auto& sccInfo = sccInfoOpt.get();
            for (std::size_t i = last_state; i > 0; i--) {
                std::size_t state = i - 1;
                // i runs from last_state ... 1
                // => state from last_state-1 to ... 0
                if (sccInfo.getSCC(state) == sccInfo.getSCC(state + 1)) {
                    // same SCC, copy denominator
                    solnDenominator.at(state) = solnDenominator.at(state+1);
                } else {
                    // state is first state in new SCC, extend denominator
                    const DenseValueType& a = triangularMatrix.getValueLeft(state, state);
                    solnDenominator.at(state) =
                            Value_ptr(new DenseValueType(std::move(a * *solnDenominator.at(state+1))));
                }
                if (debug) {
                    std::cout << "denom[" << state << "] = " << *solnDenominator.at(state) << "\n";
                }
            }

            extendForSCC.resize(sccInfo.getNumberOfSCCs(), DenseValueType(1));
        } else {
            if (handleSCCsSeparately) {
                // global denominator is the product of the
                // diagonal entries for all the *last* states in each SCC
                auto& sccInfo = sccInfoOpt.get();
                globalDenominator = DenseValueType(1);
                for (std::size_t scc = 0; scc < sccInfo.getNumberOfSCCs(); scc++) {
                    std::size_t last_in_scc = sccInfo.getLastStateInSCC(scc);
                    globalDenominator *= triangularMatrix.getValueLeft(last_in_scc, last_in_scc);
                }
            } else {
                // global denominator is the last diagonal entry in the matrix
                globalDenominator = triangularMatrix.getValueLeft(triangularMatrix.getRowCount()-1, triangularMatrix.getRowCount()-1);
            }

            if (debug) {
                std::cout << "Global denominator = " << globalDenominator;
            }
        }
    }

    template <class Matrix>
    IntermediateValueType solve(const Matrix& triangularMatrix, std::size_t row, std::vector<IntermediateValueType>& vector) {
        if (handleSCCsSeparately) {
            if (useGlobalDenominator) {
                return solve_with_global_denominator(triangularMatrix, row, vector);
            } else {
                return solve_with_scc_denominators(triangularMatrix, row, vector);
            }
        } else {
            return solve_with_global_denominator(triangularMatrix, row, vector);
        }
    }

    template <class Matrix>
    IntermediateValueType solve_sparse(const Matrix& triangularMatrix, std::size_t row, std::vector<IntermediateValueType>& vector) {
        if (handleSCCsSeparately) {
            if (useGlobalDenominator) {
                return solve_sparse_with_global_denominator(triangularMatrix, row, vector);
            } else {
                return solve_sparse_with_scc_denominators(triangularMatrix, row, vector);
            }
        } else {
            return solve_sparse_with_global_denominator(triangularMatrix, row, vector);
        }
    }

    template <class Matrix>
    IntermediateValueType solve_with_global_denominator(const Matrix& triangularMatrix, std::size_t row, std::vector<IntermediateValueType>& vector) {
        (void)vector;  // we ignore the solution vector, it will be filled in postprocessing

        DenseValueType result;

        for (std::size_t j = row + 1; j < triangularMatrix.getColumnCount() - 1; j++) {
            result -= triangularMatrix(row, j) * solnNumerator.at(j);
        }
        const DenseValueType& b = triangularMatrix(row, triangularMatrix.getColumnCount()-1);

        bool dontExtendAndCancel = false;
        if (row == triangularMatrix.getRowCount()-1 && !handleSCCsSeparately) {
            // last row, don't extend / cancel
            dontExtendAndCancel = true;
        }
        if (dontExtendAndCancel) {
            result += b;
        } else {
            // extend b with the global denominator and add to result
            result += b * globalDenominator;

            const DenseValueType& cancelValue = triangularMatrix(row, row);
            if (debug) {
                std::cout << "cancelValue = " << cancelValue << "\n";
                std::cout << "result = " << result << "\n";
            }
            // we know we can cancel by cancelValue
            result = result / cancelValue;
        }

        if (debug) std::cout << "result / cancelValue = " << result << "\n";

        // store result numerator
        solnNumerator.at(row) = std::move(result);

        // return placeholder result
        return IntermediateValueType();
    }

    template <class Matrix>
    IntermediateValueType solve_sparse_with_global_denominator(const Matrix& triangularMatrix, std::size_t row, std::vector<IntermediateValueType>& vector) {
        (void)vector;  // we ignore the solution vector, it will be filled in postprocessing

        DenseValueType result;

        std::size_t b_col = triangularMatrix.getColumnCount()-1;
        for (auto& entry : triangularMatrix.row_col_range(row, row+1, b_col-1)) {
            std::size_t j = entry.getColumn();
            result -= entry.getValue() * solnNumerator.at(j);
        }
        const DenseValueType& b = triangularMatrix.getValueRight(row, b_col);

        bool dontExtendAndCancel = false;
        if (row == triangularMatrix.getRowCount()-1 && !handleSCCsSeparately) {
            // last row, don't extend / cancel
            dontExtendAndCancel = true;
        }
        if (dontExtendAndCancel) {
            result += b;
        } else {
            // extend b with the global denominator and add to result
            result += b * globalDenominator;

            const DenseValueType& cancelValue = triangularMatrix.getValueLeft(row, row);
            if (debug) {
                std::cout << "cancelValue = " << cancelValue << "\n";
                std::cout << "result = " << result << "\n";
            }
            // we know we can cancel by cancelValue
            result = result / cancelValue;
        }


        if (debug) std::cout << "result / cancelValue = " << result << "\n";

        // store result numerator
        solnNumerator.at(row) = std::move(result);

        // return placeholder result
        return IntermediateValueType();
    }


    template <class Matrix>
    IntermediateValueType solve_with_scc_denominators(const Matrix& triangularMatrix, std::size_t row, std::vector<IntermediateValueType>& vector) {
        (void)vector;  // we ignore the solution vector, it will be filled in postprocessing

        DenseValueType result;

        auto& sccInfo = sccInfoOpt.get();
        std::size_t cur_scc = sccInfo.getSCC(row);
        bool last_state_in_scc = sccInfo.getLastStateInSCC(cur_scc) == row;
        if (last_state_in_scc) {
            const DenseValueType& factorForCurSCC = triangularMatrix(row, row);
            for (std::size_t i = cur_scc + 1; i < sccInfo.getNumberOfSCCs(); i++) {
                extendForSCC.at(i) *= factorForCurSCC;
            }
        }

        for (std::size_t j = row + 1; j < triangularMatrix.getColumnCount() - 1; j++) {
            DenseValueType a_times_x = triangularMatrix(row, j) * solnNumerator.at(j);
            // grab the value for extending for the target SCC and extend as necessary
            std::size_t target_scc = sccInfo.getSCC(j);
            if (target_scc != cur_scc) {
                const DenseValueType& extend = extendForSCC.at(sccInfo.getSCC(j));
//                assert(extend == (*solnDenominator.at(row) / *solnDenominator.at(j)));
                a_times_x *= extend;
            }
            result -= a_times_x;
        }

        // extend b with the denominator and add to result
        const DenseValueType& b = triangularMatrix(row, triangularMatrix.getColumnCount()-1);
        result += b * *solnDenominator.at(row);

        const DenseValueType& cancelValue = triangularMatrix(row, row);
        if (debug) {
            std::cout << "cancelValue = " << cancelValue << "\n";
            std::cout << "result = " << result << "\n";
        }
        // we know we can cancel by cancelValue
        result = result / cancelValue;

        if (debug) std::cout << "result / cancelValue = " << result << "\n";

        // store result numerator
        solnNumerator.at(row) = std::move(result);

        // return placeholder result
        return IntermediateValueType();
    }

    template <class Matrix>
    IntermediateValueType solve_sparse_with_scc_denominators(const Matrix& triangularMatrix, std::size_t row, std::vector<IntermediateValueType>& vector) {
        (void)vector;  // we ignore the solution vector, it will be filled in postprocessing

        DenseValueType result;

        auto& sccInfo = sccInfoOpt.get();
        std::size_t cur_scc = sccInfo.getSCC(row);
        bool last_state_in_scc = sccInfo.getLastStateInSCC(cur_scc) == row;
        if (last_state_in_scc) {
            const DenseValueType& factorForCurSCC = triangularMatrix.getValueLeft(row, row);
            for (std::size_t i = cur_scc + 1; i < sccInfo.getNumberOfSCCs(); i++) {
                extendForSCC.at(i) *= factorForCurSCC;
            }
        }

        std::size_t b_col = triangularMatrix.getColumnCount()-1;
        for (auto& entry : triangularMatrix.row_col_range(row, row+1, b_col-1)) {
            std::size_t j = entry.getColumn();
            DenseValueType a_times_x = entry.getValue() * solnNumerator.at(j);
            // grab the value for extending for the target SCC and extend as necessary
            std::size_t target_scc = sccInfo.getSCC(j);
            if (target_scc != cur_scc) {
                const DenseValueType& extend = extendForSCC.at(sccInfo.getSCC(j));
//                assert(extend == (*solnDenominator.at(row) / *solnDenominator.at(j)));
                a_times_x *= extend;
            }
            result -= a_times_x;
        }

        // extend b with the denominator and add to result
        const DenseValueType& b = triangularMatrix.getValueRight(row, b_col);
        result += b * *solnDenominator.at(row);

        const DenseValueType& cancelValue = triangularMatrix.getValueLeft(row, row);
        if (debug) {
            std::cout << "cancelValue = " << cancelValue << "\n";
            std::cout << "result = " << result << "\n";
        }
        // we know we can cancel by cancelValue
        result = result / cancelValue;

        if (debug) std::cout << "result / cancelValue = " << result << "\n";

        // store result numerator
        solnNumerator.at(row) = std::move(result);

        // return placeholder result
        return IntermediateValueType();
    }

    template <class Matrix>
    void postprocessSolution(const Matrix& triangularMatrix, std::vector<IntermediateValueType>& vector) {
        std::size_t last_row = triangularMatrix.getRowCount() - 1;
        if (useGlobalDenominator) {
            for (std::size_t i = 0; i <= last_row; i++) {
                // store actual value in solution vector
                vector.at(i) = std::move(IntermediateValueType(solnNumerator.at(i), globalDenominator));
            }
        } else {
            // we have no global denominator, instead we have stored the required denominators
            // in an array
            for (std::size_t i = 0; i <= last_row; i++) {
                // store actual value in solution vector
                vector.at(i) = std::move(IntermediateValueType(solnNumerator.at(i), *solnDenominator.at(i)));
            }
        }
    }

    static std::string getMethodName() {
        return std::string("Fractionfree Method [matrix entries: ")
                + storm::utility::NumericalTypes::getTypeName<DenseValueType>()
                + ", back-substitution: "
                + storm::utility::NumericalTypes::getTypeName<IntermediateValueType>()
                + ", final result: "
                + storm::utility::NumericalTypes::getTypeName<SparseValueType>()
                + "]";
    }

    std::string getMethodNameAndConfiguration() {
        std::string result = getMethodName();
        if (handleSCCsSeparately) {
            result += " handling SCCs separately";
        }
        if (useGlobalDenominator) {
            result += " with global denominator";
        } else {
            result += " with per-SCC denominator";
        }

        return result;
    }


private:
    bool debug = false;
    boost::optional<const storm::storage::SCCInfo&> sccInfoOpt;
    bool handleSCCsSeparately;
    bool useGlobalDenominator = true;
    storm::utility::NumericalTypesConverter<IntermediateValueType, SparseValueType> converterSparse;
    storm::utility::NumericalTypesConverter<DenseValueType, IntermediateValueType> converterIntermediate;
    boost::optional<const DenseValueType&> cancelValueForRow;
    // for backpropagation
    std::vector<DenseValueType> solnNumerator;
    // for SCC backpropagation
    std::vector<std::shared_ptr<DenseValueType>> solnDenominator;
    std::vector<DenseValueType> extendForSCC;
    DenseValueType curDenominator;
    DenseValueType globalDenominator;
};

}
}
}
