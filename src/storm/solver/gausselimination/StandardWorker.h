#pragma once

#include <vector>

namespace storm {
namespace solver {
namespace gausselimination {

template <typename ValueType>
class StandardWorker {
public:
    typedef ValueType SparseValueType;
    typedef ValueType DenseValueType;
    typedef ValueType IntermediateValueType;

    template <class Matrix>
    void prepare(Matrix& matrix, boost::optional<const storm::storage::SCCInfo&> sccInfo = boost::none) {
        (void)matrix; // IGNORE
        (void)sccInfo; // IGNORE
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
    void startTriangulation(Matrix& matrix, std::size_t pivotRow) {
        (void)matrix; // IGNORE
        (void)pivotRow; // IGNORE
    }

    template <class Matrix>
    void endTriangulation(Matrix& matrix, std::size_t k, const DenseValueType& a_kk) {
        (void)matrix; // IGNORE
        (void)k; // IGNORE
        (void)a_kk; // IGNORE
    }

    template <class Matrix>
    void startTriangulation(Matrix& matrix, std::size_t k, std::size_t i,
                            const DenseValueType& a_kk,
                            const DenseValueType& a_ik) {
        pivot = a_ik / a_kk;
        if (debug) {
            std::cout << "Pivot for row " << i << " (against row " << k << ") = " << pivot << " from " << a_ik << " and " << a_kk << "\n";
        }
    }

    template <class Matrix>
    DenseValueType doTriangulation(Matrix& matrix, std::size_t k, std::size_t i, std::size_t j, const DenseValueType& a_kk, const DenseValueType& a_ij, const DenseValueType& a_ik, const DenseValueType& a_kj) {
        if (debug) std::cout << "Triangulate for row " << i << ", col " << j << " (against row " << k << ")\n";
        return a_ij - (a_kj * pivot);
    }

    template <class Matrix>
    void startSolve(const Matrix& triangularMatrix, std::vector<IntermediateValueType>& vector) {
        (void)triangularMatrix;  // IGNORE
        (void)vector;  // IGNORE
    }

    template <class Matrix>
    void startSolve_sparse(const Matrix& triangularMatrix, std::vector<IntermediateValueType>& vector) {
        (void)triangularMatrix;  // IGNORE
        (void)vector;  // IGNORE
    }

    template <class Matrix>
    ValueType solve(const Matrix& triangularMatrix, std::size_t row, std::vector<ValueType>& vector) {
        typename Matrix::value_type result;
        for (std::size_t j = row + 1; j < triangularMatrix.getColumnCount() - 1; j++) {
            result -= triangularMatrix(row, j) * vector.at(j);
        }
        result += triangularMatrix(row, triangularMatrix.getColumnCount()-1);
        result /= triangularMatrix(row, row);

        return result;
    }

    template <class Matrix>
    ValueType solve_sparse(const Matrix& triangularMatrix, std::size_t row, std::vector<ValueType>& vector) {
        typename Matrix::value_type result;
        std::size_t b_col = triangularMatrix.getColumnCount() - 1;
        for (auto& entry : triangularMatrix.row_col_range(row, row+1, b_col-1)) {
            std::size_t j = entry.getColumn();
            result -= entry.getValue() * vector.at(j);
        }
        result += triangularMatrix.getValueRight(row, b_col);
        result /= triangularMatrix.getValueLeft(row, row);
        return result;
    }

    template <class Matrix>
    void postprocessSolution(const Matrix& triangularMatrix, std::vector<IntermediateValueType>& vector) {
        (void)triangularMatrix; // IGNORE
        (void)vector; // IGNORE
    }

    static std::string getMethodName() {
        return std::string("Standard [all values: ") + storm::utility::NumericalTypes::getTypeName<ValueType>() + "]";
    }

    std::string getMethodNameAndConfiguration() {
        return getMethodName();
    }

private:
    ValueType pivot;
    bool debug = false;
};

}
}
}
