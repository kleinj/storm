#pragma once

#include "storm/storage/DenseMatrix.h"
#include "storm/storage/Permutation.h"
#include "storm/storage/SCCInfo.h"

namespace storm {
namespace storage {

template <typename ValueType>
class PermutatedDenseMatrix {
public:
    typedef typename DenseMatrix<ValueType>::value_type value_type;
    typedef typename DenseMatrix<ValueType>::index_type index_type;

    PermutatedDenseMatrix(DenseMatrix<ValueType>& underlyingMatrix,
                          storm::storage::Permutation::ptr permutation)
        : matrix(underlyingMatrix), permutationPtr(permutation), permutation(*permutation) {
    }

    PermutatedDenseMatrix(const PermutatedDenseMatrix&) = default;
    PermutatedDenseMatrix(PermutatedDenseMatrix&&) = default;

    const storm::storage::Permutation& getPermutation() const {
        return *permutationPtr;
    }

    ValueType& getEntry(std::size_t row, std::size_t col) {
        return matrix.getEntry(permutation.getInv(row), permutation.getInv(col));
    }

    const ValueType& getEntry(std::size_t row, std::size_t col) const {
        return matrix.getEntry(permutation.getInv(row), permutation.getInv(col));
    }

    ValueType& operator()(std::size_t row, std::size_t col) {
        return getEntry(row, col);
    }

    const ValueType& operator()(std::size_t row, std::size_t col) const {
        return getEntry(row, col);
    }

    std::size_t getRowCount() const {
        return matrix.getRowCount();
    }

    std::size_t getColumnCount() const {
        return matrix.getColumnCount();
    }

    template <typename T>
    friend std::ostream& operator<<(std::ostream& out, PermutatedDenseMatrix<T> const& matrix);

private:
    DenseMatrix<ValueType>& matrix;
    storm::storage::Permutation::ptr permutationPtr;
    const storm::storage::Permutation& permutation;
};

template <typename ValueType>
inline std::ostream& operator<<(std::ostream& out, PermutatedDenseMatrix<ValueType> const& matrix) {
    for (std::size_t i = 0; i < matrix.getRowCount(); i++) {
        for (std::size_t j = 0; j < matrix.getColumnCount(); j++) {
            out << "[" << i << "," << j << "] = " << matrix(i,j) << "\n";
        }
    }
    return out;
}

}
}

