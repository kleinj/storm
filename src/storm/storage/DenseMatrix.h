#pragma once

#include <vector>
#include <cassert>

#include "storm/utility/NumericalTypesConverter.h"

namespace storm {
namespace storage {

template <typename ValueType>
class DenseMatrix {
public:
    typedef ValueType value_type;
    typedef std::size_t index_type;

    // TODO rows*cols can overflow
    DenseMatrix(std::size_t rows, std::size_t cols) : rows(rows), cols(cols), matrix(rows * cols) {
    }

    ValueType& getEntry(std::size_t row, std::size_t col) {
        assert(row < rows);
        assert(col < cols);
        std::size_t index = row * cols + col;
        return matrix.at(index);
    }

    const ValueType& getEntry(std::size_t row, std::size_t col) const {
        assert(row < rows);
        assert(col < cols);
        std::size_t index = row * cols + col;
        return matrix.at(index);
    }

    ValueType& operator()(std::size_t row, std::size_t col) {
        return getEntry(row, col);
    }

    const ValueType& operator()(std::size_t row, std::size_t col) const {
        return getEntry(row, col);
    }

    std::size_t getRowCount() const {
        return rows;
    }

    std::size_t getColumnCount() const {
        return cols;
    }

    template <typename SparseValueType, typename Converter = storm::utility::NumericalTypesConverter<SparseValueType,ValueType>>
    static DenseMatrix augmentedMatrix(const storm::storage::SparseMatrix<SparseValueType>& A, const std::vector<SparseValueType>& b) {
        assert(A.getRowCount() == A.getColumnCount());
        assert(A.getRowCount() == b.size());

        Converter converter;

        DenseMatrix matrix(A.getRowCount(), A.getColumnCount()+1);
        for (std::size_t i = 0; i< A.getRowCount(); i++) {
            for (auto& entry : A.getRow(i)) {
                matrix(i, entry.getColumn()) = converter.convert(entry.getValue());
            }

            matrix(i, A.getColumnCount()) = converter.convert(b.at(i));
        }

        return matrix;
    }

    template <typename T>
    friend std::ostream& operator<<(std::ostream& out, DenseMatrix<T> const& matrix);

private:
    std::size_t rows;
    std::size_t cols;
    std::vector<ValueType> matrix;
};

template <typename ValueType>
inline std::ostream& operator<<(std::ostream& out, DenseMatrix<ValueType> const& matrix) {
    for (std::size_t i = 0; i < matrix.getRowCount(); i++) {
        for (std::size_t j = 0; j < matrix.getColumnCount(); j++) {
            out << "[" << i << "," << j << "] = " << matrix(i,j) << "\n";
        }
    }
    return out;
}


}
}

