#pragma once

namespace storm {
namespace storage {

template <typename Matrix>
class SparseAccessForDenseMatrix {
public:
    class VirtualMatrixEntry {
    public:
        typedef typename Matrix::index_type index_type;
        typedef typename Matrix::value_type value_type;

        /*!
         * Constructs a virtual matrix entry with the given column.
         *
         * @param matrix The underlying dense matrix.
         * @param column The column of the matrix entry.

         */
        VirtualMatrixEntry(const Matrix& matrix, index_type row, index_type column) : matrix(matrix), row(row), column(column) {}

        VirtualMatrixEntry() = delete;
        VirtualMatrixEntry(VirtualMatrixEntry const& other) = default;
        VirtualMatrixEntry& operator=(VirtualMatrixEntry const& other) = default;
        VirtualMatrixEntry(VirtualMatrixEntry&& other) = default;
        VirtualMatrixEntry& operator=(VirtualMatrixEntry&& other) = default;

        /*!
         * Retrieves the column of the matrix entry.
         *
         * @return The column of the matrix entry.
         */
        index_type getColumn() const {
            return column;
        }

        /*!
         * Retrieves the value of the matrix entry.
         *
         * @return The value of the matrix entry.
         */
        value_type const& getValue() const {
            return matrix(row, column);
        }

        bool operator==(const VirtualMatrixEntry& other) const {
            return column == other.column &&
                   &matrix == &other.matrix &&
                   row == other.row;
        }

        bool operator!=(const VirtualMatrixEntry& other) const {
            return !(*this == other);
        }

        void next() {
            column++;
        }

    private:
        const Matrix& matrix;
        index_type row;
        index_type column;
    };

    typedef typename Matrix::index_type index_type;
    typedef typename Matrix::value_type value_type;
    typedef VirtualMatrixEntry entry_type;

    class const_iterator {
    public:
        typedef VirtualMatrixEntry value_type;
        typedef VirtualMatrixEntry& reference;
        typedef VirtualMatrixEntry* pointer;
        typedef int difference_type;
        typedef std::forward_iterator_tag iterator_category;

        const_iterator(const Matrix& matrix, index_type row, index_type column) :
            current(matrix, row, column) { }

        const_iterator(const Matrix& matrix, index_type row) :
            current(matrix, row, matrix.getColumnCount()) { }

        const_iterator& operator++() {
            // prefix operator
            current.next();
            return *this;
        }

        const_iterator operator++(int) {
            // postfix operator
            const_iterator i = *this;
            current.next();
            return i;
        }

        const reference operator*() {
            return current;
        }

        const pointer operator->() {
            return &current;
        }

        bool operator==(const const_iterator& other) {
            return current == other.current;
        }

        bool operator!=(const const_iterator& other) {
            return current != other.current;
        }

    private:
        VirtualMatrixEntry current; // the current entry
    };

    class column_range {
    public:
        column_range(const Matrix& matrix, index_type row, index_type first_col, index_type last_col) :
            matrix(matrix), row(row), first_col(first_col), last_col(last_col) {
        }

        const_iterator begin() const {
            return const_iterator(matrix, row, first_col);
        }

        const_iterator end() const {
            return const_iterator(matrix, row, last_col + 1);
        }

    private:
        const Matrix& matrix;
        index_type row;
        index_type first_col;
        index_type last_col;
    };

    SparseAccessForDenseMatrix(const Matrix& matrix) : matrix(matrix) {}

    /*!
     * Returns a const reference to the value at the given row / column
     * (using linear search in the sparse entries).
     * If the value is zero (and thus not stored in the matrix),
     * returns a const ref to a dummy zero value.
     */
    value_type const& getValue(index_type rowIndex, index_type colIndex) const {
        return matrix(rowIndex, colIndex);
    }

    /*!
     * Returns a const reference to the value at the given row / column.
     * If the column is not to the left or right of the existing matrix entries (and thus 0)
     * or exactly the first / last non-zero matrix entry, an exception may be thrown.
     */
    value_type const& getValueFast(index_type rowIndex, index_type colIndex) const {
        return matrix(rowIndex, colIndex);
    }

    /*!
     * Returns a const reference to the value at the given row / column.
     * If the column is not to the left of the existing matrix entries (and thus 0)
     * or exactly the first non-zero matrix entry, an exception may be thrown.
     */
    value_type const& getValueLeft(index_type rowIndex, index_type colIndex) const {
        return matrix(rowIndex, colIndex);
    }

    /*!
     * Returns a const reference to the value at the given row / column.
     * If the column is not to the right of the existing matrix entries (and thus 0)
     * or exactly the last non-zero matrix entry, an exception may be thrown.
     */
    value_type const& getValueRight(index_type rowIndex, index_type colIndex) const {
        return matrix(rowIndex, colIndex);
    }

    /*!
     * Returns an iterator to the first element in the given row,
     * optionally to the first element with column index greater than first_col.
     */
    const_iterator row_begin(index_type rowIndex, index_type first_col=0) const {
        return const_iterator(matrix, rowIndex, first_col);
    }

    /*!
     * Returns an iterator to the end (beyond the last element) in the given row.
     */
    const_iterator row_end(index_type rowIndex) const {
        return const_iterator(matrix, rowIndex);
    }

    /*!
     * Returns a column range (pair of iterators) for a subset of the columns
     * of the given row.
     * @param first_col the column index of the first column
     * @param last_col  the column index of the last column (inclusive)
     */
    column_range row_col_range(index_type rowIndex, index_type first_col, index_type last_col) const {
        return column_range(matrix, rowIndex, first_col, last_col);
    }

    /*!
     * Returns the number of rows of the matrix.
     *
     * @return The number of rows of the matrix.
     */
    index_type getRowCount() const {
        return matrix.getRowCount();
    }

    /*!
     * Returns the number of columns of the matrix.
     *
     * @return The number of columns of the matrix.
     */
    index_type getColumnCount() const {
        return matrix.getColumnCount();
    }

private:
    const Matrix& matrix;
};


}
}
