#ifndef STORM_STORAGE_FLEXIBLESPARSEMATRIX_H_
#define STORM_STORAGE_FLEXIBLESPARSEMATRIX_H_

#include <cstdint>
#include <vector>

#include "storm/storage/sparse/StateType.h"
#include "storm/storage/SparseMatrix.h"
#include "storm/storage/MatrixEntry.h"
#include "storm/storage/Permutation.h"

namespace storm {
    namespace storage {
        
        class BitVector;
        
        /*!
         * The flexible sparse matrix is used during state elimination.
         */
        template<typename ValueType>
        class FlexibleSparseMatrix {
        public:
            // TODO: make this class a bit more consistent with the big sparse matrix and improve it:
            // * add stuff like iterator, clearRow, multiplyRowWithScalar
            
            typedef uint_fast64_t index_type;
            typedef ValueType value_type;
            typedef storm::storage::MatrixEntry<index_type, value_type> entry_type;
            typedef std::vector<entry_type> row_type;
            typedef typename row_type::iterator iterator;
            typedef typename row_type::const_iterator const_iterator;

            class column_range {
            public:
                class const_iterator {
                public:
                    typedef entry_type value_type;
                    typedef entry_type& reference;
                    typedef const entry_type& const_reference;
                    typedef entry_type* pointer;
                    typedef const entry_type* const_pointer;
                    typedef int difference_type;
                    typedef std::forward_iterator_tag iterator_category;

                    const_iterator(const column_range& range, bool done=false);
                    const_iterator(const const_iterator&) = default;
                    const_iterator(const_iterator&&) = default;

                    const_iterator& operator++();
                    const_iterator operator++(int);
                    const_reference operator*();
                    const_pointer operator->();
                    bool operator==(const const_iterator& other);
                    bool operator!=(const const_iterator& other);

                private:
                    bool check_done();

                    typename FlexibleSparseMatrix::const_iterator cur;
                    typename FlexibleSparseMatrix::const_iterator end;
                    const column_range& range;
                    bool done;
                };

                column_range(const row_type& row, index_type first_col, index_type last_col);
                const_iterator begin() const;
                const_iterator end() const;

            private:
                const row_type& row;
                index_type first_col;
                index_type last_col;
            };

            /*!
             * Constructs an empty flexible sparse matrix.
             */
            FlexibleSparseMatrix();

            /*!
             * Constructs a flexible sparse matrix with rows many rows.
             * @param rows number of rows.
             */
            FlexibleSparseMatrix(index_type rows);

            /*! Default copy constructor */
            FlexibleSparseMatrix(const FlexibleSparseMatrix<ValueType>&) = default;
            /*! Default move constructor */
            FlexibleSparseMatrix(FlexibleSparseMatrix<ValueType>&&) = default;

            /*!
             * Constructs a flexible sparse matrix from a sparse matrix.
             * @param matrix Sparse matrix to construct from.
             * @param setAllValuesToOne If true, all set entries are set to one. Default is false.
             * @param revertEquationSystem If true, the matrix that will be created is the matrix (1-A), where A is the
             * provided matrix.
             */
            template<typename SourceValueType>
            FlexibleSparseMatrix(storm::storage::SparseMatrix<SourceValueType> const& matrix, bool setAllValuesToOne = false, bool revertEquationSystem = false);

            /*!
             * Static constructor of a flexible sparse matrix from a sparse matrix and a vector,
             * potentially converting between source and target types.
             * @param matrix Sparse matrix to construct from.
             * @param v      Vector for the last column (has to have matrix.rowCount() entries)
             */
            template<typename SourceValueType>
            static FlexibleSparseMatrix augmentedFromSparseMatrix(storm::storage::SparseMatrix<SourceValueType> const& matrix, const std::vector<SourceValueType>& b);

            /*!
             * Reserves space for elements in row.
             * @param row Row to reserve in.
             * @param numberOfElements Number of elements to reserve space for.
             */
            void reserveInRow(index_type row, index_type numberOfElements);

            /*!
             * Returns an object representing the given row.
             *
             * @param row The row to get.
             * @return An object representing the given row.
             */
            row_type& getRow(index_type);

            /*!
             * Returns an object representing the given row.
             *
             * @param row The row to get.
             * @return An object representing the given row.
             */
            row_type const& getRow(index_type) const;
            
            /*!
             * Returns an object representing the offset'th row in the rowgroup
             * @param rowGroup the row group
             * @param offset which row in the group
             * @return An object representing the given row.
             */
            row_type& getRow(index_type rowGroup, index_type offset);
            
            /*!
             * Returns an object representing the offset'th row in the rowgroup
             * @param rowGroup the row group
             * @param offset which row in the group
             * @return An object representing the given row.
             */
            row_type const& getRow(index_type rowGroup, index_type entryInGroup) const;
            
            /*!
             * Returns the grouping of rows of this matrix.
             *
             * @return The grouping of rows of this matrix.
             */
            std::vector<index_type> const& getRowGroupIndices() const;

            /*!
             * Returns a const reference to the value at the given row / column
             * (using linear search in the sparse entries).
             * If the value is zero (and thus not stored in the matrix),
             * returns a const ref to a dummy zero value.
             */
            value_type const& getValue(index_type rowIndex, index_type colIndex) const;

            /*!
             * Returns a const reference to the value at the given row / column.
             * If the column is not to the left or right of the existing matrix entries (and thus 0)
             * or exactly the first / last non-zero matrix entry, an exception may be thrown.
             */
            value_type const& getValueFast(index_type rowIndex, index_type colIndex) const;

            /*!
             * Returns a const reference to the value at the given row / column.
             * If the column is not to the left of the existing matrix entries (and thus 0)
             * or exactly the first non-zero matrix entry, an exception may be thrown.
             */
            value_type const& getValueLeft(index_type rowIndex, index_type colIndex) const;

            /*!
             * Returns a const reference to the value at the given row / column.
             * If the column is not to the right of the existing matrix entries (and thus 0)
             * or exactly the last non-zero matrix entry, an exception may be thrown.
             */
            value_type const& getValueRight(index_type rowIndex, index_type colIndex) const;

            /*!
             * Returns an iterator to the first element in the given row,
             * optionally to the first element with column index greater than first_col.
             */
            const_iterator row_begin(index_type rowIndex, index_type first_col=0) const;

            /*!
             * Returns an iterator to the end (beyond the last element) in the given row.
             */
            const_iterator row_end(index_type rowIndex) const;

            /*!
             * Returns a column range (pair of iterators) for a subset of the columns
             * of the given row.
             * @param first_col the column index of the first column
             * @param last_col  the column index of the last column (inclusive)
             */
            column_range row_col_range(index_type rowIndex, index_type first_col, index_type last_col) const;

            /*!
             * Returns a const reference to the zero value.
             */
            const value_type& getZeroValue() const;

            /*!
             * Returns the number of rows of the matrix.
             *
             * @return The number of rows of the matrix.
             */
            index_type getRowCount() const;

            /*!
             * Returns the number of columns of the matrix.
             *
             * @return The number of columns of the matrix.
             */
            index_type getColumnCount() const;

            /*!
             * Returns the cached number of nonzero entries in the matrix.
             *
             * @return The number of nonzero entries in the matrix.
             */
            index_type getNonzeroEntryCount() const;
            
            /*!
             * Returns the number of row groups in the matrix.
             *
             * @return The number of row groups in the matrix.
             */
            index_type getRowGroupCount() const;
            
            /*!
             * Returns the size of the given row group.
             *
             * @param group The group whose size to retrieve.
             * @return The number of rows that belong to the given row group.
             */
            index_type getRowGroupSize(index_type group) const;
            
            /*!
             * Computes the sum of the entries in a given row.
             *
             * @param row The row that is to be summed.
             * @return The sum of the selected row.
             */
            value_type getRowSum(index_type row) const;

            /*!
             * Recomputes the number of columns and the number of non-zero entries.
             */
            void updateDimensions();

            /*!
             * Checks if the matrix has no elements.
             * @return True, if the matrix is empty.
             */
            bool empty() const;
            
            /*!
             * Retrieves whether the matrix has a (possibly) trivial row grouping.
             *
             * @return True iff the matrix has a (possibly) trivial row grouping.
             */
            bool hasTrivialRowGrouping() const;

            /*!
             * Erases all entries whose row and column does not satisfy the given rowConstraint and the given columnConstraint
             *
             * @param rowConstraint A bit vector indicating which row entries to keep.
             * @param columnConstraint A bit vector indicating which column entries to keep.
             */
            void filterEntries(storm::storage::BitVector const& rowConstraint, storm::storage::BitVector const& columnConstraint);

            /*!
             * Permutes this flexible sparse matrix according to the permutation.
             */
            void permute(const Permutation& permutation);

            /*!
             * Creates a sparse matrix from the flexible sparse matrix.
             * @return The sparse matrix.
             */
            template <typename T = ValueType>
            storm::storage::SparseMatrix<T> createSparseMatrix();
                        
            /*!
             * Creates a sparse matrix from the flexible sparse matrix.
             * Only the selected rows and columns will be considered.
             * Empty rowGroups will be ignored
             *
             * @param rowConstraint A bit vector indicating which rows to keep.
             * @param columnConstraint A bit vector indicating which columns to keep.

             *
             * @return The sparse matrix.
             */
            template <typename T = ValueType>
            storm::storage::SparseMatrix<T> createSparseMatrix(storm::storage::BitVector const& rowConstraint, storm::storage::BitVector const& columnConstraint);

            /*!
             * Checks whether the given state has a self-loop with an arbitrary probability in the probability matrix.
             *
             * @param state The state for which to check whether it possesses a self-loop.
             * @return True iff the given state has a self-loop with an arbitrary probability in the probability matrix.
             */
            bool rowHasDiagonalElement(storm::storage::sparse::state_type state);

            /*!
             * Replace the row with the given index by the new row.
             */
            void replaceRow(index_type rowIndex, row_type&& new_row);

            /*!
             * Print row.
             * @param out Output stream.
             * @param rowIndex Index of row to print.
             * @return Output with printed row.
             */
            std::ostream& printRow(std::ostream& out, index_type const& rowIndex) const;
            
            template<typename TPrime>
            friend std::ostream& operator<<(std::ostream& out, FlexibleSparseMatrix<TPrime> const& matrix);

        private:
            std::vector<row_type> data;

            // The number of columns of the matrix.
            index_type columnCount;

            // The number of entries in the matrix.
            index_type nonzeroEntryCount;
            
            // A flag indicating whether the matrix has a trivial row grouping. Note that this may be true and yet
            // there may be row group indices, because they were requested from the outside.
            bool trivialRowGrouping;
            
            // A vector indicating the row groups of the matrix.
            std::vector<index_type> rowGroupIndices;

            // storage for a zero value (to return as a const reference for missing values)
            const value_type zero_val;
        };
    }
}

#endif /* STORM_STORAGE_FLEXIBLESPARSEMATRIX_H_ */
