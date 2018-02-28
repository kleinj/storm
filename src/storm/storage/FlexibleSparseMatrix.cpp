#include "storm/storage/FlexibleSparseMatrix.h"

#include "storm/storage/SparseMatrix.h"
#include "storm/storage/BitVector.h"
#include "storm/adapters/RationalFunctionAdapter.h"

#include "storm/utility/macros.h"
#include "storm/utility/constants.h"
#include "storm/utility/NumericalTypesConverter.h"
#include "storm/exceptions/InvalidArgumentException.h"
#include "storm/exceptions/InvalidAccessException.h"
#include "storm/exceptions/NotImplementedException.h"

#include <algorithm>

namespace storm {
    namespace storage {
        template<typename ValueType>
        FlexibleSparseMatrix<ValueType>::FlexibleSparseMatrix() : FlexibleSparseMatrix(0) {
            // Intentionally left empty.
        }

        template<typename ValueType>
        FlexibleSparseMatrix<ValueType>::FlexibleSparseMatrix(index_type rows) : data(rows), columnCount(0), nonzeroEntryCount(0), trivialRowGrouping(true), zero_val(storm::utility::zero<value_type>()) {
            // Intentionally left empty.
        }
        
        template<typename ValueType>
        template<typename SourceValueType>
        FlexibleSparseMatrix<ValueType>::FlexibleSparseMatrix(storm::storage::SparseMatrix<SourceValueType> const& matrix, bool setAllValuesToOne, bool revertEquationSystem) : data(matrix.getRowCount()), columnCount(matrix.getColumnCount()), nonzeroEntryCount(matrix.getNonzeroEntryCount()), trivialRowGrouping(matrix.hasTrivialRowGrouping()), zero_val(storm::utility::zero<value_type>()) {
            STORM_LOG_THROW(!revertEquationSystem || trivialRowGrouping, storm::exceptions::InvalidArgumentException, "Illegal option for creating flexible matrix.");
            
            storm::utility::NumericalTypesConverter<SourceValueType,ValueType> converter;

            if (!trivialRowGrouping) {
                rowGroupIndices = matrix.getRowGroupIndices();
            }
            for (index_type rowIndex = 0; rowIndex < matrix.getRowCount(); ++rowIndex) {
                auto const& row = matrix.getRow(rowIndex);
                reserveInRow(rowIndex, row.getNumberOfEntries());
                for (auto const& element : row) {
                    // If the probability is zero, we skip this entry.
                    if (storm::utility::isZero(element.getValue())) {
                        if (revertEquationSystem && rowIndex == element.getColumn()) {
                            getRow(rowIndex).emplace_back(element.getColumn(), storm::utility::one<ValueType>());
                        } else {
                            continue;
                        }
                    }
                    if (setAllValuesToOne) {
                        if (revertEquationSystem && element.getColumn() == rowIndex && storm::utility::isOne(element.getValue())) {
                            continue;
                        } else {
                            getRow(rowIndex).emplace_back(element.getColumn(), storm::utility::one<ValueType>());
                        }
                    } else {
                        if (revertEquationSystem) {
                            if (element.getColumn() == rowIndex) {
                                if (storm::utility::isOne(element.getValue())) {
                                    continue;
                                }
                                getRow(rowIndex).emplace_back(element.getColumn(), storm::utility::one<ValueType>() - converter.convert(element.getValue()));
                            } else {
                                getRow(rowIndex).emplace_back(element.getColumn(), -converter.convert(element.getValue()));
                            }
                        } else {
                            getRow(rowIndex).emplace_back(element.getColumn(), converter.convert(element.getValue()));
                        }
                    }
                }
            }
        }

        template<typename TargetValueType>
        template<typename SourceValueType>
        FlexibleSparseMatrix<TargetValueType>
        FlexibleSparseMatrix<TargetValueType>::augmentedFromSparseMatrix(storm::storage::SparseMatrix<SourceValueType> const& matrix, std::vector<SourceValueType> const &b) {
            STORM_LOG_THROW(matrix.hasTrivialRowGrouping(), storm::exceptions::InvalidArgumentException, "Illegal option for creating flexible matrix.");

            FlexibleSparseMatrix<TargetValueType> result(matrix.getRowCount());
            result.columnCount = matrix.getColumnCount() + 1;
            result.trivialRowGrouping = matrix.hasTrivialRowGrouping();

            index_type bColIndex = result.columnCount - 1;
            index_type nnz = matrix.getNonzeroEntryCount();

            storm::utility::NumericalTypesConverter<SourceValueType,TargetValueType> converter;

            for (index_type rowIndex = 0; rowIndex < matrix.getRowCount(); ++rowIndex) {
                typename storm::storage::SparseMatrix<SourceValueType>::const_rows row = matrix.getRow(rowIndex);
                result.reserveInRow(rowIndex, row.getNumberOfEntries() + 1);
                for (auto const& element : row) {
                    // If the probability is zero, we skip this entry.
                    if (!storm::utility::isZero(element.getValue())) {
                        result.getRow(rowIndex).emplace_back(element.getColumn(), converter.convert(element.getValue()));
                    }
                }
                if (!storm::utility::isZero(b[rowIndex])) {
                    result.getRow(rowIndex).emplace_back(bColIndex, converter.convert(b[rowIndex]));
                    ++nnz;
                }
            }

            result.nonzeroEntryCount = nnz;

            return result;
        }

        template<typename ValueType>
        void FlexibleSparseMatrix<ValueType>::reserveInRow(index_type row, index_type numberOfElements) {
            this->data[row].reserve(numberOfElements);
        }
        
        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::row_type& FlexibleSparseMatrix<ValueType>::getRow(index_type index) {
            return this->data[index];
        }
        
        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::row_type const& FlexibleSparseMatrix<ValueType>::getRow(index_type index) const {
            return this->data[index];
        }
        
        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::row_type& FlexibleSparseMatrix<ValueType>::getRow(index_type rowGroup, index_type offset) {
            STORM_LOG_ASSERT(rowGroup < this->getRowGroupCount(), "Invalid rowGroup.");
            STORM_LOG_ASSERT(offset < this->getRowGroupSize(rowGroup), "Invalid offset.");
            return getRow(rowGroupIndices[rowGroup] + offset);
        }
        
        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::row_type const& FlexibleSparseMatrix<ValueType>::getRow(index_type rowGroup, index_type offset) const {
            STORM_LOG_ASSERT(rowGroup < this->getRowGroupCount(), "Invalid rowGroup.");
            STORM_LOG_ASSERT(offset < this->getRowGroupSize(rowGroup), "Invalid offset.");
            return getRow(rowGroupIndices[rowGroup] + offset);
        }

        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::value_type const& FlexibleSparseMatrix<ValueType>::getValue(index_type rowIndex, index_type colIndex) const {
            const row_type& row = getRow(rowIndex);
            for (const_iterator it = row.begin(); it != row.end(); ++it) {
                index_type c = it->getColumn();
                if (c == colIndex) {
                    return it->getValue();
                }
                if (c > colIndex) {
                    return zero_val;
                }
            }
            return zero_val;
        }

        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::value_type const& FlexibleSparseMatrix<ValueType>::getValueFast(index_type rowIndex, index_type colIndex) const {
            const row_type& row = getRow(rowIndex);
            if (row.empty()) {
                return zero_val;
            }

            const entry_type& first_entry = row.front();
            index_type first_col = first_entry.getColumn();
            if (first_col == colIndex) {
                return first_entry.getValue();
            } else if (first_col > colIndex) {
                return zero_val;
            }

            const entry_type& last_entry = row.back();
            index_type last_col = last_entry.getColumn();
            if (last_col == colIndex) {
                return last_entry.getValue();
            } else if (last_col < colIndex) {
                return zero_val;
            }

            STORM_LOG_THROW(true, storm::exceptions::InvalidAccessException, "getValueFast failed, trying to access column " << colIndex <<", entries are between " << first_col << " and " << last_col);
        }

        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::value_type const& FlexibleSparseMatrix<ValueType>::getValueLeft(index_type rowIndex, index_type colIndex) const {
            const row_type& row = getRow(rowIndex);
            if (row.empty()) {
                return zero_val;
            }

            const entry_type& entry = row.front();
            index_type first_col = entry.getColumn();
            if (first_col == colIndex) {
                return entry.getValue();
            } else if (first_col > colIndex) {
                return zero_val;
            }

            STORM_LOG_THROW(true, storm::exceptions::InvalidAccessException, "getValueLeft failed, trying to access column " << rowIndex << ", first entry is at " << first_col);
        }

        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::value_type const& FlexibleSparseMatrix<ValueType>::getValueRight(index_type rowIndex, index_type colIndex) const {
            const row_type& row = getRow(rowIndex);
            if (row.empty()) {
                return zero_val;
            }

            const entry_type& entry = row.back();
            index_type last_col = entry.getColumn();
            if (last_col == colIndex) {
                return entry.getValue();
            } else if (last_col < colIndex) {
                return zero_val;
            }

            STORM_LOG_THROW(true, storm::exceptions::InvalidAccessException, "getValueRight failed, trying to access column " << colIndex << ", last entry is at " << last_col);
        }

        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::const_iterator FlexibleSparseMatrix<ValueType>::row_begin(index_type rowIndex, index_type first_col) const
        {
            const row_type& row = getRow(rowIndex);
            for (const_iterator it = row.begin(); it != row.end(); ++it) {
                if (it->getColumn() >= first_col) {
                    return it;
                }
            }
            return row.end();
        }

        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::const_iterator FlexibleSparseMatrix<ValueType>::row_end(index_type rowIndex) const {
            const row_type& row = getRow(rowIndex);
            return row.end();
        }

        template <typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::column_range
        FlexibleSparseMatrix<ValueType>::row_col_range(index_type rowIndex, index_type first_col, index_type last_col) const {
            return column_range(getRow(rowIndex), first_col, last_col);
        }

        template <typename ValueType>
        FlexibleSparseMatrix<ValueType>::column_range::column_range(const row_type& row, index_type first_col, index_type last_col)
            : row(row), first_col(first_col), last_col(last_col) {
        }

        template <typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::column_range::const_iterator
        FlexibleSparseMatrix<ValueType>::column_range::begin() const {
            return const_iterator(*this);
        }

        template <typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::column_range::const_iterator
        FlexibleSparseMatrix<ValueType>::column_range::end() const {
            return const_iterator(*this, true);
        }

        template <typename ValueType>
        FlexibleSparseMatrix<ValueType>::column_range::const_iterator::const_iterator(const column_range& range, bool done)
            : range(range), done(done) {
            cur = range.row.begin();
            end = range.row.end();

            if (range.first_col > range.last_col) {
                this->done = true;
                return;
            }

            while (!check_done()) {
                if (cur->getColumn() >= range.first_col) {
                    return;
                }
                ++cur;
            }

            this->done = true;
        }

        template <typename ValueType>
        bool FlexibleSparseMatrix<ValueType>::column_range::const_iterator::check_done() {
            if (done || cur == end || cur->getColumn() > range.last_col) {
                return true;
            }
            return false;
        }


        template <typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::column_range::const_iterator&
        FlexibleSparseMatrix<ValueType>::column_range::const_iterator::operator++() {
            // prefix operator
            ++cur;
            if (check_done()) {
                done = true;
            }
            return *this;
        }

        template <typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::column_range::const_iterator
        FlexibleSparseMatrix<ValueType>::column_range::const_iterator::operator++(int) {
            // postfix operator
            const_iterator i = *this;
            ++cur;
            if (check_done()) {
                done = true;
            }
            return i;
        }

        template <typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::column_range::const_iterator::const_reference
        FlexibleSparseMatrix<ValueType>::column_range::const_iterator::operator*() {
            STORM_LOG_ASSERT(!done, "Accessing invalid iterator");
            return *cur;
        }

        template <typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::column_range::const_iterator::const_pointer
        FlexibleSparseMatrix<ValueType>::column_range::const_iterator::operator->() {
            STORM_LOG_ASSERT(!done, "Accessing invalid iterator");
            return &(*cur);
        }

        template <typename ValueType>
        bool FlexibleSparseMatrix<ValueType>::column_range::const_iterator::operator==(const const_iterator& other) {
            if (done || other.done) {
                return done == other.done && &range == &other.range;
            } else {
                return cur == other.cur;
            }
        }

        template <typename ValueType>
        bool FlexibleSparseMatrix<ValueType>::column_range::const_iterator::operator!=(const const_iterator& other) {
            return !(*this == other);
        }

        template <typename ValueType>
        const ValueType& FlexibleSparseMatrix<ValueType>::getZeroValue() const {
            return zero_val;
        }

        template<typename ValueType>
        std::vector<typename FlexibleSparseMatrix<ValueType>::index_type> const& FlexibleSparseMatrix<ValueType>::getRowGroupIndices() const {
            return rowGroupIndices;
        }
        
        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::index_type FlexibleSparseMatrix<ValueType>::getRowCount() const {
            return this->data.size();
        }
        
        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::index_type FlexibleSparseMatrix<ValueType>::getColumnCount() const {
            return columnCount;
        }
        
        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::index_type FlexibleSparseMatrix<ValueType>::getNonzeroEntryCount() const {
            return nonzeroEntryCount;
        }
        
        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::index_type FlexibleSparseMatrix<ValueType>::getRowGroupCount() const {
            return rowGroupIndices.size() - 1;
        }
        
        template<typename ValueType>
        typename FlexibleSparseMatrix<ValueType>::index_type FlexibleSparseMatrix<ValueType>::getRowGroupSize(index_type group) const {
            return rowGroupIndices[group + 1] - rowGroupIndices[group];
        }
        
        template<typename ValueType>
        ValueType FlexibleSparseMatrix<ValueType>::getRowSum(index_type row) const {
            ValueType sum = storm::utility::zero<ValueType>();
            for (auto const& element : getRow(row)) {
                sum += element.getValue();
            }
            return sum;
        }

        template<typename ValueType>
        void FlexibleSparseMatrix<ValueType>::updateDimensions() {
            this->nonzeroEntryCount = 0;
            this->columnCount = 0;
            for (auto const& row : this->data) {
                for (auto const& element : row) {
                    STORM_LOG_ASSERT(!storm::utility::isZero(element.getValue()), "Entry is 0.");
                    ++this->nonzeroEntryCount;
                    this->columnCount = std::max(element.getColumn() + 1, this->columnCount);
                }
            }
        }

        template<typename ValueType>
        bool FlexibleSparseMatrix<ValueType>::empty() const {
            for (auto const& row : this->data) {
                if (!row.empty()) {
                    return false;
                }
            }
            return true;
        }
        
        template<typename ValueType>
        bool FlexibleSparseMatrix<ValueType>::hasTrivialRowGrouping() const {
            return trivialRowGrouping;
        }
        
        template<typename ValueType>
        void FlexibleSparseMatrix<ValueType>::filterEntries(storm::storage::BitVector const& rowConstraint, storm::storage::BitVector const& columnConstraint) {
            for (uint_fast64_t rowIndex = 0; rowIndex < this->data.size(); ++rowIndex) {
                auto& row = this->data[rowIndex];
                if (!rowConstraint.get(rowIndex)) {
                    row.clear();
                    row.shrink_to_fit();
                    continue;
                }
                row_type newRow;
                for (auto const& element : row) {
                    if (columnConstraint.get(element.getColumn())) {
                        newRow.push_back(element);
                    }
                }
                row = std::move(newRow);
            }
        }

        template<typename ValueType>
        void FlexibleSparseMatrix<ValueType>::permute(const Permutation& permutation) {
            STORM_LOG_THROW(trivialRowGrouping, storm::exceptions::NotImplementedException, "FlexibleSparseMatrix::permute in the presence of row groupings not implemented");

            // change column indizes
            for (auto &row : this->data) {
                for (auto &entry : row) {
                    entry.setColumn(permutation.get(entry.getColumn()));
                }
                // sort to obtain monotonic order again
                std::sort(row.begin(), row.end(), typename entry_type::column_compare());
            }

            // now, we have just permute the rows
            permutation.permute(this->data);
        }

        template<typename ValueType>
        template<typename T>
        storm::storage::SparseMatrix<T> FlexibleSparseMatrix<ValueType>::createSparseMatrix() {
            uint_fast64_t numEntries = 0;
            for (auto const& row : this->data) {
                numEntries += row.size();
            }
            
            storm::utility::NumericalTypesConverter<ValueType,T> converter;

            storm::storage::SparseMatrixBuilder<T> matrixBuilder(getRowCount(), getColumnCount(), numEntries, hasTrivialRowGrouping(), hasTrivialRowGrouping() ? 0 : getRowGroupCount());
            uint_fast64_t currRowIndex = 0;
            auto rowGroupIndexIt = getRowGroupIndices().begin();
            for (auto const& row : this->data) {
                if(!hasTrivialRowGrouping()) {
                    while (currRowIndex >= *rowGroupIndexIt) {
                        matrixBuilder.newRowGroup(currRowIndex);
                        ++rowGroupIndexIt;
                    }
                }
                for (auto const& entry : row) {
                    matrixBuilder.addNextValue(currRowIndex, entry.getColumn(), converter.convert(entry.getValue()));
                }
                ++currRowIndex;
            }
            // The matrix might end with one or more empty row groups
            if(!hasTrivialRowGrouping()) {
                while (currRowIndex >= *rowGroupIndexIt) {
                    matrixBuilder.newRowGroup(currRowIndex);
                    ++rowGroupIndexIt;
                }
            }
            return matrixBuilder.build();
        }
        
        template<typename ValueType>
        template<typename T>
        storm::storage::SparseMatrix<T> FlexibleSparseMatrix<ValueType>::createSparseMatrix(storm::storage::BitVector const& rowConstraint, storm::storage::BitVector const& columnConstraint) {
            uint_fast64_t numEntries = 0;
            for (auto const& rowIndex : rowConstraint) {
                auto const& row = data[rowIndex];
                for(auto const& entry : row) {
                    if (columnConstraint.get(entry.getColumn())) {
                        ++numEntries;
                    }
                }
            }
            uint_fast64_t numRowGroups = 0;
            if (!hasTrivialRowGrouping()) {
                auto lastRowGroupIndexIt = getRowGroupIndices().end() - 1;
                auto rowGroupIndexIt = getRowGroupIndices().begin();
                while (rowGroupIndexIt != lastRowGroupIndexIt) {
                    // Check whether the rowGroup will be nonempty
                    if(rowConstraint.getNextSetIndex(*rowGroupIndexIt) < *(++rowGroupIndexIt)) {
                        ++numRowGroups;
                    }
                }
            }
            
            std::vector<uint_fast64_t> oldToNewColumnIndexMapping(getColumnCount(), getColumnCount());
            uint_fast64_t newColumnIndex = 0;
            for (auto const& oldColumnIndex : columnConstraint) {
                oldToNewColumnIndexMapping[oldColumnIndex] = newColumnIndex++;
            }

            storm::utility::NumericalTypesConverter<ValueType,T> converter;
            
            storm::storage::SparseMatrixBuilder<T> matrixBuilder(rowConstraint.getNumberOfSetBits(), newColumnIndex, numEntries, true, !hasTrivialRowGrouping(), numRowGroups);
            uint_fast64_t currRowIndex = 0;
            auto rowGroupIndexIt = getRowGroupIndices().begin();
            for (auto const& oldRowIndex : rowConstraint) {
                if(!hasTrivialRowGrouping() && oldRowIndex >= *rowGroupIndexIt) {
                    matrixBuilder.newRowGroup(currRowIndex);
                    // Skip empty row groups
                    do {
                        ++rowGroupIndexIt;
                    } while (oldRowIndex >= *rowGroupIndexIt);
                }
                auto const& row = data[oldRowIndex];
                for (auto const& entry : row) {
                    if(columnConstraint.get(entry.getColumn())) {
                        matrixBuilder.addNextValue(currRowIndex, oldToNewColumnIndexMapping[entry.getColumn()], converter.convert(entry.getValue()));
                    }
                }
                ++currRowIndex;
            }
            return matrixBuilder.build();
        }

        template<typename ValueType>
        bool FlexibleSparseMatrix<ValueType>::rowHasDiagonalElement(storm::storage::sparse::state_type state) {
            for (auto const& entry : this->getRow(state)) {
                if (entry.getColumn() < state) {
                    continue;
                } else if (entry.getColumn() > state) {
                    return false;
                } else if (entry.getColumn() == state) {
                    return true;
                }
            }
            return false;
        }

        template<typename ValueType>
        void FlexibleSparseMatrix<ValueType>::replaceRow(index_type rowIndex, row_type&& new_row) {
            if (!new_row.empty() && new_row.back().getColumn() > columnCount+1) {
                columnCount = new_row.back().getColumn() + 1;
            }
            data[rowIndex] = new_row;
        }

        template<typename ValueType>
        std::ostream& FlexibleSparseMatrix<ValueType>::printRow(std::ostream& out, index_type const& rowIndex) const {
            index_type columnIndex = 0;
            row_type row = this->getRow(rowIndex);
            for (index_type column = 0; column < this->getColumnCount(); ++column) {
                if (columnIndex < row.size() && row[columnIndex].getColumn() == column) {
                    // Insert entry
                    out << row[columnIndex].getValue() << "\t";
                    ++columnIndex;
                } else {
                    // Insert zero
                    out << "0\t";
                }
            }
            return out;
        }
        

        template<typename ValueType>
        std::ostream& operator<<(std::ostream& out, FlexibleSparseMatrix<ValueType> const& matrix) {
            typedef typename FlexibleSparseMatrix<ValueType>::index_type FlexibleIndex;
            
            // Print column numbers in header.
            out << "\t\t";
            for (FlexibleIndex i = 0; i < matrix.getColumnCount(); ++i) {
                out << i << "\t";
            }
            out << std::endl;
            
            if (!matrix.hasTrivialRowGrouping()) {
                // Iterate over all row groups
                FlexibleIndex rowGroupCount = matrix.getRowGroupCount();
                for (FlexibleIndex rowGroup = 0; rowGroup < rowGroupCount; ++rowGroup) {
                    out << "\t---- group " << rowGroup << "/" << (rowGroupCount - 1) << " ---- " << std::endl;
                    FlexibleIndex endRow = matrix.rowGroupIndices[rowGroup + 1];
                    // Iterate over all rows.
                    for (FlexibleIndex row = matrix.rowGroupIndices[rowGroup]; row < endRow; ++row) {
                        // Print the actual row.
                        out << rowGroup << "\t(\t";
                        matrix.printRow(out, row);
                        out << "\t)\t" << rowGroup << std::endl;
                    }
                }

            } else {
                // Iterate over all rows
                for (FlexibleIndex row = 0; row < matrix.getRowCount(); ++row) {
                    // Print the actual row.
                    out << row << "\t(\t";
                    matrix.printRow(out, row);
                    out << "\t)\t" << row << std::endl;
                }
            }
            
            // Print column numbers in footer.
            out << "\t\t";
            for (FlexibleIndex i = 0; i < matrix.getColumnCount(); ++i) {
                out << i << "\t";
            }
            out << std::endl;
            return out;
        }

#define INSTANTIATE(T) \
        template class FlexibleSparseMatrix<T>; \
        template std::ostream& operator<<(std::ostream& out, FlexibleSparseMatrix<T> const& matrix);

#define INSTANTIATE_CONSTRUCTION(FROM,TO) \
        template FlexibleSparseMatrix<TO>::FlexibleSparseMatrix(storm::storage::SparseMatrix<FROM> const& matrix, bool setAllValuesToOne, bool revertEquationSystem); \
        template FlexibleSparseMatrix<TO> FlexibleSparseMatrix<TO>::augmentedFromSparseMatrix(storm::storage::SparseMatrix<FROM> const& matrix, std::vector<FROM> const &b);

#define INSTANTIATE_SPARSE_MATRIX_CONSTRUCTION(T) \
        template SparseMatrix<T> FlexibleSparseMatrix<T>::createSparseMatrix(); \
        template SparseMatrix<T> FlexibleSparseMatrix<T>::createSparseMatrix(BitVector const& rowConstraint, BitVector const& columnConstraint);

#define INSTANTIATE_ALL(T) \
        INSTANTIATE(T) \
        INSTANTIATE_CONSTRUCTION(T,T) \
        INSTANTIATE_SPARSE_MATRIX_CONSTRUCTION(T)

        // Explicitly instantiate the matrix.

        // double
        INSTANTIATE_ALL(double)

        // float
        INSTANTIATE_ALL(float)

        // int
        INSTANTIATE_ALL(int)

        // state_type
        INSTANTIATE_ALL(storm::storage::sparse::state_type)

#ifdef STORM_HAVE_CARL
#if defined(STORM_HAVE_CLN)
        INSTANTIATE_ALL(storm::ClnRationalNumber)
#endif

#if defined(STORM_HAVE_GMP)
        INSTANTIATE_ALL(GmpRationalNumber)
#endif

        INSTANTIATE_ALL(storm::RationalFunction)

        INSTANTIATE(storm::PlainRationalFunction)
        INSTANTIATE_CONSTRUCTION(storm::RationalFunction,storm::PlainRationalFunction)

        INSTANTIATE(storm::PlainRationalFunctionNoSimplify)
        INSTANTIATE_CONSTRUCTION(storm::RationalFunction,storm::PlainRationalFunctionNoSimplify)

        INSTANTIATE(storm::FactorizedRationalFunctionNoSimplify)
        INSTANTIATE_CONSTRUCTION(storm::RationalFunction,storm::FactorizedRationalFunctionNoSimplify)

        INSTANTIATE(storm::RawPolynomial)
        INSTANTIATE_CONSTRUCTION(storm::RationalFunction,storm::RawPolynomial)

        // Intervals
        INSTANTIATE_ALL(storm::Interval)

#endif

    } // namespace storage
} // namespace storm
