#pragma once

#include <utility>

namespace storm {
    namespace storage {

        template<typename IndexType, typename ValueType>
        class MatrixEntry {
        public:
            typedef IndexType index_type;
            typedef ValueType value_type;

            /*!
             * Constructs a matrix entry with the given column and value.
             *
             * @param column The column of the matrix entry.
             * @param value The value of the matrix entry.
             */
            MatrixEntry(index_type column, value_type value) : entry(column, value) {
                // Intentionally left empty.
            }

            /*!
             * Move-constructs the matrix entry fro the given column-value pair.
             *
             * @param pair The column-value pair from which to move-construct the matrix entry.
             */
            MatrixEntry(std::pair<index_type, value_type>&& pair) : entry(std::move(pair)) {
                // Intentionally left empty.
            }


            MatrixEntry() = default;
            MatrixEntry(MatrixEntry const& other) = default;
            MatrixEntry& operator=(MatrixEntry const& other) = default;
#ifndef WINDOWS
            MatrixEntry(MatrixEntry&& other) = default;
            MatrixEntry& operator=(MatrixEntry&& other) = default;
#endif

            /*!
             * Retrieves the column of the matrix entry.
             *
             * @return The column of the matrix entry.
             */
            index_type const& getColumn() const {
                return this->entry.first;
            }

            /*!
             * Sets the column of the current entry.
             *
             * @param column The column to set for this entry.
             */
            void setColumn(index_type const& column) {
                this->entry.first = column;
            }


            /*!
             * Retrieves the value of the matrix entry.
             *
             * @return The value of the matrix entry.
             */
            value_type const& getValue() const {
                return this->entry.second;
            }


            /*!
             * Sets the value of the entry in the matrix.
             *
             * @param value The value that is to be set for this entry.
             */
            void setValue(value_type const& value) {
                this->entry.second = value;
            }


            /*!
             * Retrieves a pair of column and value that characterizes this entry.
             *
             * @return A column-value pair that characterizes this entry.
             */
            std::pair<index_type, value_type> const& getColumnValuePair() const {
                return this->entry;
            }

            /*!
             * Multiplies the entry with the given factor and returns the result.
             *
             * @param factor The factor with which to multiply the entry.
             */
            MatrixEntry operator*(value_type factor) const{
                return MatrixEntry(this->getColumn(), this->getValue() * factor);
            }

            bool operator==(MatrixEntry const& other) const {
                return this->entry.first == other.entry.first && this->entry.second == other.entry.second;
            }

            bool operator!=(MatrixEntry const& other) const {
                return !(*this == other);
            }

            struct column_compare {
                bool operator() (const MatrixEntry& a, const MatrixEntry& b) {
                    return a.getColumn() < b.getColumn();
                }
            };

            template<typename IndexTypePrime, typename ValueTypePrime>
            friend std::ostream& operator<<(std::ostream& out, MatrixEntry<IndexTypePrime, ValueTypePrime> const& entry);
        private:
            // The actual matrix entry.
            std::pair<index_type, value_type> entry;
        };

        template<typename IndexTypePrime, typename ValueTypePrime>
        inline std::ostream& operator<<(std::ostream& out, MatrixEntry<IndexTypePrime, ValueTypePrime> const& entry) {
            out << "(" << entry.getColumn() << ", " << entry.getValue() << ")";
            return out;
        }

        /*!
         * Computes the hash value of a matrix entry.
         */
        template<typename IndexType, typename ValueType>
        inline std::size_t hash_value(MatrixEntry<IndexType, ValueType> const& matrixEntry) {
            std::size_t seed = 0;
            boost::hash_combine(seed, matrixEntry.getColumn());
            boost::hash_combine(seed, matrixEntry.getValue());
            return seed;
        }

    }
}
