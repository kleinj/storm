#pragma once

#include "storm/solver/LinearEquationSolver.h"
#include "storm/storage/SCCInfo.h"
#include <memory>


namespace storm {
    namespace storage {
        // fwd
        class HTMLMatrixExport;
    }

    namespace solver {

        /*!
         * A class that uses gaussian elimination to implement the LinearEquationSolver interface.
         */
        template<typename ValueType>
        class GaussEliminationLinearEquationSolver : public LinearEquationSolver<ValueType> {
        public:
            GaussEliminationLinearEquationSolver();

            virtual void setMatrix(storm::storage::SparseMatrix<ValueType> const& A) override;
            virtual void setMatrix(storm::storage::SparseMatrix<ValueType>&& A) override;

            virtual LinearEquationSolverProblemFormat getEquationProblemFormat(Environment const& env) const override;

            virtual void multiply(std::vector<ValueType>& x, std::vector<ValueType> const* b, std::vector<ValueType>& result) const override;

        protected:
            virtual bool internalSolveEquations(Environment const& env, std::vector<ValueType>& x, std::vector<ValueType> const& b) const override;

        private:
            virtual uint64_t getMatrixRowCount() const override;
            virtual uint64_t getMatrixColumnCount() const override;

            void assignAndConvertIfNecessary(std::vector<ValueType>& a, std::vector<ValueType>& b) const;

            template <typename OtherValueType>
            void assignAndConvertIfNecessary(std::vector<ValueType>& a, std::vector<OtherValueType>& b) const;

            // If the solver takes possession of the matrix, we store the moved matrix in this member, so it gets deleted
            // when the solver is destructed.
            std::unique_ptr<storm::storage::SparseMatrix<ValueType>> localA;

            // A pointer to the original sparse matrix given to this solver. If the solver takes possession of the matrix
            // the pointer refers to localA.
            storm::storage::SparseMatrix<ValueType> const* A;

            template <typename Worker>
            void solve(std::vector<ValueType>& x, std::vector<ValueType> const& b) const;

            template <typename Worker, typename Matrix>
            void solve(Matrix& augmentedMatrix, std::vector<ValueType>& x, std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport, bool debug) const;

            template <typename Worker, typename Matrix, typename OtherValueType>
            void solve(Matrix& augmentedMatrix, boost::optional<const storm::storage::SCCInfo&> sccInfo, std::vector<OtherValueType>& x, std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport, bool debug) const;

        };

        template<typename ValueType>
        class GaussEliminationLinearEquationSolverFactory : public LinearEquationSolverFactory<ValueType> {
        public:
            virtual std::unique_ptr<storm::solver::LinearEquationSolver<ValueType>> create(Environment const& env, LinearEquationSolverTask const& task = LinearEquationSolverTask::Unspecified) const override;

            virtual std::unique_ptr<LinearEquationSolverFactory<ValueType>> clone() const override;
        };
    }
}
