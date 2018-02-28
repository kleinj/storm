#pragma once

#include "storm/solver/LinearEquationSolver.h"
#include "storm/solver/SolverSelectionOptions.h"
#include "storm/utility/eigen.h"

namespace storm {
    namespace solver {
        
        /*!
         * A class that uses the Eigen library to implement the LinearEquationSolver interface.
         */
        template<typename ValueType>
        class EigenLinearEquationSolver : public LinearEquationSolver<ValueType> {
        public:
            EigenLinearEquationSolver();
            EigenLinearEquationSolver(storm::storage::SparseMatrix<ValueType> const& A);
            EigenLinearEquationSolver(storm::storage::SparseMatrix<ValueType>&& A);
            
            virtual void setMatrix(storm::storage::SparseMatrix<ValueType> const& A) override;
            virtual void setMatrix(storm::storage::SparseMatrix<ValueType>&& A) override;
            
            virtual void multiply(std::vector<ValueType>& x, std::vector<ValueType> const* b, std::vector<ValueType>& result) const override;

            virtual LinearEquationSolverProblemFormat getEquationProblemFormat(Environment const& env) const override;

        protected:
            virtual bool internalSolveEquations(Environment const& env, std::vector<ValueType>& x, std::vector<ValueType> const& b) const override;
            virtual bool internalSolve(Environment const& env, StormEigen::SparseMatrix<ValueType>& eigenMatrix, std::vector<ValueType>& x, std::vector<ValueType> const &b) const;

        private:
            EigenLinearEquationSolverMethod getMethod(Environment const& env, bool isExactMode) const;
            
            virtual uint64_t getMatrixRowCount() const override;
            virtual uint64_t getMatrixColumnCount() const override;
            
            // The (eigen) matrix associated with this equation solver.
            std::unique_ptr<StormEigen::SparseMatrix<ValueType>> eigenA;

            // If the solver takes posession of the matrix, we store the moved matrix in this member, so it gets deleted
            // when the solver is destructed.
            std::unique_ptr<storm::storage::SparseMatrix<ValueType>> localA;

            // A pointer to the original sparse matrix given to this solver. If the solver takes posession of the matrix
            // the pointer refers to localA.
            storm::storage::SparseMatrix<ValueType> const* A;
        };
        
        template<typename ValueType>
        class EigenLinearEquationSolverFactory : public LinearEquationSolverFactory<ValueType> {
        public:
            using LinearEquationSolverFactory<ValueType>::create;

            virtual std::unique_ptr<storm::solver::LinearEquationSolver<ValueType>> create(Environment const& env, LinearEquationSolverTask const& task = LinearEquationSolverTask::Unspecified) const override;

            virtual std::unique_ptr<LinearEquationSolverFactory<ValueType>> clone() const override;
        };
    }
}
