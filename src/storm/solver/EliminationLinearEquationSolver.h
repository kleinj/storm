#ifndef STORM_SOLVER_ELIMINATIONLINEAREQUATIONSOLVER_H_
#define STORM_SOLVER_ELIMINATIONLINEAREQUATIONSOLVER_H_

#include "storm/solver/LinearEquationSolver.h"

#include "storm/settings/modules/EliminationSettings.h"

namespace storm {
    namespace storage {
        // fwd
        template <typename ValueType> class SparseMatrix;
        template <typename ValueType> class FlexibleSparseMatrix;
        class HTMLMatrixExport;
    }

    namespace solver {


        /*!
         * A class that uses gaussian elimination to implement the LinearEquationSolver interface.
         */
        template<typename ValueType>
        class EliminationLinearEquationSolver : public LinearEquationSolver<ValueType> {
        public:
            EliminationLinearEquationSolver();
            EliminationLinearEquationSolver(storm::storage::SparseMatrix<ValueType> const& A);
            EliminationLinearEquationSolver(storm::storage::SparseMatrix<ValueType>&& A);
            
            virtual void setMatrix(storm::storage::SparseMatrix<ValueType> const& A) override;
            virtual void setMatrix(storm::storage::SparseMatrix<ValueType>&& A) override;
            
            virtual void multiply(std::vector<ValueType>& x, std::vector<ValueType> const* b, std::vector<ValueType>& result) const override;

            virtual LinearEquationSolverProblemFormat getEquationProblemFormat(Environment const& env) const override;

        protected:
            virtual bool internalSolveEquations(Environment const& env, std::vector<ValueType>& x, std::vector<ValueType> const& b) const override;
            bool internalSolve(Environment const& env, const storm::storage::SparseMatrix<ValueType>& transitionMatrix, const std::vector<ValueType>& b, std::vector<ValueType>& x, std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport) const;

        private:
            virtual uint64_t getMatrixRowCount() const override;
            virtual uint64_t getMatrixColumnCount() const override;
            
            // If the solver takes posession of the matrix, we store the moved matrix in this member, so it gets deleted
            // when the solver is destructed.
            std::unique_ptr<storm::storage::SparseMatrix<ValueType>> localA;

            // A pointer to the original sparse matrix given to this solver. If the solver takes posession of the matrix
            // the pointer refers to localA.
            storm::storage::SparseMatrix<ValueType> const* A;
        };
        
        template<typename ValueType>
        class EliminationLinearEquationSolverFactory : public LinearEquationSolverFactory<ValueType> {
        public:
            using LinearEquationSolverFactory<ValueType>::create;
            
            virtual std::unique_ptr<storm::solver::LinearEquationSolver<ValueType>> create(Environment const& env, LinearEquationSolverTask const& task = LinearEquationSolverTask::Unspecified) const override;
            
            virtual std::unique_ptr<LinearEquationSolverFactory<ValueType>> clone() const override;

        };
    }
}

#endif /* STORM_SOLVER_ELIMINATIONLINEAREQUATIONSOLVER_H_ */
