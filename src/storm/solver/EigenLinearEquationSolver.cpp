#include "storm/solver/EigenLinearEquationSolver.h"

#include "storm/adapters/EigenAdapter.h"

#include "storm/environment/solver/EigenSolverEnvironment.h"
#include "storm/storage/StronglyConnectedComponentDecomposition.h"
#include "storm/storage/TopologicalMatrixReordering.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/GeneralSettings.h"
#include "storm/settings/modules/EigenEquationSolverSettings.h"
#include "storm/settings/modules/CoreSettings.h"

#include "storm/utility/NumericalTypes.h"
#include "storm/utility/Stopwatch.h"
#include "storm/utility/vector.h"
#include "storm/utility/macros.h"
#include "storm/utility/GCDLog.h"
#include "storm/exceptions/InvalidSettingsException.h"

namespace storm {
    namespace solver {

        template<typename ValueType>
        EigenLinearEquationSolver<ValueType>::EigenLinearEquationSolver() : A(nullptr) {
            // Intentionally left empty.
        }

        template<typename ValueType>
        EigenLinearEquationSolver<ValueType>::EigenLinearEquationSolver(storm::storage::SparseMatrix<ValueType> const& A) : A(nullptr) {
            this->setMatrix(A);
        }

        template<typename ValueType>
        EigenLinearEquationSolver<ValueType>::EigenLinearEquationSolver(storm::storage::SparseMatrix<ValueType>&& A) : A(nullptr) {
            this->setMatrix(std::move(A));
        }
        
        template<typename ValueType>
        void EigenLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType> const& A) {
            this->A = &A;
            localA.reset();
            eigenA = storm::adapters::EigenAdapter::toEigenSparseMatrix<ValueType>(A);
            this->clearCache();
        }
        
        template<typename ValueType>
        void EigenLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType>&& A) {
            localA = std::make_unique<storm::storage::SparseMatrix<ValueType>>(std::move(A));
            this->A = localA.get();
            eigenA = storm::adapters::EigenAdapter::toEigenSparseMatrix<ValueType>(*localA);
            this->clearCache();
        }
        
        template<typename ValueType>
        EigenLinearEquationSolverMethod EigenLinearEquationSolver<ValueType>::getMethod(Environment const& env, bool isExactMode) const {
            // Adjust the method if none was specified and we are using rational numbers.
            auto method = env.solver().eigen().getMethod();
            
            if (isExactMode && method != EigenLinearEquationSolverMethod::SparseLU) {
                if (env.solver().eigen().isMethodSetFromDefault()) {
                    STORM_LOG_INFO("Selecting 'SparseLU' as the solution technique to guarantee exact results.");
                } else {
                    STORM_LOG_WARN("The selected solution method does not guarantee exact results. Falling back to SparseLU.");
                }
                method = EigenLinearEquationSolverMethod::SparseLU;
            } else if (env.solver().isForceSoundness() && method != EigenLinearEquationSolverMethod::SparseLU) {
                if (env.solver().eigen().isMethodSetFromDefault()) {
                    STORM_LOG_INFO("Selecting 'SparseLU' as the solution technique to guarantee sound results. If you want to override this, please explicitly specify a different method.");
                    method = EigenLinearEquationSolverMethod::SparseLU;
                } else {
                    STORM_LOG_WARN("The selected solution method does not guarantee sound results.");
                }
            }
            return method;
        }

        template<typename ValueType>
        bool EigenLinearEquationSolver<ValueType>::internalSolveEquations(Environment const& env, std::vector<ValueType>& x, std::vector<ValueType> const& b) const {
            bool topologicalOrdering = env.solver().eigen().isTopologicalOrderingSet();

            if (topologicalOrdering) {
                STORM_LOG_INFO("Obtaining SCC decomposition and topological ordering...");

                storm::utility::Stopwatch timerSCC(true);
                typedef storm::storage::StronglyConnectedComponentDecomposition<ValueType> SCCDecomposition;
                SCCDecomposition sccs(*A);
                auto sorted = storm::storage::TopologicalMatrixReordering::topologicallySorted(*A, sccs);

                auto& sortedMatrix = std::get<0>(sorted);
                const auto& sccInfo = *std::get<1>(sorted).sccInfo;
                const auto& permutation = *std::get<1>(sorted).permutation;

                std::vector<ValueType> b_(b);
                permutation.permute(b_);

                timerSCC.stop();
                STORM_LOG_INFO("SCC decomposition and topological ordering: " << timerSCC);
                STORM_LOG_INFO(sccInfo.getStatistics());

                auto sortedEigenA = storm::adapters::EigenAdapter::toEigenSparseMatrix<ValueType>(sortedMatrix);

                bool rv = internalSolve(env, *sortedEigenA, x, b_);

                STORM_LOG_INFO("Reorder solution vector to correspond to model ordering...");
                // un-permutate solution vector
                permutation.permuteInv(x);

                return rv;
            } else {
                return internalSolve(env, *eigenA, x, b);
            }
        }

        #ifdef STORM_HAVE_CARL
        // Specialization for storm::RationalNumber
        template<>
        bool EigenLinearEquationSolver<storm::RationalNumber>::internalSolve(Environment const& env, StormEigen::SparseMatrix<storm::RationalNumber>& eigenMatrix, std::vector<storm::RationalNumber>& x, std::vector<storm::RationalNumber> const& b) const {
            auto solutionMethod = getMethod(env, true);
            STORM_LOG_WARN_COND(solutionMethod == EigenLinearEquationSolverMethod::SparseLU, "Switching method to SparseLU.");
            STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with with rational numbers using LU factorization (Eigen library).");

#ifdef STORM_HAVE_GCD_LOG
            storm::utility::GCDLog::enable();
	    STORM_LOG_INFO("Enabling GCD logging...\n");
#endif
	    
            // Map the input vectors to Eigen's format.
            auto eigenX = StormEigen::Matrix<storm::RationalNumber, StormEigen::Dynamic, 1>::Map(x.data(), x.size());
            auto eigenB = StormEigen::Matrix<storm::RationalNumber, StormEigen::Dynamic, 1>::Map(b.data(), b.size());
            
            StormEigen::SparseLU<StormEigen::SparseMatrix<storm::RationalNumber>, StormEigen::COLAMDOrdering<int>> solver;
            solver.compute(eigenMatrix);
            solver._solve_impl(eigenB, eigenX);

#ifdef STORM_HAVE_GCD_LOG
            storm::utility::GCDLog::disable();
	    STORM_LOG_INFO("Disabling GCD logging...\n");
#endif
	    
            return solver.info() == StormEigen::ComputationInfo::Success;
        }
        
        // Specialization for storm::RationalFunction
        template<>
        bool EigenLinearEquationSolver<storm::RationalFunction>::internalSolve(Environment const& env, StormEigen::SparseMatrix<storm::RationalFunction>& eigenMatrix, std::vector<storm::RationalFunction>& x, std::vector<storm::RationalFunction> const& b) const {
            auto solutionMethod = getMethod(env, true);
            STORM_LOG_WARN_COND(solutionMethod == EigenLinearEquationSolverMethod::SparseLU, "Switching method to SparseLU.");
            STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with rational functions using LU factorization (Eigen library).");

#ifdef STORM_HAVE_GCD_LOG
            storm::utility::GCDLog::enable();
	    STORM_LOG_INFO("Enabling GCD logging...\n");
#endif	    
            // Map the input vectors to Eigen's format.
            auto eigenX = StormEigen::Matrix<storm::RationalFunction, StormEigen::Dynamic, 1>::Map(x.data(), x.size());
            auto eigenB = StormEigen::Matrix<storm::RationalFunction, StormEigen::Dynamic, 1>::Map(b.data(), b.size());

            storm::utility::Stopwatch timer(true);

            StormEigen::SparseLU<StormEigen::SparseMatrix<storm::RationalFunction>, StormEigen::COLAMDOrdering<int>> solver;
            solver.compute(eigenMatrix);
            solver._solve_impl(eigenB, eigenX);

#ifdef STORM_HAVE_GCD_LOG
            storm::utility::GCDLog::disable();
	    STORM_LOG_INFO("Disabling GCD logging...\n");
#endif

            timer.stop();
            STORM_LOG_INFO("Time for parametric computation: " << timer);

            if (storm::settings::getModule<storm::settings::modules::CoreSettings>().isResultStatsSet()) {
                std::cout << storm::utility::NumericalTypes::getStats(x);
            }

            return solver.info() == StormEigen::ComputationInfo::Success;
        }
#endif

        
        template<typename ValueType>
        bool EigenLinearEquationSolver<ValueType>::internalSolve(Environment const& env, StormEigen::SparseMatrix<ValueType>& eigenMatrix, std::vector<ValueType>& x, std::vector<ValueType> const &b) const {
            // Map the input vectors to Eigen's format.
            auto eigenX = StormEigen::Matrix<ValueType, StormEigen::Dynamic, 1>::Map(x.data(), x.size());
            auto eigenB = StormEigen::Matrix<ValueType, StormEigen::Dynamic, 1>::Map(b.data(), b.size());

            auto solutionMethod = getMethod(env, false);
            if (solutionMethod == EigenLinearEquationSolverMethod::SparseLU) {
                STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with sparse LU factorization (Eigen library).");
                StormEigen::SparseLU<StormEigen::SparseMatrix<ValueType>, StormEigen::COLAMDOrdering<int>> solver;
                solver.compute(eigenMatrix);
                solver._solve_impl(eigenB, eigenX);
            } else {
                bool converged = false;
                uint64_t numberOfIterations = 0;
                uint64_t maxIter = env.solver().eigen().getMaximalNumberOfIterations();
                uint64_t restartThreshold = env.solver().eigen().getRestartThreshold();
                ValueType precision = storm::utility::convertNumber<ValueType>(env.solver().eigen().getPrecision());
                EigenLinearEquationSolverPreconditioner preconditioner = env.solver().eigen().getPreconditioner();
                if (solutionMethod == EigenLinearEquationSolverMethod::Bicgstab) {
                    if (preconditioner == EigenLinearEquationSolverPreconditioner::Ilu) {
                        STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with BiCGSTAB with Ilu preconditioner (Eigen library).");

                        StormEigen::BiCGSTAB<StormEigen::SparseMatrix<ValueType>, StormEigen::IncompleteLUT<ValueType>> solver;
                        solver.compute(eigenMatrix);
                        solver.setTolerance(precision);
                        solver.setMaxIterations(maxIter);
                        eigenX = solver.solveWithGuess(eigenB, eigenX);
                        converged = solver.info() == StormEigen::ComputationInfo::Success;
                        numberOfIterations = solver.iterations();
                    } else if (preconditioner == EigenLinearEquationSolverPreconditioner::Diagonal) {
                        STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with BiCGSTAB with Diagonal preconditioner (Eigen library).");

                        StormEigen::BiCGSTAB<StormEigen::SparseMatrix<ValueType>, StormEigen::DiagonalPreconditioner<ValueType>> solver;
                        solver.setTolerance(precision);
                        solver.setMaxIterations(maxIter);
                        solver.compute(eigenMatrix);
                        eigenX = solver.solveWithGuess(eigenB, eigenX);
                        converged = solver.info() == StormEigen::ComputationInfo::Success;
                        numberOfIterations = solver.iterations();
                    } else {
                        STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with BiCGSTAB with identity preconditioner (Eigen library).");

                        StormEigen::BiCGSTAB<StormEigen::SparseMatrix<ValueType>, StormEigen::IdentityPreconditioner> solver;
                        solver.setTolerance(precision);
                        solver.setMaxIterations(maxIter);
                        solver.compute(eigenMatrix);
                        eigenX = solver.solveWithGuess(eigenB, eigenX);
                        numberOfIterations = solver.iterations();
                        converged = solver.info() == StormEigen::ComputationInfo::Success;
                    }
                } else if (solutionMethod == EigenLinearEquationSolverMethod::DGmres) {
                    if (preconditioner == EigenLinearEquationSolverPreconditioner::Ilu) {
                        STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with DGMRES with Ilu preconditioner (Eigen library).");

                        StormEigen::DGMRES<StormEigen::SparseMatrix<ValueType>, StormEigen::IncompleteLUT<ValueType>> solver;
                        solver.setTolerance(precision);
                        solver.setMaxIterations(maxIter);
                        solver.set_restart(restartThreshold);
                        solver.compute(eigenMatrix);
                        eigenX = solver.solveWithGuess(eigenB, eigenX);
                        converged = solver.info() == StormEigen::ComputationInfo::Success;
                        numberOfIterations = solver.iterations();
                    } else if (preconditioner == EigenLinearEquationSolverPreconditioner::Diagonal) {
                        STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with DGMRES with Diagonal preconditioner (Eigen library).");

                        StormEigen::DGMRES<StormEigen::SparseMatrix<ValueType>, StormEigen::DiagonalPreconditioner<ValueType>> solver;
                        solver.setTolerance(precision);
                        solver.setMaxIterations(maxIter);
                        solver.set_restart(restartThreshold);
                        solver.compute(eigenMatrix);
                        eigenX = solver.solveWithGuess(eigenB, eigenX);
                        converged = solver.info() == StormEigen::ComputationInfo::Success;
                        numberOfIterations = solver.iterations();
                    } else {
                        STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with DGMRES with identity preconditioner (Eigen library).");

                        StormEigen::DGMRES<StormEigen::SparseMatrix<ValueType>, StormEigen::IdentityPreconditioner> solver;
                        solver.setTolerance(precision);
                        solver.setMaxIterations(maxIter);
                        solver.set_restart(restartThreshold);
                        solver.compute(eigenMatrix);
                        eigenX = solver.solveWithGuess(eigenB, eigenX);
                        converged = solver.info() == StormEigen::ComputationInfo::Success;
                        numberOfIterations = solver.iterations();
                    }
                } else if (solutionMethod == EigenLinearEquationSolverMethod::Gmres) {
                    if (preconditioner == EigenLinearEquationSolverPreconditioner::Ilu) {
                        STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with GMRES with Ilu preconditioner (Eigen library).");

                        StormEigen::GMRES<StormEigen::SparseMatrix<ValueType>, StormEigen::IncompleteLUT<ValueType>> solver;
                        solver.setTolerance(precision);
                        solver.setMaxIterations(maxIter);
                        solver.set_restart(restartThreshold);
                        solver.compute(eigenMatrix);
                        eigenX = solver.solveWithGuess(eigenB, eigenX);
                        converged = solver.info() == StormEigen::ComputationInfo::Success;
                        numberOfIterations = solver.iterations();
                    } else if (preconditioner == EigenLinearEquationSolverPreconditioner::Diagonal) {
                        STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with GMRES with Diagonal preconditioner (Eigen library).");

                        StormEigen::GMRES<StormEigen::SparseMatrix<ValueType>, StormEigen::DiagonalPreconditioner<ValueType>> solver;
                        solver.setTolerance(precision);
                        solver.setMaxIterations(maxIter);
                        solver.set_restart(restartThreshold);
                        solver.compute(eigenMatrix);
                        eigenX = solver.solveWithGuess(eigenB, eigenX);
                        converged = solver.info() == StormEigen::ComputationInfo::Success;
                        numberOfIterations = solver.iterations();
                    } else {
                        STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with GMRES with identity preconditioner (Eigen library).");

                        StormEigen::GMRES<StormEigen::SparseMatrix<ValueType>, StormEigen::IdentityPreconditioner> solver;
                        solver.setTolerance(precision);
                        solver.setMaxIterations(maxIter);
                        solver.set_restart(restartThreshold);
                        solver.compute(eigenMatrix);
                        eigenX = solver.solveWithGuess(eigenB, eigenX);
                        converged = solver.info() == StormEigen::ComputationInfo::Success;
                        numberOfIterations = solver.iterations();
                    }
                }
                
                // Make sure that all results conform to the (global) bounds.
                storm::utility::vector::clip(x, this->lowerBound, this->upperBound);
                
                // Check if the solver converged and issue a warning otherwise.
                if (converged) {
                    STORM_LOG_INFO("Iterative solver converged after " << numberOfIterations << " iterations.");
                    return true;
                } else {
                    STORM_LOG_WARN("Iterative solver did not converge.");
                    return false;
                }
            }
            
            return true;
        }
        
        template<typename ValueType>
        void EigenLinearEquationSolver<ValueType>::multiply(std::vector<ValueType>& x, std::vector<ValueType> const* b, std::vector<ValueType>& result) const {
            // Typedef the map-type so we don't have to spell it out.
            typedef decltype(StormEigen::Matrix<ValueType, StormEigen::Dynamic, 1>::Map(b->data(), b->size())) MapType;

            auto eigenX = StormEigen::Matrix<ValueType, StormEigen::Dynamic, 1>::Map(x.data(), x.size());
            auto eigenResult = StormEigen::Matrix<ValueType, StormEigen::Dynamic, 1>::Map(result.data(), result.size());

            std::unique_ptr<MapType> eigenB;
            if (b != nullptr) {
                eigenB = std::make_unique<MapType>(StormEigen::Matrix<ValueType, StormEigen::Dynamic, 1>::Map(b->data(), b->size()));
            }
            
            if (&x != &result) {
                if (b != nullptr) {
                    eigenResult.noalias() = *eigenA * eigenX + *eigenB;
                } else {
                    eigenResult.noalias() = *eigenA * eigenX;
                }
            } else {
                if (b != nullptr) {
                    eigenResult = *eigenA * eigenX + *eigenB;
                } else {
                    eigenResult = *eigenA * eigenX;
                }
            }
        }
        
        template<typename ValueType>
        LinearEquationSolverProblemFormat EigenLinearEquationSolver<ValueType>::getEquationProblemFormat(Environment const& env) const {
            return LinearEquationSolverProblemFormat::EquationSystem;
        }
        
        template<typename ValueType>
        uint64_t EigenLinearEquationSolver<ValueType>::getMatrixRowCount() const {
            return eigenA->rows();
        }
        
        template<typename ValueType>
        uint64_t EigenLinearEquationSolver<ValueType>::getMatrixColumnCount() const {
            return eigenA->cols();
        }

        template<typename ValueType>
        std::unique_ptr<storm::solver::LinearEquationSolver<ValueType>> EigenLinearEquationSolverFactory<ValueType>::create(Environment const& env, LinearEquationSolverTask const& task) const {
            return std::make_unique<storm::solver::EigenLinearEquationSolver<ValueType>>();
        }
        
        template<typename ValueType>
        std::unique_ptr<LinearEquationSolverFactory<ValueType>> EigenLinearEquationSolverFactory<ValueType>::clone() const {
            return std::make_unique<EigenLinearEquationSolverFactory<ValueType>>(*this);
        }
        
        template class EigenLinearEquationSolver<double>;
        template class EigenLinearEquationSolverFactory<double>;

#ifdef STORM_HAVE_CARL
        template class EigenLinearEquationSolver<storm::RationalNumber>;
        template class EigenLinearEquationSolver<storm::RationalFunction>;
        
        template class EigenLinearEquationSolverFactory<storm::RationalNumber>;
        template class EigenLinearEquationSolverFactory<storm::RationalFunction>;
#endif
    }
}
