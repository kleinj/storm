#include "storm/solver/EliminationLinearEquationSolver.h"


#include <numeric>
#include "storm/storage/HTMLMatrixExport.h"
#include "storm/storage/StronglyConnectedComponentDecomposition.h"
#include "storm/storage/TopologicalMatrixReordering.h"
#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/EliminationSettings.h"
// #include "storm/settings/modules/ParametricSettings.h"
#include "storm/settings/modules/CoreSettings.h"

#include "storm/solver/stateelimination/StatePriorityQueue.h"
#include "storm/solver/stateelimination/PrioritizedStateEliminator.h"
#include "storm/utility/NumericalTypes.h"
#include "storm/utility/Stopwatch.h"
#include "storm/utility/graph.h"
#include "storm/utility/macros.h"
#include "storm/utility/vector.h"
#include "storm/utility/stateelimination.h"
#include "storm/utility/GCDLog.h"

namespace storm {
    namespace solver {
        
        using namespace stateelimination;
        using namespace storm::utility::stateelimination;
        
        template<typename ValueType>
        EliminationLinearEquationSolver<ValueType>::EliminationLinearEquationSolver() {
            // Intentionally left empty.
        }
        
        template<typename ValueType>
        EliminationLinearEquationSolver<ValueType>::EliminationLinearEquationSolver(storm::storage::SparseMatrix<ValueType> const& A) : localA(nullptr), A(nullptr) {
            this->setMatrix(A);
        }
        
        template<typename ValueType>
        EliminationLinearEquationSolver<ValueType>::EliminationLinearEquationSolver(storm::storage::SparseMatrix<ValueType>&& A) : localA(nullptr), A(nullptr) {
            this->setMatrix(std::move(A));
        }
        
        template<typename ValueType>
        void EliminationLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType> const& A) {
            this->A = &A;
            localA.reset();
            this->clearCache();
        }
        
        template<typename ValueType>
        void EliminationLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType>&& A) {
            localA = std::make_unique<storm::storage::SparseMatrix<ValueType>>(std::move(A));
            this->A = localA.get();
            this->clearCache();
        }

        template<typename ValueType>
        bool EliminationLinearEquationSolver<ValueType>::internalSolveEquations(Environment const& env, std::vector<ValueType>& x, std::vector<ValueType> const& b) const {
            // FIXME: This solver will not work for all input systems. More concretely, the current implementation will
            // not work for systems that have a 1 on the diagonal. This is not a restriction of this technique in general
            // but arbitrary matrices require pivoting, which is not currently implemented.
            
            STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with elimination");

            std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport;
            if (storm::settings::getModule<storm::settings::modules::EliminationSettings>().isExportHTMLSet()) {
                std::string filename = storm::settings::getModule<storm::settings::modules::EliminationSettings>().getExportHTMLFilename();
                STORM_LOG_INFO("Logging state elimination steps to " << filename);
                htmlExport.reset(new storm::storage::HTMLMatrixExport(filename));
                htmlExport->setPrintComplexityOption(storm::settings::getModule<storm::settings::modules::EliminationSettings>().isExportHTMLComplexitySet());
                htmlExport->setPrintStatsOption(storm::settings::getModule<storm::settings::modules::EliminationSettings>().isExportHTMLStatsSet());

                htmlExport->printHeader("State Elimination", htmlExport->getStandardCSS());
                htmlExport->print(std::string("<p>Order: ")
                     + storm::settings::getModule<storm::settings::modules::EliminationSettings>().getEliminationOrderAsString()
                     + ", Method: " + storm::settings::getModule<storm::settings::modules::EliminationSettings>().getEliminationMethodAsString()
                     + "</p>\n");
            }

            storm::storage::SparseMatrix<ValueType> const& transitionMatrix = localA ? *localA : *A;
            // Initialize the solution to the right-hand side of the equation system.
            x = b;

            bool topologicalOrdering = storm::settings::getModule<storm::settings::modules::EliminationSettings>().isTopologicalOrderingSet();

            if (topologicalOrdering) {
                if (htmlExport) {
                    htmlExport->slideBegin();
                    htmlExport->print("<h2>Unsorted (before topological ordering)</h2>");
                    htmlExport->printMatrix(transitionMatrix, x);
                    htmlExport->print("<hr>");
                    htmlExport->slideEnd();
                }

                STORM_LOG_INFO("Obtaining SCC decomposition and topological ordering...");

                storm::utility::Stopwatch timerSCC(true);
                typedef storm::storage::StronglyConnectedComponentDecomposition<ValueType> SCCDecomposition;
                SCCDecomposition sccs(transitionMatrix);
                auto sorted = storm::storage::TopologicalMatrixReordering::topologicallySorted(transitionMatrix, sccs);

                auto& sortedMatrix = std::get<0>(sorted);
                const auto& sccInfo = *std::get<1>(sorted).sccInfo;
                const auto& permutation = *std::get<1>(sorted).permutation;

                permutation.permute(x);
                std::vector<ValueType> b_(b);
                permutation.permute(b_);

                timerSCC.stop();
                STORM_LOG_INFO("SCC decomposition and topological ordering: " << timerSCC);
                STORM_LOG_INFO(sccInfo.getStatistics());

                if (htmlExport) {
                    // set scc info for sortedMatrix
                    htmlExport->setSCCInfo(sccInfo);
                }

                bool rv = internalSolve(env, sortedMatrix, b_, x, htmlExport);

                STORM_LOG_INFO("Reorder solution vector to correspond to model ordering...");
                // un-permutate solution vector
                permutation.permuteInv(x);

                if (htmlExport) {
                    htmlExport->resetSCCInfo();
                    htmlExport->slideBegin();
                    htmlExport->print("<h2>Solution (after reverting topological ordering)</h2>");
                    htmlExport->printVector(x);
                    htmlExport->slideEnd();
                }

                return rv;
            } else {
                return internalSolve(env, transitionMatrix, b, x, htmlExport);
            }
        }

        template <typename ValueType>
        bool EliminationLinearEquationSolver<ValueType>::internalSolve(Environment const& env, const storm::storage::SparseMatrix<ValueType>& transitionMatrix, const std::vector<ValueType>& b, std::vector<ValueType>& x, std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport) const {
            storm::storage::SparseMatrix<ValueType> backwardTransitions = transitionMatrix.transpose();

            // Translate the matrix and its transpose into the flexible format.
            storm::storage::FlexibleSparseMatrix<ValueType> flexibleMatrix(transitionMatrix, false);
            storm::storage::FlexibleSparseMatrix<ValueType> flexibleBackwardTransitions(backwardTransitions, true);

            if (htmlExport) {
                htmlExport->slideBegin();
                htmlExport->print("<h2>Initial</h2>");
                htmlExport->printMatrix(flexibleMatrix, x);
                htmlExport->print("<hr>");
                htmlExport->slideEnd();
            }

            storm::utility::Stopwatch timer(true);

            boost::optional<std::vector<uint_fast64_t>> distanceBasedPriorities;
            
            // TODO: get the order from the environment
            storm::settings::modules::EliminationSettings::EliminationOrder order = storm::settings::getModule<storm::settings::modules::EliminationSettings>().getEliminationOrder();
            
            if (eliminationOrderNeedsDistances(order)) {
                // Since we have no initial states at this point, we determine a representative of every BSCC regarding
                // the backward transitions, because this means that every row is reachable from this set of rows, which
                // we require to make sure we cover every row.
                storm::storage::BitVector initialRows = storm::utility::graph::getBsccCover(backwardTransitions);
                distanceBasedPriorities = getDistanceBasedPriorities(transitionMatrix, backwardTransitions, initialRows, b, eliminationOrderNeedsForwardDistances(order), eliminationOrderNeedsReversedDistances(order));
            }
            
            std::shared_ptr<StatePriorityQueue> priorityQueue = createStatePriorityQueue<ValueType>(distanceBasedPriorities, flexibleMatrix, flexibleBackwardTransitions, b, storm::storage::BitVector(x.size(), true));
            
#ifdef STORM_HAVE_GCD_LOG
            storm::utility::GCDLog::enable();
	    STORM_LOG_INFO("Enabling GCD logging...\n");
#endif


            // Create a state eliminator to perform the actual elimination.
            PrioritizedStateEliminator<ValueType> eliminator(flexibleMatrix, flexibleBackwardTransitions, priorityQueue, x);
            
            // Eliminate all states.
            while (priorityQueue->hasNext()) {
                auto state = priorityQueue->pop();
                eliminator.eliminateState(state, false);

                if (htmlExport) {
                    htmlExport->slideBegin();
                    htmlExport->print("<h2>After elimination of state " + std::to_string(state) +"</h2>");
                    htmlExport->printMatrix(flexibleMatrix, x, boost::none, state);
                    htmlExport->print("<hr>");
                    htmlExport->slideEnd();
                }
            }
            
#ifdef STORM_HAVE_GCD_LOG
            storm::utility::GCDLog::disable();
	    STORM_LOG_INFO("Disabling GCD logging...\n");
#endif

            timer.stop();
            STORM_LOG_INFO("Time for parametric computation: " << timer);

            if (storm::settings::getModule<storm::settings::modules::CoreSettings>().isResultStatsSet()) {
                std::cout << storm::utility::NumericalTypes::getStats(x);
            }

            if (htmlExport) {
                htmlExport->slideBegin();
                htmlExport->print("<h2>Solution</h2>");
                htmlExport->printVector(x);
                htmlExport->slideEnd();
            }

            return true;
        }
        
        template<typename ValueType>
        void EliminationLinearEquationSolver<ValueType>::multiply(std::vector<ValueType>& x, std::vector<ValueType> const* b, std::vector<ValueType>& result) const {
            if (&x != &result) {
                A->multiplyWithVector(x, result);
                if (b != nullptr) {
                    storm::utility::vector::addVectors(result, *b, result);
                }
            } else {
                // If the two vectors are aliases, we need to create a temporary.
                std::vector<ValueType> tmp(result.size());
                A->multiplyWithVector(x, tmp);
                if (b != nullptr) {
                    storm::utility::vector::addVectors(tmp, *b, result);
                }
            }
        }
        
        template<typename ValueType>
        LinearEquationSolverProblemFormat EliminationLinearEquationSolver<ValueType>::getEquationProblemFormat(Environment const& env) const {
            return LinearEquationSolverProblemFormat::FixedPointSystem;
        }
        
        template<typename ValueType>
        uint64_t EliminationLinearEquationSolver<ValueType>::getMatrixRowCount() const {
            return this->A->getRowCount();
        }
        
        template<typename ValueType>
        uint64_t EliminationLinearEquationSolver<ValueType>::getMatrixColumnCount() const {
            return this->A->getColumnCount();
        }
        
        template<typename ValueType>
        std::unique_ptr<storm::solver::LinearEquationSolver<ValueType>> EliminationLinearEquationSolverFactory<ValueType>::create(Environment const& env, LinearEquationSolverTask const& task) const {
            return std::make_unique<storm::solver::EliminationLinearEquationSolver<ValueType>>();
        }
        
        template<typename ValueType>
        std::unique_ptr<LinearEquationSolverFactory<ValueType>> EliminationLinearEquationSolverFactory<ValueType>::clone() const {
            return std::make_unique<EliminationLinearEquationSolverFactory<ValueType>>(*this);
        }
        
        template class EliminationLinearEquationSolver<double>;
        template class EliminationLinearEquationSolverFactory<double>;

#ifdef STORM_HAVE_CARL
        
        template class EliminationLinearEquationSolver<storm::RationalNumber>;
        template class EliminationLinearEquationSolver<storm::RationalFunction>;
        
        template class EliminationLinearEquationSolverFactory<storm::RationalNumber>;
        template class EliminationLinearEquationSolverFactory<storm::RationalFunction>;
#endif
    }
}

