#include <memory>
#include <type_traits>

#include "storm/solver/GaussEliminationLinearEquationSolver.h"

#include "gausselimination/FractionfreeWorker.h"
#include "gausselimination/GaussElimination.h"
#include "gausselimination/RationalFunctionComplexity.h"
#include "gausselimination/StandardWorker.h"
#include "gausselimination/StandardWorkerOtherValueType.h"
#include "storm/storage/DenseMatrix.h"
#include "storm/storage/PermutatedDenseMatrix.h"
#include "storm/storage/TopologicalMatrixReordering.h"
#include "storm/utility/NumericalTypesConverter.h"
#include "storm/settings/modules/GaussEliminationSettings.h"
#include "storm/settings/modules/CoreSettings.h"
#include "storm/settings/SettingsManager.h"

#include "storm/utility/NumericalTypes.h"

#include "storm/storage/StronglyConnectedComponentDecomposition.h"

#include "storm/utility/Stopwatch.h"

//#include "storm/utility/graph.h"
#include "storm/utility/macros.h"
#include "storm/utility/vector.h"
//#include "storm/utility/stateelimination.h"
#include "storm/utility/GCDLog.h"

namespace gauss = storm::solver::gausselimination;

namespace storm {
namespace solver {

template <typename ValueType>
GaussEliminationLinearEquationSolver<ValueType>::GaussEliminationLinearEquationSolver() : localA(nullptr), A(nullptr) {
}

template<typename ValueType>
void GaussEliminationLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType> const& A) {
    this->A = &A;
    localA.reset();
    this->clearCache();
}

template<typename ValueType>
void GaussEliminationLinearEquationSolver<ValueType>::setMatrix(storm::storage::SparseMatrix<ValueType>&& A) {
    localA = std::make_unique<storm::storage::SparseMatrix<ValueType>>(std::move(A));
    this->A = localA.get();
    this->clearCache();
}

template<typename ValueType>
uint64_t GaussEliminationLinearEquationSolver<ValueType>::getMatrixRowCount() const {
    return this->A->getRowCount();
}

template<typename ValueType>
uint64_t GaussEliminationLinearEquationSolver<ValueType>::getMatrixColumnCount() const {
    return this->A->getColumnCount();
}

template<typename ValueType>
std::unique_ptr<storm::solver::LinearEquationSolver<ValueType>> GaussEliminationLinearEquationSolverFactory<ValueType>::create(Environment const& env, LinearEquationSolverTask const& task) const {
    return std::make_unique<GaussEliminationLinearEquationSolver<ValueType>>();
}

template<typename ValueType>
std::unique_ptr<LinearEquationSolverFactory<ValueType>> GaussEliminationLinearEquationSolverFactory<ValueType>::clone() const {
    return std::make_unique<GaussEliminationLinearEquationSolverFactory<ValueType>>(*this);
}

template<typename ValueType>
LinearEquationSolverProblemFormat GaussEliminationLinearEquationSolver<ValueType>::getEquationProblemFormat(Environment const& env) const {
    return LinearEquationSolverProblemFormat::EquationSystem;
}
  
template<typename ValueType>
void GaussEliminationLinearEquationSolver<ValueType>::multiply(std::vector<ValueType>& x, std::vector<ValueType> const* b, std::vector<ValueType>& result) const {
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

template <typename ValueType>
template <typename OtherValueType>
void GaussEliminationLinearEquationSolver<ValueType>::assignAndConvertIfNecessary(std::vector<ValueType>& a, std::vector<OtherValueType>& b) const {
    storm::utility::NumericalTypesConverter<OtherValueType, ValueType> converter;

    storm::utility::Stopwatch timer(true);
    for (std::size_t i = 0; i < b.size(); i++) {
        a[i] = std::move(converter.convert(b[i]));
    }
    timer.stop();
    STORM_LOG_INFO("Converting solution vector from " << storm::utility::NumericalTypes::getTypeName<OtherValueType>() << " to " << storm::utility::NumericalTypes::getTypeName<ValueType>() << ": " << timer);

    if (storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().isReportReductionsSet()) {
        STORM_LOG_INFO(gauss::RationalFunctionComplexity::reportMaxReduction(a, b));
    }

    b.clear();
}

template <typename ValueType>
void GaussEliminationLinearEquationSolver<ValueType>::assignAndConvertIfNecessary(std::vector<ValueType>& a, std::vector<ValueType>& b) const {
    a.swap(b);
    b.clear();
}

template <typename SourceValueType, typename TargetValueType>
static storm::storage::DenseMatrix<TargetValueType> generateAugmentedDenseMatrix(const storm::storage::SparseMatrix<SourceValueType>& A, const std::vector<SourceValueType>& b) {
    return storm::storage::DenseMatrix<TargetValueType>::augmentedMatrix(A, b);
}


template <typename SourceValueType, typename TargetValueType>
static storm::storage::FlexibleSparseMatrix<TargetValueType> generateAugmentedFlexibleSparseMatrix(const storm::storage::SparseMatrix<SourceValueType>& A, const std::vector<SourceValueType>& b) {
    storm::storage::FlexibleSparseMatrix<TargetValueType> result;

    storm::utility::NumericalTypesConverter<SourceValueType,TargetValueType> converter;
    return storm::storage::FlexibleSparseMatrix<TargetValueType>::augmentedFromSparseMatrix(A, b);
}

template <typename ValueType>
template <typename Worker>
void GaussEliminationLinearEquationSolver<ValueType>::solve(std::vector<ValueType>& x, std::vector<ValueType> const& b) const {
    typedef typename Worker::DenseValueType DenseValueType;
    typedef typename Worker::IntermediateValueType IntermediateValueType;

    STORM_LOG_INFO("Solving linear equation system (" << x.size() << " rows) with Gauss elimination");

    bool debug = storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().isDebug();

    if (debug) {
        std::cout << "Input matrix:\n" << *A;
        std::cout << "\nb:\n" << storm::utility::vector::toString(b) << "\n";
    }

    std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport;
    if (storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().isExportHTMLSet()) {
        std::string filename = storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().getExportHTMLFilename();
        STORM_LOG_INFO("Logging Gaussian elimination steps to " << filename);
        htmlExport.reset(new storm::storage::HTMLMatrixExport(filename));
        htmlExport->setPrintComplexityOption(storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().isExportHTMLComplexitySet());
        htmlExport->setPrintStatsOption(storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().isExportHTMLStatsSet());

        htmlExport->printHeader("Gaussian Elimination", htmlExport->getStandardCSS());
        htmlExport->print(std::string("<p>Method: ") + Worker::getMethodName() + "</p>\n");
    }

    switch (storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().getMatrix()) {
    case storm::settings::modules::GaussEliminationSettings::Matrix::Dense: {
        STORM_LOG_INFO("Converting sparse matrix to dense augmented matrix...");
        storm::utility::Stopwatch timer(true);
        storm::storage::DenseMatrix<DenseValueType> matrix = generateAugmentedDenseMatrix<ValueType, DenseValueType>(*A, b);
        timer.stop();
        STORM_LOG_INFO("Converting sparse matrix to dense augmented matrix: " << timer);
        timer.reset();

        if (debug) {
            std::cout << "Augmented matrix:\n" << matrix;
        }

        solve<Worker>(matrix, x, htmlExport, debug);
        break;
    }
    case storm::settings::modules::GaussEliminationSettings::Matrix::Sparse: {
        STORM_LOG_INFO("Converting sparse matrix to flexible sparse augmented matrix...");
        storm::utility::Stopwatch timer(true);
        storm::storage::FlexibleSparseMatrix<DenseValueType> matrix = generateAugmentedFlexibleSparseMatrix<ValueType,DenseValueType>(*A, b);
        timer.stop();
        STORM_LOG_INFO("Converting sparse matrix to flexible sparse augmented matrix: " << timer);
        timer.reset();

        if (debug) {
            std::cout << "Augmented matrix:\n" << matrix;
        }

        solve<Worker>(matrix, x, htmlExport, debug);
        break;
    }
    default:
        STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "");
    }
}

template <typename ValueType>
template <typename Worker, typename Matrix>
void GaussEliminationLinearEquationSolver<ValueType>::solve(Matrix& augmentedMatrix, std::vector<ValueType>& x, std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport, bool debug) const {
    typedef typename Worker::DenseValueType DenseValueType;
    typedef typename Worker::IntermediateValueType IntermediateValueType;

    std::vector<IntermediateValueType> x_(x.size());

    bool topologicalOrdering = storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().isTopologicalOrderingSet();

    if (topologicalOrdering) {
        if (htmlExport) {
            htmlExport->slideBegin();
            htmlExport->print("<h2>Unsorted (before topological ordering)</h2>");
            htmlExport->printMatrix(augmentedMatrix);
            htmlExport->print("<hr>");
            htmlExport->slideEnd();
        }

        STORM_LOG_INFO("Obtaining SCC decomposition and topological ordering...");

        storm::utility::Stopwatch timerSCC(true);
        typedef storm::storage::StronglyConnectedComponentDecomposition<ValueType> SCCDecomposition;
        SCCDecomposition sccs(*A);
        auto sorted = storm::storage::TopologicalMatrixReordering::topologicallySorted(augmentedMatrix, sccs);

        auto& sortedMatrix = std::get<0>(sorted);
        const auto& sccInfo = *std::get<1>(sorted).sccInfo;
        const auto& permutation = *std::get<1>(sorted).permutation;

        timerSCC.stop();
        STORM_LOG_INFO("SCC decomposition and topological ordering: " << timerSCC);
        STORM_LOG_INFO(sccInfo.getStatistics());

        if (htmlExport) {
            // set scc info for sortedMatrix
            htmlExport->setSCCInfo(sccInfo);
        }

#ifdef STORM_HAVE_GCD_LOG
        storm::utility::GCDLog::enable();
        STORM_LOG_INFO("Enabling GCD logging...\n");
#endif
	
        storm::utility::Stopwatch timerCalc(true);
        solve<Worker>(sortedMatrix, sccInfo, x_, htmlExport, debug);

#ifdef STORM_HAVE_GCD_LOG
        storm::utility::GCDLog::disable();
        STORM_LOG_INFO("Disabling GCD logging...\n");
#endif
	
        timerCalc.stop();
        STORM_LOG_INFO("Time for parametric computation: " << timerCalc);

        STORM_LOG_INFO("Reorder solution vector to correspond to model ordering...");
        // un-permutate solution vector
        permutation.permuteInv(x_);

        if (htmlExport) {
            htmlExport->resetSCCInfo();

            htmlExport->slideBegin();
            htmlExport->print(std::string("<h2>Solution [") + storm::utility::NumericalTypes::getTypeName<IntermediateValueType>() + "] (after reverting topological ordering)</h2>");
            htmlExport->printVector(x_);
            htmlExport->slideEnd();
        }
    } else {

#ifdef STORM_HAVE_GCD_LOG
        storm::utility::GCDLog::enable();
        STORM_LOG_INFO("Enabling GCD logging...\n");
#endif

        storm::utility::Stopwatch timerCalc(true);	
        solve<Worker>(augmentedMatrix, boost::optional<const storm::storage::SCCInfo&>(), x_, htmlExport, debug);

#ifdef STORM_HAVE_GCD_LOG
        storm::utility::GCDLog::disable();
        STORM_LOG_INFO("Disabling GCD logging...\n");
#endif

        timerCalc.stop();
        STORM_LOG_INFO("Time for parametric computation: " << timerCalc);

        if (htmlExport) {
            htmlExport->slideBegin();
            htmlExport->print(std::string("<h2>Solution [") + storm::utility::NumericalTypes::getTypeName<IntermediateValueType>() + "]</h2>");
            htmlExport->printVector(x_);
            htmlExport->slideEnd();
        }
    }

    if (storm::settings::getModule<storm::settings::modules::CoreSettings>().isResultStatsSet()) {
        std::cout << storm::utility::NumericalTypes::getStats(x_);
    }

    // convert x_ solution vector to the standard ValueType used by Storm
    assignAndConvertIfNecessary(x, x_);

    if (!(std::is_same<ValueType, IntermediateValueType>::value)) {
        if (htmlExport) {
            htmlExport->slideBegin();
            htmlExport->print(std::string("<h2>Solution [") + storm::utility::NumericalTypes::getTypeName<ValueType>() + "] (after conversion to standard Storm data type)</h2>");
            htmlExport->printVector(x);
            htmlExport->slideEnd();
        }
    }

    if (htmlExport) {
        htmlExport->printFooter();
        htmlExport->close();
    }

}

template <typename ValueType>
template <typename Worker, typename Matrix, typename OtherValueType>
void GaussEliminationLinearEquationSolver<ValueType>::solve(Matrix& augmentedMatrix, boost::optional<const storm::storage::SCCInfo&> sccInfo, std::vector<OtherValueType>& x, std::shared_ptr<storm::storage::HTMLMatrixExport> htmlExport, bool debug) const {
    Worker worker;

    worker.prepare(augmentedMatrix, sccInfo);
    STORM_LOG_INFO("Gaussian Elemination via " << worker.getMethodNameAndConfiguration());

    STORM_LOG_INFO("Converting matrix to triangular form...");
    storm::utility::Stopwatch timer(true);
    gauss::GaussElimination<Matrix, Worker>::eliminate(augmentedMatrix, sccInfo, worker, debug, htmlExport);
    timer.stop();
    STORM_LOG_INFO("Converting matrix to triangular form: " << timer);

    if (debug) {
        std::cout << "\nTriangulated matrix:\n" << augmentedMatrix;
    }

    STORM_LOG_INFO("Obtaining solution vector by back-substitution...");
    timer.reset();
    timer.start();
    gauss::GaussElimination<Matrix, Worker>::solve(augmentedMatrix, worker, x, debug);
    timer.stop();
    STORM_LOG_INFO("Obtaining solution vector by back-substitution: "  << timer);
    if (debug) std::cout << "\nx (after solving):\n" << storm::utility::vector::toString(x) << "\n";

    if (htmlExport) {
        htmlExport->slideBegin();
        htmlExport->print(std::string("<h2>Solution [") + storm::utility::NumericalTypes::getTypeName<OtherValueType>() + "]</h2>");
        htmlExport->printVector(x);
        htmlExport->slideEnd();
    }
}


template <>
bool GaussEliminationLinearEquationSolver<double>::internalSolveEquations(Environment const& env, std::vector<double>& x, std::vector<double> const& b) const {
    // specialization for double
    STORM_LOG_INFO("Performing standard Gauss elmination using double values...");
    solve<gauss::StandardWorker<double>>(x, b);

    return true;
}



template <typename ValueType>
bool GaussEliminationLinearEquationSolver<ValueType>::internalSolveEquations(Environment const& env, std::vector<ValueType>& x, std::vector<ValueType> const& b) const {
    // ValueType is either storm::RationalNumberCoefficient or storm::RationalFunction here

    switch (storm::settings::getModule<storm::settings::modules::GaussEliminationSettings>().getMethod()) {
    case storm::settings::modules::GaussEliminationSettings::Method::Standard:
        solve<gauss::StandardWorker<ValueType>>(x, b);
        break;
    case storm::settings::modules::GaussEliminationSettings::Method::Factorized: {
        if (std::is_same<ValueType, storm::RationalFunction>::value) {
            solve<gauss::StandardWorker<ValueType>>(x, b);
        } else {
            solve<gauss::StandardWorkerOtherValueType<ValueType, storm::RationalFunction>>(x, b);
        }
        break;
    }
    case storm::settings::modules::GaussEliminationSettings::Method::FactorizedNoGcd: {
        solve<gauss::StandardWorkerOtherValueType<ValueType, storm::FactorizedRationalFunctionNoSimplify>>(x, b);
        break;
    }
    case storm::settings::modules::GaussEliminationSettings::Method::Plain: {
        solve<gauss::StandardWorkerOtherValueType<ValueType, storm::PlainRationalFunction>>(x, b);
        break;
    }
    case storm::settings::modules::GaussEliminationSettings::Method::PlainNoGcd: {
        solve<gauss::StandardWorkerOtherValueType<ValueType, storm::PlainRationalFunctionNoSimplify>>(x, b);
        break;
    }
    case storm::settings::modules::GaussEliminationSettings::Method::Fractionfree: {
        typedef storm::RawPolynomial DenseValueType;
        typedef storm::PlainRationalFunctionNoSimplify IntermediateValueType;

        solve<gauss::FractionfreeWorker<ValueType, IntermediateValueType, DenseValueType>>(x, b);
        break;
    }
    default:
        STORM_LOG_THROW(false, storm::exceptions::IllegalArgumentException, "Unsupported method");
    }

    return true;
}



// instantiations

template class GaussEliminationLinearEquationSolver<double>;
template class GaussEliminationLinearEquationSolverFactory<double>;

#ifdef STORM_HAVE_CARL
//template class GaussEliminationLinearEquationSolver<storm::RationalFunctionCoefficient>;
template class GaussEliminationLinearEquationSolver<storm::RationalFunction>;

//template class GaussEliminationLinearEquationSolverFactory<storm::RationalFunctionCoefficient>;
template class GaussEliminationLinearEquationSolverFactory<storm::RationalFunction>;
#endif


}
}

