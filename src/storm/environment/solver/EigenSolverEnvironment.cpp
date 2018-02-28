#include "storm/environment/solver/EigenSolverEnvironment.h"

#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/EigenEquationSolverSettings.h"
#include "storm/utility/constants.h"
#include "storm/utility/macros.h"

namespace storm {
    
    EigenSolverEnvironment::EigenSolverEnvironment() {
        auto const& eigenSettings = storm::settings::getModule<storm::settings::modules::EigenEquationSolverSettings>();
        
        method = eigenSettings.getLinearEquationSystemMethod();
        methodSetFromDefault = eigenSettings.isLinearEquationSystemMethodSetFromDefault();
        preconditioner = eigenSettings.getPreconditioningMethod();
        restartThreshold = eigenSettings.getRestartIterationCount();
        maxIterationCount = eigenSettings.getMaximalIterationCount();
        precision = storm::utility::convertNumber<storm::RationalNumber>(eigenSettings.getPrecision());
        topologicalOrdering = eigenSettings.isTopologicalOrderingSet();

    }

    EigenSolverEnvironment::~EigenSolverEnvironment() {
        // Intentionally left empty
    }
    
    storm::solver::EigenLinearEquationSolverMethod const& EigenSolverEnvironment::getMethod() const {
        return method;
    }
    
    bool EigenSolverEnvironment::isMethodSetFromDefault() const {
        return methodSetFromDefault;
    }
    
    void EigenSolverEnvironment::setMethod(storm::solver::EigenLinearEquationSolverMethod value) {
        methodSetFromDefault = false;
        method = value;
    }
    
    storm::solver::EigenLinearEquationSolverPreconditioner const& EigenSolverEnvironment::getPreconditioner() const {
        return preconditioner;
    }
    
    void EigenSolverEnvironment::setPreconditioner(storm::solver::EigenLinearEquationSolverPreconditioner value) {
        preconditioner = value;
    }
    
    uint64_t const& EigenSolverEnvironment::getRestartThreshold() const {
        return maxIterationCount;
    }
    
    void EigenSolverEnvironment::setRestartThreshold(uint64_t value) {
        maxIterationCount = value;
    }
    
    uint64_t const& EigenSolverEnvironment::getMaximalNumberOfIterations() const {
        return maxIterationCount;
    }
    
    void EigenSolverEnvironment::setMaximalNumberOfIterations(uint64_t value) {
        maxIterationCount = value;
    }
    
    storm::RationalNumber const& EigenSolverEnvironment::getPrecision() const {
        return precision;
    }
    
    void EigenSolverEnvironment::setPrecision(storm::RationalNumber value) {
        precision = value;
    }

    bool EigenSolverEnvironment::isTopologicalOrderingSet() const {
        return topologicalOrdering;
    }

    void EigenSolverEnvironment::setTopologicalOrdering(bool value) {
        topologicalOrdering = value;
    }

}
