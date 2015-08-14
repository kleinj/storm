#include "src/solver/GurobiLpSolver.h"

#ifdef STORM_HAVE_GUROBI
#include <numeric>

#include "src/storage/expressions/LinearCoefficientVisitor.h"

#include "src/settings/SettingsManager.h"
#include "src/settings/modules/DebugSettings.h"
#include "src/settings/modules/GurobiSettings.h"

#include "src/utility/macros.h"
#include "src/storage/expressions/Expression.h"
#include "src/storage/expressions/ExpressionManager.h"

#include "src/exceptions/InvalidStateException.h"
#include "src/exceptions/InvalidAccessException.h"
#include "src/exceptions/InvalidArgumentException.h"

namespace storm {
    namespace solver {
        
        GurobiLpSolver::GurobiLpSolver(std::string const& name, ModelSense const& modelSense) : LpSolver(modelSense), env(nullptr), model(nullptr), nextVariableIndex(0) {
            // Create the environment.
            int error = GRBloadenv(&env, "");
            if (error || env == nullptr) {
				LOG4CPLUS_ERROR(logger, "Could not initialize Gurobi (" << GRBgeterrormsg(env) << ", error code " << error << ").");
				throw storm::exceptions::InvalidStateException() << "Could not initialize Gurobi environment (" << GRBgeterrormsg(env) << ", error code " << error << ").";
            }
            
            // Set some general properties of the environment.
            setGurobiEnvironmentProperties();
            
            // Create the model.
            error = GRBnewmodel(env, &model, name.c_str(), 0, nullptr, nullptr, nullptr, nullptr, nullptr);
            if (error) {
				LOG4CPLUS_ERROR(logger, "Could not initialize Gurobi model (" << GRBgeterrormsg(env) << ", error code " << error << ").");
				throw storm::exceptions::InvalidStateException() << "Could not initialize Gurobi model (" << GRBgeterrormsg(env) << ", error code " << error << ").";
            }
        }
        
        GurobiLpSolver::GurobiLpSolver(std::string const& name) : GurobiLpSolver(name, ModelSense::Minimize) {
            // Intentionally left empty.
        }
        
        GurobiLpSolver::GurobiLpSolver(ModelSense const& modelSense) : GurobiLpSolver("", modelSense) {
            // Intentionally left empty.
        }
        
        GurobiLpSolver::GurobiLpSolver() : GurobiLpSolver("", ModelSense::Minimize) {
            // Intentionally left empty.
        }
        
        GurobiLpSolver::~GurobiLpSolver() {
            // Dispose of the objects allocated inside Gurobi.
            GRBfreemodel(model);
            GRBfreeenv(env);
        }
        
        void GurobiLpSolver::setGurobiEnvironmentProperties() const {
			int error = 0;

			// Enable the following line to only print the output of Gurobi if the debug flag is set.
            error = GRBsetintparam(env, "OutputFlag", storm::settings::debugSettings().isDebugSet() || storm::settings::gurobiSettings().isOutputSet() ? 1 : 0);
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to set Gurobi Parameter OutputFlag (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            
            // Enable the following line to restrict Gurobi to one thread only.
            error = GRBsetintparam(env, "Threads", storm::settings::gurobiSettings().getNumberOfThreads());
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to set Gurobi Parameter Threads (" << GRBgeterrormsg(env) << ", error code " << error << ").");

            // Enable the following line to force Gurobi to be as precise about the binary variables as required by the given precision option.
            error = GRBsetdblparam(env, "IntFeasTol", storm::settings::gurobiSettings().getIntegerTolerance());
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to set Gurobi Parameter IntFeasTol (" << GRBgeterrormsg(env) << ", error code " << error << ").");
        }
        
        void GurobiLpSolver::update() const {
            int error = GRBupdatemodel(model);
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to update Gurobi model (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            
            // Since the model changed, we erase the optimality flag.
            this->currentModelHasBeenOptimized = false;
        }
        
        storm::expressions::Variable GurobiLpSolver::addBoundedContinuousVariable(std::string const& name, double lowerBound, double upperBound, double objectiveFunctionCoefficient) {
            storm::expressions::Variable newVariable = manager->declareVariable(name, manager->getRationalType());
            this->addVariable(newVariable, GRB_CONTINUOUS, lowerBound, upperBound, objectiveFunctionCoefficient);
            return newVariable;
        }
        
        storm::expressions::Variable GurobiLpSolver::addLowerBoundedContinuousVariable(std::string const& name, double lowerBound, double objectiveFunctionCoefficient) {
            storm::expressions::Variable newVariable = manager->declareVariable(name, manager->getRationalType());
            this->addVariable(newVariable, GRB_CONTINUOUS, lowerBound, GRB_INFINITY, objectiveFunctionCoefficient);
            return newVariable;
        }

        storm::expressions::Variable GurobiLpSolver::addUpperBoundedContinuousVariable(std::string const& name, double upperBound, double objectiveFunctionCoefficient) {
            storm::expressions::Variable newVariable = manager->declareVariable(name, manager->getRationalType());
            this->addVariable(newVariable, GRB_CONTINUOUS, -GRB_INFINITY, upperBound, objectiveFunctionCoefficient);
            return newVariable;
        }

        storm::expressions::Variable GurobiLpSolver::addUnboundedContinuousVariable(std::string const& name, double objectiveFunctionCoefficient) {
            storm::expressions::Variable newVariable = manager->declareVariable(name, manager->getRationalType());
            this->addVariable(newVariable, GRB_CONTINUOUS, -GRB_INFINITY, GRB_INFINITY, objectiveFunctionCoefficient);
            return newVariable;
        }
        
        storm::expressions::Variable GurobiLpSolver::addBoundedIntegerVariable(std::string const& name, double lowerBound, double upperBound, double objectiveFunctionCoefficient) {
            storm::expressions::Variable newVariable = manager->declareVariable(name, manager->getIntegerType());
            this->addVariable(newVariable, GRB_INTEGER, lowerBound, upperBound, objectiveFunctionCoefficient);
            return newVariable;
        }
        
        storm::expressions::Variable GurobiLpSolver::addLowerBoundedIntegerVariable(std::string const& name, double lowerBound, double objectiveFunctionCoefficient) {
            storm::expressions::Variable newVariable = manager->declareVariable(name, manager->getIntegerType());
            this->addVariable(newVariable, GRB_INTEGER, lowerBound, GRB_INFINITY, objectiveFunctionCoefficient);
            return newVariable;
        }
        
        storm::expressions::Variable GurobiLpSolver::addUpperBoundedIntegerVariable(std::string const& name, double upperBound, double objectiveFunctionCoefficient) {
            storm::expressions::Variable newVariable = manager->declareVariable(name, manager->getIntegerType());
            this->addVariable(newVariable, GRB_INTEGER, -GRB_INFINITY, upperBound, objectiveFunctionCoefficient);
            return newVariable;
        }
        
        storm::expressions::Variable GurobiLpSolver::addUnboundedIntegerVariable(std::string const& name, double objectiveFunctionCoefficient) {
            storm::expressions::Variable newVariable = manager->declareVariable(name, manager->getIntegerType());
            this->addVariable(newVariable, GRB_INTEGER, -GRB_INFINITY, GRB_INFINITY, objectiveFunctionCoefficient);
            return newVariable;
        }
        
        storm::expressions::Variable GurobiLpSolver::addBinaryVariable(std::string const& name, double objectiveFunctionCoefficient) {
            storm::expressions::Variable newVariable = manager->declareVariable(name, manager->getIntegerType());
            this->addVariable(newVariable, GRB_BINARY, 0, 1, objectiveFunctionCoefficient);
            return newVariable;
        }

        void GurobiLpSolver::addVariable(storm::expressions::Variable const& variable, char variableType, double lowerBound, double upperBound, double objectiveFunctionCoefficient) {
            // Check for valid variable type.
            STORM_LOG_ASSERT(variableType == GRB_CONTINUOUS || variableType == GRB_INTEGER || variableType == GRB_BINARY, "Illegal type '" << variableType << "' for Gurobi variable.");
            
            // Finally, create the actual variable.
            int error = 0;
            error = GRBaddvar(model, 0, nullptr, nullptr, objectiveFunctionCoefficient, lowerBound, upperBound, variableType, variable.getName().c_str());
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Could not create binary Gurobi variable (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            this->variableToIndexMap.emplace(variable, nextVariableIndex);
            ++nextVariableIndex;
        }
                
        void GurobiLpSolver::addConstraint(std::string const& name, storm::expressions::Expression const& constraint) {            
            STORM_LOG_THROW(constraint.isRelationalExpression(), storm::exceptions::InvalidArgumentException, "Illegal constraint is not a relational expression.");
            STORM_LOG_THROW(constraint.getOperator() != storm::expressions::OperatorType::NotEqual, storm::exceptions::InvalidArgumentException, "Illegal constraint uses inequality operator.");
            
            storm::expressions::LinearCoefficientVisitor::VariableCoefficients leftCoefficients = storm::expressions::LinearCoefficientVisitor().getLinearCoefficients(constraint.getOperand(0));
            storm::expressions::LinearCoefficientVisitor::VariableCoefficients rightCoefficients = storm::expressions::LinearCoefficientVisitor().getLinearCoefficients(constraint.getOperand(1));
            leftCoefficients.separateVariablesFromConstantPart(rightCoefficients);
            
            // Now we need to transform the coefficients to the vector representation.
            std::vector<int> variables;
            std::vector<double> coefficients;
            for (auto const& variableCoefficientPair : leftCoefficients) {
                auto variableIndexPair = this->variableToIndexMap.find(variableCoefficientPair.first);
                variables.push_back(variableIndexPair->second);
                coefficients.push_back(leftCoefficients.getCoefficient(variableIndexPair->first));
            }
            
            // Determine the type of the constraint and add it properly.
            int error = 0;
            switch (constraint.getOperator()) {
                case storm::expressions::OperatorType::Less:
                    error = GRBaddconstr(model, variables.size(), variables.data(), coefficients.data(), GRB_LESS_EQUAL, rightCoefficients.getConstantPart() - storm::settings::gurobiSettings().getIntegerTolerance(), name == "" ? nullptr : name.c_str());
                    break;
                case storm::expressions::OperatorType::LessOrEqual:
                    error = GRBaddconstr(model, variables.size(), variables.data(), coefficients.data(), GRB_LESS_EQUAL, rightCoefficients.getConstantPart(), name == "" ? nullptr : name.c_str());
                    break;
                case storm::expressions::OperatorType::Greater:
                    error = GRBaddconstr(model, variables.size(), variables.data(), coefficients.data(), GRB_GREATER_EQUAL, rightCoefficients.getConstantPart() + storm::settings::gurobiSettings().getIntegerTolerance(), name == "" ? nullptr : name.c_str());
                    break;
                case storm::expressions::OperatorType::GreaterOrEqual:
                    error = GRBaddconstr(model, variables.size(), variables.data(), coefficients.data(), GRB_GREATER_EQUAL, rightCoefficients.getConstantPart(), name == "" ? nullptr : name.c_str());
                    break;
                case storm::expressions::OperatorType::Equal:
                    error = GRBaddconstr(model, variables.size(), variables.data(), coefficients.data(), GRB_EQUAL, rightCoefficients.getConstantPart(), name == "" ? nullptr : name.c_str());
                    break;
                default:
                    STORM_LOG_ASSERT(false, "Illegal operator in LP solver constraint.");
            }
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Could not assert constraint (" << GRBgeterrormsg(env) << ", error code " << error << ").");
        }
        
        void GurobiLpSolver::optimize() const {
            // First incorporate all recent changes.
            this->update();
         
            // Set the most recently set model sense.
            int error = GRBsetintattr(model, "ModelSense", this->getModelSense() == ModelSense::Minimize ? 1 : -1);
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to set Gurobi model sense (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            
            // Then we actually optimize the model.
            error = GRBoptimize(model);
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to optimize Gurobi model (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            
            this->currentModelHasBeenOptimized = true;
        }
        
        bool GurobiLpSolver::isInfeasible() const {
            if (!this->currentModelHasBeenOptimized) {
                throw storm::exceptions::InvalidStateException() << "Illegal call to GurobiLpSolver::isInfeasible: model has not been optimized.";
            }

            int optimalityStatus = 0;
            
            int error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimalityStatus);
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to retrieve optimization status of Gurobi model (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            
            // By default, Gurobi may tell us only that the model is either infeasible or unbounded. To decide which one
            // it is, we need to perform an extra step.
            if (optimalityStatus == GRB_INF_OR_UNBD) {
                error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_DUALREDUCTIONS, 0);
                STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to set Gurobi parameter (" << GRBgeterrormsg(env) << ", error code " << error << ").");
                
                this->optimize();
                
                error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimalityStatus);
                STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to retrieve optimization status of Gurobi model (" << GRBgeterrormsg(env) << ", error code " << error << ").");
                
                error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_DUALREDUCTIONS, 1);
                STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to set Gurobi parameter (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            }
            
            return optimalityStatus == GRB_INFEASIBLE;
        }
        
        bool GurobiLpSolver::isUnbounded() const {
            if (!this->currentModelHasBeenOptimized) {
                throw storm::exceptions::InvalidStateException() << "Illegal call to GurobiLpSolver::isUnbounded: model has not been optimized.";
            }

            int optimalityStatus = 0;
            
            int error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimalityStatus);
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to retrieve optimization status of Gurobi model (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            
            // By default, Gurobi may tell us only that the model is either infeasible or unbounded. To decide which one
            // it is, we need to perform an extra step.
            if (optimalityStatus == GRB_INF_OR_UNBD) {
                error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_DUALREDUCTIONS, 0);
                STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to set Gurobi parameter (" << GRBgeterrormsg(env) << ", error code " << error << ").");
                
                this->optimize();

                error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimalityStatus);
                STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to retrieve optimization status of Gurobi model (" << GRBgeterrormsg(env) << ", error code " << error << ").");

                error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_DUALREDUCTIONS, 1);
                STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to set Gurobi parameter (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            }
            
            return optimalityStatus == GRB_UNBOUNDED;
        }
        
        bool GurobiLpSolver::isOptimal() const {
            if (!this->currentModelHasBeenOptimized) {
                return false;
            }
            int optimalityStatus = 0;
            
            int error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimalityStatus);
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to retrieve optimization status of Gurobi model (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            
            return optimalityStatus == GRB_OPTIMAL;
        }
        
        double GurobiLpSolver::getContinuousValue(storm::expressions::Variable const& variable) const {
            if (!this->isOptimal()) {
                STORM_LOG_THROW(!this->isInfeasible(), storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from infeasible model (" << GRBgeterrormsg(env) << ").");
                STORM_LOG_THROW(!this->isUnbounded(), storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from unbounded model (" << GRBgeterrormsg(env) << ").");
                STORM_LOG_THROW(false, storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from unoptimized model (" << GRBgeterrormsg(env) << ").");
            }
            
            auto variableIndexPair = this->variableToIndexMap.find(variable);
            STORM_LOG_THROW(variableIndexPair != this->variableToIndexMap.end(), storm::exceptions::InvalidAccessException, "Accessing value of unknown variable '" << variable.getName() << "'.");
            
            double value = 0;
            int error = GRBgetdblattrelement(model, GRB_DBL_ATTR_X, variableIndexPair->second, &value);
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to get Gurobi solution (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            
            return value;
        }
        
        int_fast64_t GurobiLpSolver::getIntegerValue(storm::expressions::Variable const& variable) const {
            if (!this->isOptimal()) {
                STORM_LOG_THROW(!this->isInfeasible(), storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from infeasible model (" << GRBgeterrormsg(env) << ").");
                STORM_LOG_THROW(!this->isUnbounded(), storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from unbounded model (" << GRBgeterrormsg(env) << ").");
                STORM_LOG_THROW(false, storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from unoptimized model (" << GRBgeterrormsg(env) << ").");
            }
            
            auto variableIndexPair = this->variableToIndexMap.find(variable);
            STORM_LOG_THROW(variableIndexPair != this->variableToIndexMap.end(), storm::exceptions::InvalidAccessException, "Accessing value of unknown variable '" << variable.getName() << "'.");
            
            double value = 0;
            int error = GRBgetdblattrelement(model, GRB_DBL_ATTR_X, variableIndexPair->second, &value);
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to get Gurobi solution (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            STORM_LOG_THROW(std::abs(static_cast<int>(value) - value) <= storm::settings::gurobiSettings().getIntegerTolerance(), storm::exceptions::InvalidStateException, "Illegal value for integer variable in Gurobi solution (" << value << ").");
            
            return static_cast<int_fast64_t>(value);
        }
        
        bool GurobiLpSolver::getBinaryValue(storm::expressions::Variable const& variable) const {
            if (!this->isOptimal()) {
                STORM_LOG_THROW(!this->isInfeasible(), storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from infeasible model (" << GRBgeterrormsg(env) << ").");
                STORM_LOG_THROW(!this->isUnbounded(), storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from unbounded model (" << GRBgeterrormsg(env) << ").");
                STORM_LOG_THROW(false, storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from unoptimized model (" << GRBgeterrormsg(env) << ").");
            }

            auto variableIndexPair = this->variableToIndexMap.find(variable);
            STORM_LOG_THROW(variableIndexPair != this->variableToIndexMap.end(), storm::exceptions::InvalidAccessException, "Accessing value of unknown variable '" << variable.getName() << "'.");
            
            double value = 0;
            int error = GRBgetdblattrelement(model, GRB_DBL_ATTR_X, variableIndexPair->second, &value);
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to get Gurobi solution (" << GRBgeterrormsg(env) << ", error code " << error << ").");

            if (value > 0.5) {
                STORM_LOG_THROW(std::abs(static_cast<int>(value) - 1) <= storm::settings::gurobiSettings().getIntegerTolerance(), storm::exceptions::InvalidStateException, "Illegal value for integer variable in Gurobi solution (" << value << ").");
            } else {
                STORM_LOG_THROW(value <= storm::settings::gurobiSettings().getIntegerTolerance(), storm::exceptions::InvalidStateException, "Illegal value for integer variable in Gurobi solution (" << value << ").");
            }
            
            return static_cast<bool>(value);
        }
        
        double GurobiLpSolver::getObjectiveValue() const {
            if (!this->isOptimal()) {
                STORM_LOG_THROW(!this->isInfeasible(), storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from infeasible model (" << GRBgeterrormsg(env) << ").");
                STORM_LOG_THROW(!this->isUnbounded(), storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from unbounded model (" << GRBgeterrormsg(env) << ").");
                STORM_LOG_THROW(false, storm::exceptions::InvalidAccessException, "Unable to get Gurobi solution from unoptimized model (" << GRBgeterrormsg(env) << ").");
            }
            
            double value = 0;
            int error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &value);
            STORM_LOG_THROW(error == 0, storm::exceptions::InvalidStateException, "Unable to get Gurobi solution (" << GRBgeterrormsg(env) << ", error code " << error << ").");
            
            return value;
        }
        
        void GurobiLpSolver::writeModelToFile(std::string const& filename) const {
            int error = GRBwrite(model, filename.c_str());
            if (error) {
				LOG4CPLUS_ERROR(logger, "Unable to write Gurobi model (" << GRBgeterrormsg(env) << ", error code " << error << ") to file.");
				throw storm::exceptions::InvalidStateException() << "Unable to write Gurobi model (" << GRBgeterrormsg(env) << ", error code " << error << ") to file.";
            }
        }
    }
}

#endif