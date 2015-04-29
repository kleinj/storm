#include "gtest/gtest.h"
#include "storm-config.h"

#include "src/logic/Formulas.h"
#include "src/utility/solver.h"
#include "src/modelchecker/prctl/HybridDtmcPrctlModelChecker.h"
#include "src/modelchecker/results/HybridQuantitativeCheckResult.h"
#include "src/modelchecker/results/SymbolicQualitativeCheckResult.h"
#include "src/modelchecker/results/SymbolicQuantitativeCheckResult.h"
#include "src/parser/PrismParser.h"
#include "src/builder/DdPrismModelBuilder.h"
#include "src/models/symbolic/Dtmc.h"
#include "src/settings/SettingsManager.h"

TEST(GmmxxHybridDtmcPrctlModelCheckerTest, Die) {
    storm::prism::Program program = storm::parser::PrismParser::parse(STORM_CPP_TESTS_BASE_PATH "/functional/builder/die.pm");

    // Build the die model with its reward model.
    typename storm::builder::DdPrismModelBuilder<storm::dd::DdType::CUDD>::Options options;
    options.buildRewards = true;
    options.rewardModelName = "coin_flips";
    std::shared_ptr<storm::models::symbolic::Model<storm::dd::DdType::CUDD>> model = storm::builder::DdPrismModelBuilder<storm::dd::DdType::CUDD>::translateProgram(program, options);
    EXPECT_EQ(13, model->getNumberOfStates());
    EXPECT_EQ(20, model->getNumberOfTransitions());

    ASSERT_EQ(model->getType(), storm::models::ModelType::Dtmc);
    
    std::shared_ptr<storm::models::symbolic::Dtmc<storm::dd::DdType::CUDD>> dtmc = model->as<storm::models::symbolic::Dtmc<storm::dd::DdType::CUDD>>();
    
    storm::modelchecker::HybridDtmcPrctlModelChecker<storm::dd::DdType::CUDD, double> checker(*dtmc, std::unique_ptr<storm::utility::solver::LinearEquationSolverFactory<double>>(new storm::utility::solver::GmmxxLinearEquationSolverFactory<double>()));
    
    auto labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("one");
    auto eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    std::unique_ptr<storm::modelchecker::CheckResult> result = checker.check(*eventuallyFormula);
    result->filter(storm::modelchecker::SymbolicQualitativeCheckResult<storm::dd::DdType::CUDD>(model->getReachableStates(), model->getInitialStates()));
    storm::modelchecker::HybridQuantitativeCheckResult<storm::dd::DdType::CUDD>& quantitativeResult1 = result->asHybridQuantitativeCheckResult<storm::dd::DdType::CUDD>();
    
    EXPECT_NEAR(1.0/6.0, quantitativeResult1.getMin(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    EXPECT_NEAR(1.0/6.0, quantitativeResult1.getMax(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("two");
    eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    result = checker.check(*eventuallyFormula);
    result->filter(storm::modelchecker::SymbolicQualitativeCheckResult<storm::dd::DdType::CUDD>(model->getReachableStates(), model->getInitialStates()));
    storm::modelchecker::HybridQuantitativeCheckResult<storm::dd::DdType::CUDD>& quantitativeResult2 = result->asHybridQuantitativeCheckResult<storm::dd::DdType::CUDD>();
    
    EXPECT_NEAR(1.0/6.0, quantitativeResult2.getMin(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    EXPECT_NEAR(1.0/6.0, quantitativeResult2.getMax(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("three");
    eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    result = checker.check(*eventuallyFormula);
    result->filter(storm::modelchecker::SymbolicQualitativeCheckResult<storm::dd::DdType::CUDD>(model->getReachableStates(), model->getInitialStates()));
    storm::modelchecker::HybridQuantitativeCheckResult<storm::dd::DdType::CUDD>& quantitativeResult3 = result->asHybridQuantitativeCheckResult<storm::dd::DdType::CUDD>();
    
    EXPECT_NEAR(1.0/6.0, quantitativeResult3.getMin(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    EXPECT_NEAR(1.0/6.0, quantitativeResult3.getMax(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    
    auto done = std::make_shared<storm::logic::AtomicLabelFormula>("done");
    auto reachabilityRewardFormula = std::make_shared<storm::logic::ReachabilityRewardFormula>(done);
    
    result = checker.check(*reachabilityRewardFormula);
    result->filter(storm::modelchecker::SymbolicQualitativeCheckResult<storm::dd::DdType::CUDD>(model->getReachableStates(), model->getInitialStates()));
    storm::modelchecker::HybridQuantitativeCheckResult<storm::dd::DdType::CUDD>& quantitativeResult4 = result->asHybridQuantitativeCheckResult<storm::dd::DdType::CUDD>();
    
    EXPECT_NEAR(11.0/3.0, quantitativeResult4.getMin(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    EXPECT_NEAR(11.0/3.0, quantitativeResult4.getMax(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
}

TEST(GmmxxHybridDtmcPrctlModelCheckerTest, Crowds) {
    storm::prism::Program program = storm::parser::PrismParser::parse(STORM_CPP_TESTS_BASE_PATH "/functional/builder/crowds-5-5.pm");
    
    std::shared_ptr<storm::models::symbolic::Model<storm::dd::DdType::CUDD>> model = storm::builder::DdPrismModelBuilder<storm::dd::DdType::CUDD>::translateProgram(program);
    EXPECT_EQ(8607, model->getNumberOfStates());
    EXPECT_EQ(15113, model->getNumberOfTransitions());
    
    ASSERT_EQ(model->getType(), storm::models::ModelType::Dtmc);
    
    std::shared_ptr<storm::models::symbolic::Dtmc<storm::dd::DdType::CUDD>> dtmc = model->as<storm::models::symbolic::Dtmc<storm::dd::DdType::CUDD>>();
    
    storm::modelchecker::HybridDtmcPrctlModelChecker<storm::dd::DdType::CUDD, double> checker(*dtmc, std::unique_ptr<storm::utility::solver::LinearEquationSolverFactory<double>>(new storm::utility::solver::GmmxxLinearEquationSolverFactory<double>()));
    
    auto labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("observe0Greater1");
    auto eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    std::unique_ptr<storm::modelchecker::CheckResult> result = checker.check(*eventuallyFormula);
    result->filter(storm::modelchecker::SymbolicQualitativeCheckResult<storm::dd::DdType::CUDD>(model->getReachableStates(), model->getInitialStates()));
    storm::modelchecker::HybridQuantitativeCheckResult<storm::dd::DdType::CUDD>& quantitativeResult1 = result->asHybridQuantitativeCheckResult<storm::dd::DdType::CUDD>();
    
    EXPECT_NEAR(0.3328800375801578281, quantitativeResult1.getMin(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    EXPECT_NEAR(0.3328800375801578281, quantitativeResult1.getMax(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("observeIGreater1");
    eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    result = checker.check(*eventuallyFormula);
    result->filter(storm::modelchecker::SymbolicQualitativeCheckResult<storm::dd::DdType::CUDD>(model->getReachableStates(), model->getInitialStates()));
    storm::modelchecker::HybridQuantitativeCheckResult<storm::dd::DdType::CUDD>& quantitativeResult2 = result->asHybridQuantitativeCheckResult<storm::dd::DdType::CUDD>();
    
    EXPECT_NEAR(0.1522194965, quantitativeResult2.getMin(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    EXPECT_NEAR(0.1522194965, quantitativeResult2.getMax(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("observeOnlyTrueSender");
    eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    result = checker.check(*eventuallyFormula);
    result->filter(storm::modelchecker::SymbolicQualitativeCheckResult<storm::dd::DdType::CUDD>(model->getReachableStates(), model->getInitialStates()));
    storm::modelchecker::HybridQuantitativeCheckResult<storm::dd::DdType::CUDD>& quantitativeResult3 = result->asHybridQuantitativeCheckResult<storm::dd::DdType::CUDD>();
    
    EXPECT_NEAR(0.32153724292835045, quantitativeResult3.getMin(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    EXPECT_NEAR(0.32153724292835045, quantitativeResult3.getMax(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
}

TEST(GmmxxHybridDtmcPrctlModelCheckerTest, SynchronousLeader) {
    storm::prism::Program program = storm::parser::PrismParser::parse(STORM_CPP_TESTS_BASE_PATH "/functional/builder/leader-3-5.pm");
    
    // Build the die model with its reward model.
    typename storm::builder::DdPrismModelBuilder<storm::dd::DdType::CUDD>::Options options;
    options.buildRewards = true;
    options.rewardModelName = "num_rounds";
    std::shared_ptr<storm::models::symbolic::Model<storm::dd::DdType::CUDD>> model = storm::builder::DdPrismModelBuilder<storm::dd::DdType::CUDD>::translateProgram(program, options);
    EXPECT_EQ(273, model->getNumberOfStates());
    EXPECT_EQ(397, model->getNumberOfTransitions());
    
    ASSERT_EQ(model->getType(), storm::models::ModelType::Dtmc);
    
    std::shared_ptr<storm::models::symbolic::Dtmc<storm::dd::DdType::CUDD>> dtmc = model->as<storm::models::symbolic::Dtmc<storm::dd::DdType::CUDD>>();
    
    storm::modelchecker::HybridDtmcPrctlModelChecker<storm::dd::DdType::CUDD, double> checker(*dtmc, std::unique_ptr<storm::utility::solver::LinearEquationSolverFactory<double>>(new storm::utility::solver::GmmxxLinearEquationSolverFactory<double>()));
    
    auto labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("elected");
    auto eventuallyFormula = std::make_shared<storm::logic::EventuallyFormula>(labelFormula);
    
    std::unique_ptr<storm::modelchecker::CheckResult> result = checker.check(*eventuallyFormula);
    result->filter(storm::modelchecker::SymbolicQualitativeCheckResult<storm::dd::DdType::CUDD>(model->getReachableStates(), model->getInitialStates()));
    storm::modelchecker::SymbolicQuantitativeCheckResult<storm::dd::DdType::CUDD>& quantitativeResult1 = result->asSymbolicQuantitativeCheckResult<storm::dd::DdType::CUDD>();
    
    EXPECT_NEAR(1.0, quantitativeResult1.getMin(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    EXPECT_NEAR(1.0, quantitativeResult1.getMax(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("elected");
    auto trueFormula = std::make_shared<storm::logic::BooleanLiteralFormula>(true);
    auto boundedUntilFormula = std::make_shared<storm::logic::BoundedUntilFormula>(trueFormula, labelFormula, 20);
    
    result = checker.check(*boundedUntilFormula);
    result->filter(storm::modelchecker::SymbolicQualitativeCheckResult<storm::dd::DdType::CUDD>(model->getReachableStates(), model->getInitialStates()));
    storm::modelchecker::HybridQuantitativeCheckResult<storm::dd::DdType::CUDD>& quantitativeResult2 = result->asHybridQuantitativeCheckResult<storm::dd::DdType::CUDD>();
    
    EXPECT_NEAR(0.99999989760000074, quantitativeResult2.getMin(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    EXPECT_NEAR(0.99999989760000074, quantitativeResult2.getMax(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    
    labelFormula = std::make_shared<storm::logic::AtomicLabelFormula>("elected");
    auto reachabilityRewardFormula = std::make_shared<storm::logic::ReachabilityRewardFormula>(labelFormula);
    
    result = checker.check(*reachabilityRewardFormula);
    result->filter(storm::modelchecker::SymbolicQualitativeCheckResult<storm::dd::DdType::CUDD>(model->getReachableStates(), model->getInitialStates()));
    storm::modelchecker::HybridQuantitativeCheckResult<storm::dd::DdType::CUDD>& quantitativeResult3 = result->asHybridQuantitativeCheckResult<storm::dd::DdType::CUDD>();
    
    EXPECT_NEAR(1.0416666666666643, quantitativeResult3.getMin(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
    EXPECT_NEAR(1.0416666666666643, quantitativeResult3.getMax(), storm::settings::gmmxxEquationSolverSettings().getPrecision());
}
