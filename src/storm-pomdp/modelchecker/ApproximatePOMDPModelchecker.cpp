#include "ApproximatePOMDPModelchecker.h"

#include <boost/algorithm/string.hpp>


#include "storm/utility/ConstantsComparator.h"
#include "storm/models/sparse/Dtmc.h"
#include "storm/modelchecker/prctl/SparseDtmcPrctlModelChecker.h"
#include "storm/utility/vector.h"
#include "storm/modelchecker/results/CheckResult.h"
#include "storm/modelchecker/results/ExplicitQualitativeCheckResult.h"
#include "storm/modelchecker/results/ExplicitQuantitativeCheckResult.h"
#include "storm/api/properties.h"
#include "storm-parsers/api/storm-parsers.h"

namespace storm {
    namespace pomdp {
        namespace modelchecker {
            template<typename ValueType, typename RewardModelType>
            ApproximatePOMDPModelchecker<ValueType, RewardModelType>::ApproximatePOMDPModelchecker() {
                //Intentionally left empty
            }

            template<typename ValueType, typename RewardModelType>
            /*std::unique_ptr<POMDPCheckResult>*/ void
            ApproximatePOMDPModelchecker<ValueType, RewardModelType>::computeReachabilityProbability(
                    storm::models::sparse::Pomdp<ValueType, RewardModelType> const &pomdp,
                    std::set<uint32_t> targetObservations, bool min, uint64_t gridResolution) {
                //TODO add timing
                uint64_t maxIterations = 100;
                bool finished = false;
                uint64_t iteration = 0;

                std::vector<storm::pomdp::Belief<ValueType>> beliefList;
                std::vector<bool> beliefIsTarget;
                uint64_t nextId = 0;
                // Initial belief always has ID 0
                storm::pomdp::Belief<ValueType> initialBelief = getInitialBelief(pomdp, nextId);
                ++nextId;
                beliefList.push_back(initialBelief);
                beliefIsTarget.push_back(
                        targetObservations.find(initialBelief.observation) != targetObservations.end());


                std::vector<storm::pomdp::Belief<ValueType>> beliefGrid;
                constructBeliefGrid(pomdp, targetObservations, gridResolution, beliefList, beliefGrid, beliefIsTarget,
                                    nextId);
                nextId = beliefList.size();

                // ID -> Value
                std::map<uint64_t, ValueType> result;
                std::map<uint64_t, ValueType> result_backup;
                // ID -> ActionIndex
                std::map<uint64_t, uint64_t> chosenActions;

                // ID -> Observation -> Probability
                std::map<uint64_t, std::vector<std::map<uint32_t, ValueType>>> observationProbabilities;
                // current ID -> action -> next ID
                std::map<uint64_t, std::vector<std::map<uint32_t, uint64_t>>> nextBelieves;

                for (size_t i = 0; i < beliefGrid.size(); ++i) {
                    auto currentBelief = beliefGrid[i];
                    bool isTarget = beliefIsTarget[currentBelief.id];
                    if (isTarget) {
                        result.emplace(std::make_pair(currentBelief.id, storm::utility::one<ValueType>()));
                        result_backup.emplace(std::make_pair(currentBelief.id, storm::utility::one<ValueType>()));
                    } else {
                        result.emplace(std::make_pair(currentBelief.id, storm::utility::zero<ValueType>()));
                        result_backup.emplace(std::make_pair(currentBelief.id, storm::utility::zero<ValueType>()));

                        std::vector<std::map<uint32_t, ValueType>> observationProbabilitiesInAction;
                        std::vector<std::map<uint32_t, uint64_t>> nextBelievesInAction;

                        uint64_t numChoices = pomdp.getNumberOfChoices(
                                pomdp.getStatesWithObservation(currentBelief.observation).front());
                        for (uint64_t action = 0; action < numChoices; ++action) {
                            std::map<uint32_t, ValueType> actionObservationProbabilities = computeObservationProbabilitiesAfterAction(
                                    pomdp, currentBelief, action);
                            std::map<uint32_t, uint64_t> actionObservationBelieves;
                            for (auto iter = actionObservationProbabilities.begin();
                                 iter != actionObservationProbabilities.end(); ++iter) {
                                uint32_t observation = iter->first;
                                actionObservationBelieves[observation] = getBeliefAfterActionAndObservation(pomdp,
                                                                                                            beliefList,
                                                                                                            beliefIsTarget,
                                                                                                            targetObservations,
                                                                                                            currentBelief,
                                                                                                            action,
                                                                                                            observation,
                                                                                                            nextId);
                                nextId = beliefList.size();
                            }
                            observationProbabilitiesInAction.push_back(actionObservationProbabilities);
                            nextBelievesInAction.push_back(actionObservationBelieves);
                        }
                        observationProbabilities.emplace(
                                std::make_pair(currentBelief.id, observationProbabilitiesInAction));
                        nextBelieves.emplace(std::make_pair(currentBelief.id, nextBelievesInAction));
                    }
                }
                // Value Iteration
                auto cc = storm::utility::ConstantsComparator<ValueType>();
                while (!finished && iteration < maxIterations) {
                    STORM_LOG_DEBUG("Iteration " << iteration + 1);
                    bool improvement = false;
                    for (size_t i = 0; i < beliefGrid.size(); ++i) {
                        storm::pomdp::Belief<ValueType> currentBelief = beliefGrid[i];
                        bool isTarget = beliefIsTarget[currentBelief.id];
                        if (!isTarget) {
                            // we can take any state with the observation as they have the same number of choices
                            uint64_t numChoices = pomdp.getNumberOfChoices(
                                    pomdp.getStatesWithObservation(currentBelief.observation).front());
                            // Initialize the values for the value iteration
                            ValueType chosenValue = min ? storm::utility::infinity<ValueType>()
                                                        : -storm::utility::infinity<ValueType>();
                            uint64_t chosenActionIndex = std::numeric_limits<uint64_t>::infinity();
                            ValueType currentValue;

                            for (uint64_t action = 0; action < numChoices; ++action) {
                                currentValue = storm::utility::zero<ValueType>(); // simply change this for rewards?
                                for (auto iter = observationProbabilities[currentBelief.id][action].begin();
                                     iter != observationProbabilities[currentBelief.id][action].end(); ++iter) {
                                    uint32_t observation = iter->first;
                                    storm::pomdp::Belief<ValueType> nextBelief = beliefList[nextBelieves[currentBelief.id][action][observation]];
                                    // compute subsimplex and lambdas according to the Lovejoy paper to approximate the next belief
                                    std::pair<std::vector<std::vector<ValueType>>, std::vector<ValueType>> temp = computeSubSimplexAndLambdas(
                                            nextBelief.probabilities, gridResolution);
                                    std::vector<std::vector<ValueType>> subSimplex = temp.first;
                                    std::vector<ValueType> lambdas = temp.second;

                                    auto sum = storm::utility::zero<ValueType>();
                                    for (size_t j = 0; j < lambdas.size(); ++j) {
                                        if (lambdas[j] != storm::utility::zero<ValueType>()) {
                                            sum += lambdas[j] * result_backup.at(
                                                    getBeliefIdInVector(beliefGrid, observation, subSimplex[j]));
                                        }
                                    }
                                    currentValue += iter->second * sum;
                                }

                                // Update the selected actions TODO make this nicer
                                if ((min && cc.isLess(storm::utility::zero<ValueType>(), chosenValue - currentValue)) ||
                                    (!min &&
                                     cc.isLess(storm::utility::zero<ValueType>(), currentValue - chosenValue))) {
                                    chosenValue = currentValue;
                                    chosenActionIndex = action;
                                } else if ((min && cc.isEqual(storm::utility::zero<ValueType>(),
                                                              chosenValue - currentValue)) ||
                                           (!min &&
                                            cc.isEqual(storm::utility::zero<ValueType>(),
                                                       currentValue - chosenValue))) {
                                    chosenValue = currentValue;
                                    chosenActionIndex = action;
                                }
                                // TODO tie breaker?
                            }
                            result[currentBelief.id] = chosenValue;
                            chosenActions[currentBelief.id] = chosenActionIndex;
                            // Check if the iteration brought an improvement
                            if ((min && cc.isLess(storm::utility::zero<ValueType>(),
                                                  result_backup[currentBelief.id] - result[currentBelief.id])) ||
                                (!min && cc.isLess(storm::utility::zero<ValueType>(),
                                                   result[currentBelief.id] - result_backup[currentBelief.id]))) {
                                improvement = true;
                            }
                        }
                    }
                    finished = !improvement;
                    // back up
                    for (auto iter = result.begin(); iter != result.end(); ++iter) {
                        result_backup[iter->first] = result[iter->first];
                    }
                    ++iteration;
                }

                STORM_PRINT("Grid approximation took " << iteration << " iterations" << std::endl);

                beliefGrid.push_back(initialBelief);
                beliefIsTarget.push_back(
                        targetObservations.find(initialBelief.observation) != targetObservations.end());

                std::pair<std::vector<std::vector<ValueType>>, std::vector<ValueType>> temp = computeSubSimplexAndLambdas(
                        initialBelief.probabilities, gridResolution);
                std::vector<std::vector<ValueType>> initSubSimplex = temp.first;
                std::vector<ValueType> initLambdas = temp.second;

                auto overApprox = storm::utility::zero<ValueType>();
                for (size_t j = 0; j < initLambdas.size(); ++j) {
                    if (initLambdas[j] != storm::utility::zero<ValueType>()) {
                        overApprox += initLambdas[j] *
                                      result_backup[getBeliefIdInVector(beliefGrid, initialBelief.observation,
                                                                        initSubSimplex[j])];
                    }
                }

                // Now onto the under-approximation

                std::set<uint64_t> visitedBelieves;
                std::deque<uint64_t> believesToBeExpanded;
                std::map<uint64_t, uint64_t> beliefStateMap;
                std::vector<std::map<uint64_t, ValueType>> transitions;
                std::vector<uint64_t> targetStates;

                uint64_t stateId = 0;
                beliefStateMap[initialBelief.id] = stateId;
                ++stateId;

                // Expand the believes TODO capsuling
                visitedBelieves.insert(initialBelief.id);
                believesToBeExpanded.push_back(initialBelief.id);
                while (!believesToBeExpanded.empty()) {
                    auto currentBeliefId = believesToBeExpanded.front();
                    std::map<uint64_t, ValueType> transitionsInState;
                    STORM_LOG_DEBUG("Exploring Belief " << beliefList[currentBeliefId].observation << "||"
                                                        << beliefList[currentBeliefId].probabilities);
                    if (beliefIsTarget[currentBeliefId]) {
                        // add a self-loop to target states and save them
                        transitionsInState[beliefStateMap[currentBeliefId]] = storm::utility::one<ValueType>();
                        targetStates.push_back(beliefStateMap[currentBeliefId]);
                    } else {
                        if (chosenActions.find(currentBeliefId) == chosenActions.end()) {
                            // If the current Belief is not part of the grid, we have not computed the action to choose yet
                            chosenActions[currentBeliefId] = extractBestAction(pomdp, beliefList, beliefIsTarget,
                                                                               targetObservations,
                                                                               observationProbabilities,
                                                                               nextBelieves, result, gridResolution,
                                                                               currentBeliefId, beliefList.size(), min);
                        }
                        for (auto iter = observationProbabilities[currentBeliefId][chosenActions[currentBeliefId]].begin();
                             iter !=
                             observationProbabilities[currentBeliefId][chosenActions[currentBeliefId]].end(); ++iter) {
                            uint32_t observation = iter->first;
                            uint64_t nextBeliefId = nextBelieves[currentBeliefId][chosenActions[currentBeliefId]][observation];
                            if (visitedBelieves.insert(nextBeliefId).second) {
                                beliefStateMap[nextBeliefId] = stateId;
                                ++stateId;
                                believesToBeExpanded.push_back(nextBeliefId);
                            }
                            transitionsInState[beliefStateMap[nextBeliefId]] = iter->second;
                        }
                    }
                    transitions.push_back(transitionsInState);
                    believesToBeExpanded.pop_front();
                }

                for (size_t i = 0; i < transitions.size(); ++i) {
                    for (auto const &transition : transitions[i]) {
                        STORM_LOG_DEBUG(
                                "Transition: " << i << " -- " << transition.second << "--> " << transition.first);
                    }
                }
                storm::models::sparse::StateLabeling labeling(transitions.size());
                labeling.addLabel("init");
                labeling.addLabel("target");
                labeling.addLabelToState("init", 0);
                for (auto targetState : targetStates) {
                    labeling.addLabelToState("target", targetState);
                }
                storm::storage::sparse::ModelComponents<ValueType, RewardModelType> modelComponents(
                        buildTransitionMatrix(transitions), labeling);

                storm::models::sparse::Dtmc<ValueType, RewardModelType> underApproxDtmc(modelComponents);
                auto model = std::make_shared<storm::models::sparse::Dtmc<ValueType, RewardModelType>>(underApproxDtmc);
                model->printModelInformationToStream(std::cout);

                std::string propertyString = min ? "Pmin=? [F \"target\"]" : "Pmax=? [F \"target\"]";
                std::vector<storm::jani::Property> propertyVector = storm::api::parseProperties(propertyString);
                std::shared_ptr<storm::logic::Formula const> property = storm::api::extractFormulasFromProperties(
                        propertyVector).front();

                std::unique_ptr<storm::modelchecker::CheckResult> res(
                        storm::api::verifyWithSparseEngine<ValueType>(model, storm::api::createTask<ValueType>(property,
                                                                                                               true)));
                STORM_LOG_ASSERT(res, "Result does not exist.");
                res->filter(storm::modelchecker::ExplicitQualitativeCheckResult(model->getInitialStates()));
                ValueType resultValue = res->asExplicitQuantitativeCheckResult<ValueType>().getValueMap().begin()->second;


                STORM_PRINT("Over-Approximation Result: " << overApprox << std::endl);
                STORM_PRINT("Under-Approximation Result: " << resultValue << std::endl);
            }

            template<typename ValueType, typename RewardModelType>
            storm::storage::SparseMatrix<ValueType>
            ApproximatePOMDPModelchecker<ValueType, RewardModelType>::buildTransitionMatrix(
                    std::vector<std::map<uint64_t, ValueType>> transitions) {
                uint_fast64_t currentRow = 0;
                uint64_t nrEntries = 0;
                for (auto const &map : transitions) {
                    nrEntries += map.size();
                }
                storm::storage::SparseMatrixBuilder<ValueType> smb(transitions.size(), transitions.size(), nrEntries);
                for (auto const &map : transitions) {
                    for (auto const &transition : map) {
                        smb.addNextValue(currentRow, transition.first, transition.second);
                    }
                    ++currentRow;
                }
                return smb.build();
            }

            template<typename ValueType, typename RewardModelType>
            uint64_t ApproximatePOMDPModelchecker<ValueType, RewardModelType>::extractBestAction(
                    storm::models::sparse::Pomdp<ValueType, RewardModelType> const &pomdp,
                    std::vector<storm::pomdp::Belief<ValueType>> &beliefList,
                    std::vector<bool> &beliefIsTarget,
                    std::set<uint32_t> &targetObservations,
                    std::map<uint64_t, std::vector<std::map<uint32_t, ValueType>>> &observationProbabilities,
                    std::map<uint64_t, std::vector<std::map<uint32_t, uint64_t>>> &nextBelieves,
                    std::map<uint64_t, ValueType> &result,
                    uint64_t gridResolution, uint64_t currentBeliefId, uint64_t nextId, bool min) {
                auto cc = storm::utility::ConstantsComparator<ValueType>();
                storm::pomdp::Belief<ValueType> currentBelief = beliefList[currentBeliefId];
                std::vector<std::map<uint32_t, ValueType>> observationProbabilitiesInAction;
                std::vector<std::map<uint32_t, uint64_t>> nextBelievesInAction;

                uint64_t numChoices = pomdp.getNumberOfChoices(
                        pomdp.getStatesWithObservation(currentBelief.observation).front());
                for (uint64_t action = 0; action < numChoices; ++action) {
                    std::map<uint32_t, ValueType> actionObservationProbabilities = computeObservationProbabilitiesAfterAction(
                            pomdp, currentBelief, action);
                    std::map<uint32_t, uint64_t> actionObservationBelieves;
                    for (auto iter = actionObservationProbabilities.begin();
                        iter != actionObservationProbabilities.end(); ++iter) {
                        uint32_t observation = iter->first;
                        actionObservationBelieves[observation] = getBeliefAfterActionAndObservation(pomdp,
                                                                                                    beliefList,
                                                                                                    beliefIsTarget,
                                                                                                    targetObservations,
                                                                                                    currentBelief,
                                                                                                    action,
                                                                                                    observation,
                                                                                                    nextId);
                        nextId = beliefList.size();
                    }
                    observationProbabilitiesInAction.push_back(actionObservationProbabilities);
                    nextBelievesInAction.push_back(actionObservationBelieves);
                }
                //STORM_LOG_DEBUG("ID " << currentBeliefId << " add " << observationProbabilitiesInAction);
                observationProbabilities.emplace(std::make_pair(currentBeliefId, observationProbabilitiesInAction));
                nextBelieves.emplace(std::make_pair(currentBeliefId, nextBelievesInAction));

                // choose the action which results in the value computed by the over-approximation
                ValueType chosenValue = min ? storm::utility::infinity<ValueType>()
                                            : -storm::utility::infinity<ValueType>();
                uint64_t chosenActionIndex = std::numeric_limits<uint64_t>::infinity();
                ValueType currentValue;

                for (uint64_t action = 0; action < numChoices; ++action) {
                    currentValue = storm::utility::zero<ValueType>(); // simply change this for rewards?
                    for (auto iter = observationProbabilities[currentBelief.id][action].begin();
                         iter != observationProbabilities[currentBelief.id][action].end(); ++iter) {
                        uint32_t observation = iter->first;
                        storm::pomdp::Belief<ValueType> nextBelief = beliefList[nextBelieves[currentBelief.id][action][observation]];

                        // compute subsimplex and lambdas according to the Lovejoy paper to approximate the next belief
                        auto temp = computeSubSimplexAndLambdas(
                                nextBelief.probabilities, gridResolution);
                        std::vector<std::vector<ValueType>> subSimplex = temp.first;
                        std::vector<ValueType> lambdas = temp.second;

                        auto sum = storm::utility::zero<ValueType>();
                        for (size_t j = 0; j < lambdas.size(); ++j) {
                            if (lambdas[j] != storm::utility::zero<ValueType>()) {
                                sum += lambdas[j] * result.at(
                                        getBeliefIdInVector(beliefList, observation, subSimplex[j]));
                            }
                        }
                        currentValue += iter->second * sum;
                    }

                    // Update the selected actions TODO make this nicer
                    if ((min && cc.isLess(storm::utility::zero<ValueType>(), chosenValue - currentValue)) ||
                        (!min &&
                         cc.isLess(storm::utility::zero<ValueType>(), currentValue - chosenValue))) {
                        chosenValue = currentValue;
                        chosenActionIndex = action;
                    } else if ((min && cc.isEqual(storm::utility::zero<ValueType>(), chosenValue - currentValue)) ||
                               (!min &&
                                cc.isEqual(storm::utility::zero<ValueType>(), currentValue - chosenValue))) {
                        chosenValue = currentValue;
                        chosenActionIndex = action;
                    }
                }
                return chosenActionIndex;
            }


            template<typename ValueType, typename RewardModelType>
            uint64_t ApproximatePOMDPModelchecker<ValueType, RewardModelType>::getBeliefIdInVector(
                    std::vector<storm::pomdp::Belief<ValueType>> &grid, uint32_t observation,
                    std::vector<ValueType> probabilities) {
                for (auto const &belief : grid) {
                    if (belief.observation == observation) {
                        if (belief.probabilities == probabilities) {
                            return belief.id;
                        }
                    }
                }
                return -1;
            }

            template<typename ValueType, typename RewardModelType>
            storm::pomdp::Belief<ValueType> ApproximatePOMDPModelchecker<ValueType, RewardModelType>::getInitialBelief(
                    storm::models::sparse::Pomdp<ValueType, RewardModelType> const &pomdp, uint64_t id) {
                STORM_LOG_ASSERT(pomdp.getInitialStates().getNumberOfSetBits() < 2,
                                 "POMDP contains more than one initial state");
                STORM_LOG_ASSERT(pomdp.getInitialStates().getNumberOfSetBits() == 1,
                                 "POMDP does not contain an initial state");
                std::vector<ValueType> distribution(pomdp.getNumberOfStates(), storm::utility::zero<ValueType>());
                uint32_t observation = 0;
                for (uint64_t state = 0; state < pomdp.getNumberOfStates(); ++state) {
                    if (pomdp.getInitialStates()[state] == 1) {
                        distribution[state] = storm::utility::one<ValueType>();
                        observation = pomdp.getObservation(state);
                    }
                }
                return storm::pomdp::Belief<ValueType>{id, observation, distribution};
            }

            template<typename ValueType, typename RewardModelType>
            void ApproximatePOMDPModelchecker<ValueType, RewardModelType>::constructBeliefGrid(
                    storm::models::sparse::Pomdp<ValueType, RewardModelType> const &pomdp,
                    std::set<uint32_t> target_observations, uint64_t gridResolution,
                    std::vector<storm::pomdp::Belief<ValueType>> &beliefList,
                    std::vector<storm::pomdp::Belief<ValueType>> &grid, std::vector<bool> &beliefIsKnown,
                    uint64_t nextId) {
                bool isTarget;
                uint64_t newId = nextId;

                for (uint32_t observation = 0; observation < pomdp.getNrObservations(); ++observation) {
                    std::vector<uint64_t> statesWithObservation = pomdp.getStatesWithObservation(observation);
                    isTarget = target_observations.find(observation) != target_observations.end();

                    // TODO this can probably be condensed
                    if (statesWithObservation.size() == 1) {
                        // If there is only one state with the observation, we can directly add the corresponding belief
                        std::vector<ValueType> distribution(pomdp.getNumberOfStates(),
                                                            storm::utility::zero<ValueType>());
                        distribution[statesWithObservation.front()] = storm::utility::one<ValueType>();
                        storm::pomdp::Belief<ValueType> belief = {newId, observation, distribution};
                        STORM_LOG_TRACE(
                                "Add Belief " << std::to_string(newId) << " [(" << std::to_string(observation) << "),"
                                              << distribution << "]");
                        beliefList.push_back(belief);
                        grid.push_back(belief);
                        beliefIsKnown.push_back(isTarget);
                        ++newId;
                    } else {
                        // Otherwise we have to enumerate all possible distributions with regards to the grid
                        // helper is used to derive the distribution of the belief
                        std::vector<ValueType> helper(statesWithObservation.size(), ValueType(0));
                        helper[0] = storm::utility::convertNumber<ValueType>(gridResolution);
                        bool done = false;
                        uint64_t index = 0;

                        while (!done) {
                            std::vector<ValueType> distribution(pomdp.getNumberOfStates(),
                                                                storm::utility::zero<ValueType>());
                            for (size_t i = 0; i < statesWithObservation.size() - 1; ++i) {
                                distribution[statesWithObservation[i]] = (helper[i] - helper[i + 1]) /
                                                                         storm::utility::convertNumber<ValueType>(
                                                                                 gridResolution);
                            }
                            distribution[statesWithObservation.back()] =
                                    helper[statesWithObservation.size() - 1] /
                                    storm::utility::convertNumber<ValueType>(gridResolution);

                            storm::pomdp::Belief<ValueType> belief = {newId, observation, distribution};
                            STORM_LOG_TRACE(
                                    "Add Belief " << std::to_string(newId) << " [(" << std::to_string(observation)
                                                  << ")," << distribution << "]");
                            beliefList.push_back(belief);
                            grid.push_back(belief);
                            beliefIsKnown.push_back(isTarget);
                            if (helper[statesWithObservation.size() - 1] ==
                                storm::utility::convertNumber<ValueType>(gridResolution)) {
                                // If the last entry of helper is the gridResolution, we have enumerated all necessary distributions
                                done = true;
                            } else {
                                // Update helper by finding the index to increment
                                index = statesWithObservation.size() - 1;
                                while (helper[index] == helper[index - 1]) {
                                    --index;
                                }
                                STORM_LOG_ASSERT(index > 0, "Error in BeliefGrid generation - index wrong");
                                // Increment the value at the index
                                ++helper[index];
                                // Reset all indices greater than the changed one to 0
                                ++index;
                                while (index < statesWithObservation.size()) {
                                    helper[index] = 0;
                                    ++index;
                                }
                            }
                            ++newId;
                        }
                    }
                }
            }

            template<typename ValueType, typename RewardModelType>
            std::pair<std::vector<std::vector<ValueType>>, std::vector<ValueType>>
            ApproximatePOMDPModelchecker<ValueType, RewardModelType>::computeSubSimplexAndLambdas(
                    std::vector<ValueType> probabilities, uint64_t resolution) {
                // This is the Freudenthal Triangulation as described in Lovejoy (a whole lotta math)
                // Variable names are based on the paper

                std::vector<ValueType> x(probabilities.size(), storm::utility::zero<ValueType>());
                std::vector<ValueType> v(probabilities.size(), storm::utility::zero<ValueType>());
                std::vector<ValueType> d(probabilities.size(), storm::utility::zero<ValueType>());

                for (size_t i = 0; i < probabilities.size(); ++i) {
                    for (size_t j = i; j < probabilities.size(); ++j) {
                        x[i] += storm::utility::convertNumber<ValueType>(resolution) * probabilities[j];
                    }
                    v[i] = storm::utility::floor(x[i]);
                    d[i] = x[i] - v[i];
                }

                auto p = storm::utility::vector::getSortedIndices(d);

                std::vector<std::vector<ValueType>> qs;
                for (size_t i = 0; i < probabilities.size(); ++i) {
                    std::vector<ValueType> q(probabilities.size(), storm::utility::zero<ValueType>());
                    if (i == 0) {
                        for (size_t j = 0; j < probabilities.size(); ++j) {
                            q[j] = v[j];
                        }
                        qs.push_back(q);
                    } else {
                        for (size_t j = 0; j < probabilities.size(); ++j) {
                            if (j == p[i - 1]) {
                                q[j] = qs[i - 1][j] + storm::utility::one<ValueType>();
                            } else {
                                q[j] = qs[i - 1][j];
                            }
                        }
                        qs.push_back(q);
                    }
                }

                std::vector<std::vector<ValueType>> subSimplex;
                for (auto q : qs) {
                    std::vector<ValueType> node;
                    for (size_t i = 0; i < probabilities.size(); ++i) {
                        if (i != probabilities.size() - 1) {
                            node.push_back((q[i] - q[i + 1]) / storm::utility::convertNumber<ValueType>(resolution));
                        } else {
                            node.push_back(q[i] / storm::utility::convertNumber<ValueType>(resolution));
                        }
                    }
                    subSimplex.push_back(node);
                }

                std::vector<ValueType> lambdas(probabilities.size(), storm::utility::zero<ValueType>());
                auto sum = storm::utility::zero<ValueType>();
                for (size_t i = 1; i < probabilities.size(); ++i) {
                    lambdas[i] = d[p[i - 1]] - d[p[i]];
                    sum += d[p[i - 1]] - d[p[i]];
                }
                lambdas[0] = storm::utility::one<ValueType>() - sum;

                //TODO add assertion that we are close enough
                return std::make_pair(subSimplex, lambdas);
            }


            template<typename ValueType, typename RewardModelType>
            std::map<uint32_t, ValueType>
            ApproximatePOMDPModelchecker<ValueType, RewardModelType>::computeObservationProbabilitiesAfterAction(
                    storm::models::sparse::Pomdp<ValueType, RewardModelType> const &pomdp,
                    storm::pomdp::Belief<ValueType> belief,
                    uint64_t actionIndex) {
                std::map<uint32_t, ValueType> res;
                // the id is not important here as we immediately discard the belief (very hacky, I don't like it either)
                std::vector<ValueType> postProbabilities = getBeliefAfterAction(pomdp, belief, actionIndex,
                                                                                0).probabilities;
                for (uint64_t state = 0; state < pomdp.getNumberOfStates(); ++state) {
                    uint32_t observation = pomdp.getObservation(state);
                    if (postProbabilities[state] != storm::utility::zero<ValueType>()) {
                        if (res.count(observation) == 0) {
                            res[observation] = postProbabilities[state];
                        } else {
                            res[observation] += postProbabilities[state];
                        }
                    }
                }
                return res;
            }

            template<typename ValueType, typename RewardModelType>
            storm::pomdp::Belief<ValueType>
            ApproximatePOMDPModelchecker<ValueType, RewardModelType>::getBeliefAfterAction(
                    storm::models::sparse::Pomdp<ValueType, RewardModelType> const &pomdp,
                    storm::pomdp::Belief<ValueType> belief,
                    uint64_t actionIndex, uint64_t id) {
                std::vector<ValueType> distributionAfter(pomdp.getNumberOfStates(), storm::utility::zero<ValueType>());
                uint32_t observation = 0;
                for (uint64_t state = 0; state < pomdp.getNumberOfStates(); ++state) {
                    if (belief.probabilities[state] != storm::utility::zero<ValueType>()) {
                        auto row = pomdp.getTransitionMatrix().getRow(
                                pomdp.getChoiceIndex(storm::storage::StateActionPair(state, actionIndex)));
                        for (auto const &entry : row) {
                            observation = pomdp.getObservation(entry.getColumn());
                            distributionAfter[entry.getColumn()] += belief.probabilities[state] * entry.getValue();
                        }
                    }
                }
                return storm::pomdp::Belief<ValueType>{id, observation, distributionAfter};
            }

            template<typename ValueType, typename RewardModelType>
            uint64_t ApproximatePOMDPModelchecker<ValueType, RewardModelType>::getBeliefAfterActionAndObservation(
                    storm::models::sparse::Pomdp<ValueType, RewardModelType> const &pomdp,
                    std::vector<storm::pomdp::Belief<ValueType>> &beliefList, std::vector<bool> &beliefIsTarget,
                    std::set<uint32_t> &targetObservations, storm::pomdp::Belief<ValueType> belief,
                    uint64_t actionIndex, uint32_t observation, uint64_t id) {
                std::vector<ValueType> distributionAfter(pomdp.getNumberOfStates(), storm::utility::zero<ValueType>());
                for (uint64_t state = 0; state < pomdp.getNumberOfStates(); ++state) {
                    if (belief.probabilities[state] != storm::utility::zero<ValueType>()) {
                        auto row = pomdp.getTransitionMatrix().getRow(
                                pomdp.getChoiceIndex(storm::storage::StateActionPair(state, actionIndex)));
                        for (auto const &entry : row) {
                            if (pomdp.getObservation(entry.getColumn()) == observation) {
                                distributionAfter[entry.getColumn()] += belief.probabilities[state] * entry.getValue();
                            }
                        }
                    }
                }
                // We have to normalize the distribution
                auto sum = storm::utility::zero<ValueType>();
                for (ValueType const &entry : distributionAfter) {
                    sum += entry;
                }
                for (size_t i = 0; i < pomdp.getNumberOfStates(); ++i) {
                    distributionAfter[i] /= sum;
                }
                if (getBeliefIdInVector(beliefList, observation, distributionAfter) != uint64_t(-1)) {
                    return getBeliefIdInVector(beliefList, observation, distributionAfter);
                } else {
                    beliefList.push_back(storm::pomdp::Belief<ValueType>{id, observation, distributionAfter});
                    beliefIsTarget.push_back(targetObservations.find(observation) != targetObservations.end());
                    return id;
                }
            }


            template
            class ApproximatePOMDPModelchecker<double>;

#ifdef STORM_HAVE_CARL

            //template class ApproximatePOMDPModelchecker<storm::RationalFunction>;
            template
            class ApproximatePOMDPModelchecker<storm::RationalNumber>;

#endif
        }
    }
}