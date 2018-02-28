#pragma once

#include <string>
#include <vector>
#include <memory>

namespace storm {
namespace storage {

/**
 * Data structure containing info about SCCs. Assumes that
 * the state indices contained in each SCC are consecutive.
 */
class SCCInfo {
public:
    typedef std::shared_ptr<SCCInfo> ptr;

    SCCInfo() = default;
    SCCInfo(const SCCInfo& other) = default;
    SCCInfo(SCCInfo&& other) = default;

    /** construction: notify about start of new SCC */
    void startSCC(bool isTrivial) {
        curIndexInSCC = 0;
        if (isTrivial) {
            trivialSCCCount++;
        }
    }

    /** construction: notify about end of current SCC */
    void endSCC() {
        curSCC++;
    }

    /** construction: add a state to the current SCC */
    void addState(std::size_t state) {
        assert(state == stateToSCC.size());
        stateToSCC.push_back(curSCC);
        if (curIndexInSCC == 0) {
            // first state in this SCC
            sccToFirstStateInSCC.push_back(state);
        }
        stateToIndexInSCC.push_back(curIndexInSCC++);
    }

    /** Return the SCC index for the given state */
    std::size_t getSCC(std::size_t state) const {
        return stateToSCC.at(state);
    }

    /** Return the index inside its SCC for the given state */
    std::size_t getIndexInSCC(std::size_t state) const {
        return stateToIndexInSCC.at(state);
    }

    /** Return the number of distinct SCCs */
    std::size_t getNumberOfSCCs() const {
        return sccToFirstStateInSCC.size();
    }

    /** Return the number of states in the given SCC */
    std::size_t getNumberOfStatesInSCC(std::size_t scc_index) const {
        return getLastStateInSCC(scc_index) - getFirstStateInSCC(scc_index) + 1;
    }

    /** Return the first state in the SCC with the given index */
    std::size_t getFirstStateInSCC(std::size_t scc_index) const {
        return sccToFirstStateInSCC.at(scc_index);
    }

    /** Return the last state in the SCC with the given index */
    std::size_t getLastStateInSCC(std::size_t scc_index) const {
        if (scc_index + 1 == getNumberOfSCCs()) {
            // last SCC, return the last state index
            return stateToSCC.size()-1;
        } else {
            // get first state in the *next* SCC, subtract one
            // to get the last one in this SCC
            return getFirstStateInSCC(scc_index + 1) - 1;
        }
    }

    std::string getStatistics() const {
        std::string result("Number of SCCs: " + std::to_string(getNumberOfSCCs()));

        std::size_t singletons = 0;
        for (std::size_t i = 0; i < getNumberOfSCCs(); i++) {
            if (getNumberOfStatesInSCC(i) == 1) {
                singletons++;
            }
        }
        result += ", of which " + std::to_string(singletons) + " are singleton SCCs, of which " + std::to_string(trivialSCCCount) + " are trivial.";
        return result;
    }

private:
    std::vector<std::size_t> stateToSCC;
    std::vector<std::size_t> stateToIndexInSCC;
    std::vector<std::size_t> states;
    std::vector<std::size_t> sccToFirstStateInSCC;

    std::size_t curSCC = 0;
    std::size_t curIndexInSCC = 0;
    std::size_t trivialSCCCount = 0;
};


}
}
