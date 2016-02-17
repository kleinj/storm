#include "src/storage/TotalScheduler.h"
#include "src/exceptions/InvalidArgumentException.h"

namespace storm {
    namespace storage {
        TotalScheduler::TotalScheduler(uint_fast64_t numberOfStates) : choices(numberOfStates) {
            // Intentionally left empty.
        }
        
        TotalScheduler::TotalScheduler(std::vector<uint_fast64_t> const& choices) : choices(choices) {
            // Intentionally left empty.
        }
        
        TotalScheduler::TotalScheduler(std::vector<uint_fast64_t>&& choices) : choices(std::move(choices)) {
            // Intentionally left empty.
        }
        
        void TotalScheduler::setChoice(uint_fast64_t state, uint_fast64_t choice) {
            if (state > choices.size()) {
                throw storm::exceptions::InvalidArgumentException() << "Invalid call to TotalScheduler::setChoice: scheduler cannot not define a choice for state " << state << ".";
            }
            choices[state] = choice;
        }
        
        bool TotalScheduler::isChoiceDefined(uint_fast64_t state) const {
            return state < choices.size();
        }
        
        uint_fast64_t TotalScheduler::getChoice(uint_fast64_t state) const {
            if (state >= choices.size()) {
                throw storm::exceptions::InvalidArgumentException() << "Invalid call to TotalScheduler::getChoice: scheduler does not define a choice for state " << state << ".";
            }

            return choices[state];
        }
        
        std::ostream& operator<<(std::ostream& out, TotalScheduler const& scheduler) {
            out << "total scheduler (defined on " << scheduler.choices.size() << " states) [ ";
            for (uint_fast64_t state = 0; state < scheduler.choices.size() - 1; ++state) {
                out << state << " -> " << scheduler.choices[state] << ", ";
            }
            if (scheduler.choices.size() > 0) {
                out << (scheduler.choices.size() - 1) << " -> " << scheduler.choices[scheduler.choices.size() - 1] << " ]";
            }
            return out;
        }

    }
}