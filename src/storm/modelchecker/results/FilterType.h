#pragma once
#include <string>

namespace storm {
    namespace modelchecker {
        
        enum class StateFilter { ARGMIN, ARGMAX };
        
        enum class FilterType { MIN, MAX, SUM, AVG, COUNT, FORALL, EXISTS, ARGMIN, ARGMAX, VALUES, PRINT };
        
        std::string toString(FilterType);
        bool isStateFilter(FilterType);
    }
}
