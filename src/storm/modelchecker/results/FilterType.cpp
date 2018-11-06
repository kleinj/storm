#include "FilterType.h"


namespace storm {
    namespace modelchecker {
        
        std::string toString(FilterType ft) {
            switch(ft) {
                case FilterType::ARGMAX:
                    return "the argmax";
                case FilterType::ARGMIN:
                    return "the argmin";
                case FilterType::AVG:
                    return "the average";
                case FilterType::COUNT:
                    return "the number of";
                case FilterType::EXISTS:
                    return "whether there exists a state in";
                case FilterType::FORALL:
                    return "whether for all states in";
                case FilterType::MAX:
                    return "the maximum";
                case FilterType::MIN:
                    return "the minumum";
                case FilterType::SUM:
                    return "the sum";
                case FilterType::VALUES:
                    return "the values";
                case FilterType::PRINT:
                    return "printing";
            }
        }
    }
}
