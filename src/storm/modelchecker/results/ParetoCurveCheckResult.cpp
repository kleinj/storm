#include "storm/modelchecker/results/ParetoCurveCheckResult.h"

#include "storm/adapters/RationalFunctionAdapter.h"
#include "storm/utility/vector.h"
#include "storm/exceptions/InvalidOperationException.h"

namespace storm {
    namespace modelchecker {

        template<typename ValueType>
        ParetoCurveCheckResult<ValueType>::ParetoCurveCheckResult() {
            // Intentionally left empty.
        }

        template<typename ValueType>
        ParetoCurveCheckResult<ValueType>::ParetoCurveCheckResult(std::vector<point_type> const& points, polytope_type const& underApproximation, polytope_type const& overApproximation) : points(points), underApproximation(underApproximation), overApproximation(overApproximation) {
            // Intentionally left empty.
        }

        template<typename ValueType>
        ParetoCurveCheckResult<ValueType>::ParetoCurveCheckResult(std::vector<point_type>&& points, polytope_type&& underApproximation, polytope_type&& overApproximation) : points(points), underApproximation(underApproximation), overApproximation(overApproximation) {
            // Intentionally left empty.
        }

        template<typename ValueType>
        bool ParetoCurveCheckResult<ValueType>::isParetoCurveCheckResult() const {
            return true;
        }
        
        template<typename ValueType>
        std::vector<typename ParetoCurveCheckResult<ValueType>::point_type> const& ParetoCurveCheckResult<ValueType>::getPoints() const {
            return points;
        }
        
        template<typename ValueType>
        typename ParetoCurveCheckResult<ValueType>::polytope_type const& ParetoCurveCheckResult<ValueType>::getUnderApproximation() const {
            return underApproximation;
        }
        
        template<typename ValueType>
        typename ParetoCurveCheckResult<ValueType>::polytope_type const& ParetoCurveCheckResult<ValueType>::getOverApproximation() const {
            return overApproximation;
        }
        
        template<typename ValueType>
        std::ostream& ParetoCurveCheckResult<ValueType>::writeToStream(std::ostream& out) const {
            out << std::endl;
            out << "Underapproximation of achievable values: " << underApproximation->toString() << std::endl;
            out << "Overapproximation of achievable values: " << overApproximation->toString() << std::endl;
            out << points.size() << " pareto optimal points found (Note that these points are safe, i.e., contained in the underapproximation, but there is no guarantee for optimality):" << std::endl;
            for(auto const& p : points) {
                out << "   (";
                for(auto it = p.begin(); it != p.end(); ++it){
                    if(it != p.begin()){
                        out << ", ";
                    }
                    out << std::setw(10) << *it;
                }
                out << " )" << std::endl;
            }
            return out;
        }

        template<typename ValueType>
        std::ostream& ParetoCurveCheckResult<ValueType>::writeAsSparseVectorToStream(std::ostream& out) const {
            STORM_LOG_THROW(false, storm::exceptions::InvalidOperationException, "Can not yet print explicit qualitative results as vector.");
        }
      
        template class ParetoCurveCheckResult<double>;
        
#ifdef STORM_HAVE_CARL
        template class ParetoCurveCheckResult<storm::RationalNumber>;
#endif
    }
}
