/* 
 * File:   ParameterRegion.h
 * Author: tim
 *
 * Created on August 10, 2015, 1:51 PM
 */

#ifndef STORM_MODELCHECKER_REGION_PARAMETERREGION_H
#define	STORM_MODELCHECKER_REGION_PARAMETERREGION_H

#include <map>

#include "src/modelchecker/region/RegionCheckResult.h"
#include "src/utility/region.h"


namespace storm {
    namespace modelchecker{
        namespace region {
            template<typename ParametricType>
            class ParameterRegion{
            public:
                typedef typename storm::utility::region::VariableType<ParametricType> VariableType;
                typedef typename storm::utility::region::CoefficientType<ParametricType> CoefficientType;
                typedef typename std::map<VariableType, CoefficientType> VariableSubstitutionType;

                ParameterRegion(VariableSubstitutionType const& lowerBounds, VariableSubstitutionType const& upperBounds);
                ParameterRegion(VariableSubstitutionType&& lowerBounds, VariableSubstitutionType&& upperBounds);
                virtual ~ParameterRegion();

                std::set<VariableType> getVariables() const;
                CoefficientType const& getLowerBound(VariableType const& variable) const;
                CoefficientType const& getUpperBound(VariableType const& variable) const;
                const VariableSubstitutionType getUpperBounds() const;
                const VariableSubstitutionType getLowerBounds() const;

                /*!
                 * Returns a vector of all possible combinations of lower and upper bounds of the given variables.
                 * The first entry of the returned vector will map every variable to its lower bound
                 * The second entry will map every variable to its lower bound, except the first one (i.e. *getVariables.begin())
                 * ...
                 * The last entry will map every variable to its upper bound
                 * 
                 * If the given set of variables is empty, the returned vector will contain an empty map
                 */
                std::vector<VariableSubstitutionType> getVerticesOfRegion(std::set<VariableType> const& consideredVariables) const;

                /*!
                 * Returns some point that lies within this region
                 */
                VariableSubstitutionType getSomePoint() const;

                RegionCheckResult getCheckResult() const;
                void setCheckResult(RegionCheckResult checkResult);

                /*!
                 * Retrieves a point in the region for which is considered property is not satisfied.
                 * If such a point is not known, the returned map is empty.
                 */
                VariableSubstitutionType getViolatedPoint() const;

                /*!
                 * Sets a point in the region for which the considered property is not satisfied. 
                 */
                void setViolatedPoint(VariableSubstitutionType const& violatedPoint);

                /*!
                 * Retrieves a point in the region for which is considered property is satisfied.
                 * If such a point is not known, the returned map is empty.
                 */
                VariableSubstitutionType getSatPoint() const;

                /*!
                 * Sets a point in the region for which the considered property is satisfied. 
                 */
                void setSatPoint(VariableSubstitutionType const& satPoint);

                //returns the region as string in the format 0.3<=p<=0.4,0.2<=q<=0.5;
                std::string toString() const;

                /*
                 * Can be used to parse a single parameter with its bounds from a string of the form "0.3<=p<=0.5".
                 * The numbers are parsed as doubles and then converted to SparseDtmcRegionModelChecker::CoefficientType.
                 * The results will be inserted in the given maps
                 * 
                 */
                static void parseParameterBounds( 
                        VariableSubstitutionType& lowerBounds,
                        VariableSubstitutionType& upperBounds,
                        std::string const& parameterBoundsString
                );

                /*
                 * Can be used to parse a single region from a string of the form "0.3<=p<=0.5,0.4<=q<=0.7".
                 * The numbers are parsed as doubles and then converted to SparseDtmcRegionModelChecker::CoefficientType.
                 * 
                 */
                static ParameterRegion parseRegion(
                        std::string const& regionString
                );

                /*
                 * Can be used to parse a vector of region from a string of the form "0.3<=p<=0.5,0.4<=q<=0.7;0.1<=p<=0.3,0.2<=q<=0.4".
                 * The numbers are parsed as doubles and then converted to SparseDtmcRegionModelChecker::CoefficientType.
                 * 
                 */
                static std::vector<ParameterRegion> parseMultipleRegions(
                        std::string const& regionsString
                );


                /*
                 * Retrieves the regions that are specified in the settings.
                 */
                static std::vector<ParameterRegion> getRegionsFromSettings();

                private:

                void init();

                VariableSubstitutionType const lowerBounds;
                VariableSubstitutionType const upperBounds;
                std::set<VariableType> variables;
                RegionCheckResult checkResult;
                VariableSubstitutionType satPoint;
                VariableSubstitutionType violatedPoint;
            };
        } //namespace region
    }
}

#endif	/* STORM_MODELCHECKER_REGION_PARAMETERREGION_H */

