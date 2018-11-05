#ifndef STORM_STORAGE_DD_INTERNALCUDDDDMANAGER_H_
#define STORM_STORAGE_DD_INTERNALCUDDDDMANAGER_H_

#include <boost/optional.hpp>

#include "storm/storage/dd/DdType.h"
#include "storm/storage/dd/InternalDdManager.h"

#include "storm/storage/dd/cudd/InternalCuddBdd.h"
#include "storm/storage/dd/cudd/InternalCuddAdd.h"

#include "cuddObj.hh"

namespace storm {
    namespace dd {
        template<DdType LibraryType, typename ValueType>
        class InternalAdd;

        template<DdType LibraryType>
        class InternalBdd;

        template<>
        class InternalDdManager<DdType::CUDD> {
        public:
            friend class InternalBdd<DdType::CUDD>;

            template<DdType LibraryType, typename ValueType>
            friend class InternalAdd;

            // helper class to inhibit dynamic reordering as long as the object is in scope
            // increments reorderingInhibitionCounter on construction for the given manager,
            // decrements on destruction
            class DynamicReorderingInhibitor {
            public:
                DynamicReorderingInhibitor(InternalDdManager<DdType::CUDD>& manager);
                ~DynamicReorderingInhibitor();

                // inhibit implicit copy, assignment and move
                DynamicReorderingInhibitor(DynamicReorderingInhibitor const&) = delete;
                DynamicReorderingInhibitor(DynamicReorderingInhibitor&&) = delete;
                DynamicReorderingInhibitor& operator=(DynamicReorderingInhibitor const&) = delete;

            private:
                InternalDdManager<DdType::CUDD>& manager;
            };

	    
            /*!
             * Creates a new internal manager for CUDD DDs.
             */
            InternalDdManager();
            
            /*!
             * Destroys the CUDD manager.
             */
            ~InternalDdManager();
            
            /*!
             * Retrieves a BDD representing the constant one function.
             *
             * @return A BDD representing the constant one function.
             */
            InternalBdd<DdType::CUDD> getBddOne() const;
            
            /*!
             * Retrieves an ADD representing the constant one function.
             *
             * @return An ADD representing the constant one function.
             */
            template<typename ValueType>
            InternalAdd<DdType::CUDD, ValueType> getAddOne() const;
            
            /*!
             * Retrieves a BDD representing the constant zero function.
             *
             * @return A BDD representing the constant zero function.
             */
            InternalBdd<DdType::CUDD> getBddZero() const;
            
            /*!
             * Retrieves a BDD that maps to true iff the encoding is less or equal than the given bound.
             *
             * @return A BDD with encodings corresponding to values less or equal than the bound.
             */
            InternalBdd<DdType::CUDD> getBddEncodingLessOrEqualThan(uint64_t bound, InternalBdd<DdType::CUDD> const& cube, uint64_t numberOfDdVariables) const;
            
            /*!
             * Retrieves an ADD representing the constant zero function.
             *
             * @return An ADD representing the constant zero function.
             */
            template<typename ValueType>
            InternalAdd<DdType::CUDD, ValueType> getAddZero() const;

            /*!
             * Retrieves an ADD representing an undefined value.
             *
             * @return An ADD representing an undefined value.
             */
            template<typename ValueType>
            InternalAdd<DdType::CUDD, ValueType> getAddUndefined() const;

            /*!
             * Retrieves an ADD representing the constant function with the given value.
             *
             * @return An ADD representing the constant function with the given value.
             */
            template<typename ValueType>
            InternalAdd<DdType::CUDD, ValueType> getConstant(ValueType const& value) const;
            
            /*!
             * Creates new layered DD variables and returns the cubes as a result.
             *
             * @param position An optional position at which to insert the new variable. This may only be given, if the
             * manager supports ordered insertion.
             * @return The cubes belonging to the DD variables.
             */
            std::vector<InternalBdd<DdType::CUDD>> createDdVariables(uint64_t numberOfLayers, boost::optional<uint_fast64_t> const& position = boost::none);
            
            /*!
             * Checks whether this manager supports the ordered insertion of variables, i.e. inserting variables at
             * positions between already existing variables.
             *
             * @return True iff the manager supports ordered insertion.
             */
            bool supportsOrderedInsertion() const;
            
            /*!
             * Sets whether or not dynamic reordering is allowed for the DDs managed by this manager.
             *
             * @param value If set to true, dynamic reordering is allowed and forbidden otherwise.
             */
            void allowDynamicReordering(bool value);
            
            /*!
             * Retrieves whether dynamic reordering is currently allowed.
             *
             * @return True iff dynamic reordering is currently allowed.
             */
            bool isDynamicReorderingAllowed() const;
            
            /*!
             * Triggers a reordering of the DDs managed by this manager.
             */
            void triggerReordering();

            /*!
             * Requests an object that inhibits dynamic reordering as long as it is
             * in scope.
             */
            std::unique_ptr<DynamicReorderingInhibitor> getDynamicReorderingInhibitor() const;
            
            /*!
             * Performs a debug check if available.
             */
            void debugCheck() const;
            
            /*!
             * Retrieves the number of DD variables managed by this manager.
             *
             * @return The number of managed variables.
             */
            uint_fast64_t getNumberOfDdVariables() const;

            /*!
             * Retrieves the underlying CUDD manager.
             *
             * @return The underlying CUDD manager.
             */
            cudd::Cudd& getCuddManager();

            /*!
             * Retrieves the underlying CUDD manager.
             *
             * @return The underlying CUDD manager.
             */
            cudd::Cudd const& getCuddManager() const;

        private:
            // Communicate to CUDD whether currently dynamic reordering is allowed.
            void setDynamicReorderingState();

            // The manager responsible for the DDs created/modified with this DdManager.
            cudd::Cudd cuddManager;
            
            // The technique that is used for dynamic reordering.
            Cudd_ReorderingType reorderingTechnique;
            
            // Keeps track of the number of registered DD variables.
            uint_fast64_t numberOfDdVariables;

            // Flag: Is reordering allowed in principle?
            bool allowReorder;

            // Counter to inhibit reordering: 0 = may occur, >0 = one or more DynamicReorderingInhibitors are active
            uint32_t reorderingInhibitionCounter;
        };        
    }
}

#endif /* STORM_STORAGE_DD_INTERNALCUDDDDMANAGER_H_ */
