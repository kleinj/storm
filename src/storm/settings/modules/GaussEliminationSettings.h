#pragma once

#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/ModuleSettings.h"

namespace storm {
    namespace settings {
        namespace modules {
            
            /*!
             * This class represents the settings for the Gauss elimination-based procedures.
             */
            class GaussEliminationSettings : public ModuleSettings {
            public:

                /*!
                 * An enum that contains all available sub-methods.
                 */
                enum class Method { Standard, Factorized, Plain, FactorizedNoGcd, PlainNoGcd, Fractionfree};

                /*!
                 * An enum that contains all available sub-methods.
                 */
                enum class Matrix { Sparse, Dense };

                /*!
                 * Creates a new set of parametric model checking settings.
                 */
                GaussEliminationSettings();

                /*!
                 * Retrieves whether we should print debug information.
                 *
                 * @return True iff the option was set.
                 */
                bool isDebug() const;

                bool isExportHTMLSet() const;
                std::string getExportHTMLFilename() const;
                bool isExportHTMLComplexitySet() const;
                bool isExportHTMLStatsSet() const;

                bool isTopologicalOrderingSet() const;
                bool isIndividualSCCsSet() const;
                bool isSCCDenominatorsSet() const;
                bool isReportReductionsSet() const;

                Method getMethod() const;
                Matrix getMatrix() const;

                const static std::string moduleName;
                
            private:
                const static std::string debugOptionName;
                const static std::string noTopologicalOrderingOptionName;
                const static std::string noIndividualSCCsOptionName;
                const static std::string SCCDenominatorsOptionName;
                const static std::string reportReductionsOptionName;
                const static std::string methodOptionName;
                const static std::string matrixOptionName;
                const static std::string exportHTMLOptionName;
                const static std::string exportHTMLComplexityOptionName;
                const static std::string exportHTMLStatsOptionName;
            };
            
        } // namespace modules
    } // namespace settings
} // namespace storm
