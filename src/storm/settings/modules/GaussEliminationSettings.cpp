#include "storm/settings/modules/GaussEliminationSettings.h"

#include "storm/settings/Option.h"
#include "storm/settings/OptionBuilder.h"
#include "storm/settings/ArgumentBuilder.h"
#include "storm/settings/Argument.h"

#include "storm/utility/macros.h"
#include "storm/exceptions/IllegalArgumentValueException.h"

namespace storm {
    namespace settings {
        namespace modules {

            const std::string GaussEliminationSettings::moduleName = "gausselimination";
            const std::string GaussEliminationSettings::methodOptionName = "method";
            const std::string GaussEliminationSettings::matrixOptionName = "matrix";
            const std::string GaussEliminationSettings::debugOptionName = "debug";
            const std::string GaussEliminationSettings::noTopologicalOrderingOptionName = "notopo";
            const std::string GaussEliminationSettings::noIndividualSCCsOptionName = "noindividualsccs";
            const std::string GaussEliminationSettings::SCCDenominatorsOptionName = "sccdenominators";
            const std::string GaussEliminationSettings::reportReductionsOptionName = "reportreductions";
            const std::string GaussEliminationSettings::exportHTMLOptionName = "exporthtml";
            const std::string GaussEliminationSettings::exportHTMLComplexityOptionName = "exporthtmlcomplexity";
            const std::string GaussEliminationSettings::exportHTMLStatsOptionName = "exporthtmlstats";

            GaussEliminationSettings::GaussEliminationSettings() : ModuleSettings(moduleName) {
                std::vector<std::string> methods = {"standard", "factorized", "plain", "factorizednogcd", "plainnogcd", "fractionfree"};
                this->addOption(storm::settings::OptionBuilder(moduleName, methodOptionName, true, "The elimination technique to use.")
                                .addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The name of the elimination technique to use.").addValidatorString(ArgumentValidatorFactory::createMultipleChoiceValidator(methods)).setDefaultValueString("standard").build()).build());
                std::vector<std::string> matrizes = {"sparse", "dense"};
                this->addOption(storm::settings::OptionBuilder(moduleName, matrixOptionName, true, "The matrix type to use during Gaussian elimination.")
                                .addArgument(storm::settings::ArgumentBuilder::createStringArgument("type", "The type of the matrix to use during Gaussian elimination.").addValidatorString(ArgumentValidatorFactory::createMultipleChoiceValidator(matrizes)).setDefaultValueString("sparse").build()).build());

                this->addOption(storm::settings::OptionBuilder(moduleName, noTopologicalOrderingOptionName, true, "Don't perform topological ordering").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, noIndividualSCCsOptionName, true, "Don't handle each SCC separately").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, SCCDenominatorsOptionName, true, "Use per-SCC denominator").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, reportReductionsOptionName, true, "Report reductions due to gcd-simplifications").build());

                this->addOption(storm::settings::OptionBuilder(moduleName, debugOptionName, true, "Debug mode").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, exportHTMLOptionName, "", "Export the various matrizes to HTML file")
                                .addArgument(storm::settings::ArgumentBuilder::createStringArgument("filename", "the name of the HTML file.").build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, exportHTMLComplexityOptionName, true, "Show complexity in HTML export").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, exportHTMLStatsOptionName, true, "Show statistics in HTML export").build());
            }

            bool GaussEliminationSettings::isDebug() const {
                return this->getOption(debugOptionName).getHasOptionBeenSet();
            }

            bool GaussEliminationSettings::isTopologicalOrderingSet() const {
                return !this->getOption(noTopologicalOrderingOptionName).getHasOptionBeenSet();
            }

            bool GaussEliminationSettings::isIndividualSCCsSet() const {
                return !this->getOption(noIndividualSCCsOptionName).getHasOptionBeenSet();
            }

            bool GaussEliminationSettings::isSCCDenominatorsSet() const {
                return this->getOption(SCCDenominatorsOptionName).getHasOptionBeenSet();
            }

            bool GaussEliminationSettings::isReportReductionsSet() const {
                return this->getOption(reportReductionsOptionName).getHasOptionBeenSet();
            }

            bool GaussEliminationSettings::isExportHTMLSet() const {
                return this->getOption(exportHTMLOptionName).getHasOptionBeenSet();
            }

            std::string GaussEliminationSettings::getExportHTMLFilename() const {
                return this->getOption(exportHTMLOptionName).getArgumentByName("filename").getValueAsString();
            }

            bool GaussEliminationSettings::isExportHTMLComplexitySet() const {
                return this->getOption(exportHTMLComplexityOptionName).getHasOptionBeenSet();
            }

            bool GaussEliminationSettings::isExportHTMLStatsSet() const {
                return this->getOption(exportHTMLStatsOptionName).getHasOptionBeenSet();
            }

            GaussEliminationSettings::Method GaussEliminationSettings::getMethod() const {
                std::string methodAsString = this->getOption(methodOptionName).getArgumentByName("name").getValueAsString();
                if (methodAsString == "standard") {
                    return Method::Standard;
                } else if (methodAsString == "factorized") {
                    return Method::Factorized;
                } else if (methodAsString == "plain") {
                    return Method::Plain;
                } else if (methodAsString == "factorizednogcd") {
                    return Method::FactorizedNoGcd;
                } else if (methodAsString == "plainnogcd") {
                    return Method::PlainNoGcd;
                } else if (methodAsString == "fractionfree") {
                    return Method::Fractionfree;
                } else {
                    STORM_LOG_THROW(false, storm::exceptions::IllegalArgumentValueException, "Illegal method selected.");
                }
            }

            GaussEliminationSettings::Matrix GaussEliminationSettings::getMatrix() const {
                std::string matrixAsString = this->getOption(matrixOptionName).getArgumentByName("type").getValueAsString();
                if (matrixAsString == "sparse") {
                    return Matrix::Sparse;
                } else if (matrixAsString == "dense") {
                    return Matrix::Dense;
                } else {
                    STORM_LOG_THROW(false, storm::exceptions::IllegalArgumentValueException, "Illegal matrix type selected.");
                }
            }


        } // namespace modules
    } // namespace settings
} // namespace storm
