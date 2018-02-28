#include "storm/settings/modules/EliminationSettings.h"

#include "storm/settings/Option.h"
#include "storm/settings/OptionBuilder.h"
#include "storm/settings/ArgumentBuilder.h"
#include "storm/settings/Argument.h"

#include "storm/utility/macros.h"
#include "storm/exceptions/IllegalArgumentValueException.h"

namespace storm {
    namespace settings {
        namespace modules {
            
            const std::string EliminationSettings::moduleName = "elimination";
            const std::string EliminationSettings::eliminationMethodOptionName = "method";
            const std::string EliminationSettings::eliminationOrderOptionName = "order";
            const std::string EliminationSettings::entryStatesLastOptionName = "entrylast";
            const std::string EliminationSettings::maximalSccSizeOptionName = "sccsize";
            const std::string EliminationSettings::useDedicatedModelCheckerOptionName = "use-dedicated-mc";
            const std::string EliminationSettings::topologicalOrderingOptionName = "topo";
            const std::string EliminationSettings::exportHTMLOptionName = "exporthtml";
            const std::string EliminationSettings::exportHTMLComplexityOptionName = "exporthtmlcomplexity";
            const std::string EliminationSettings::exportHTMLStatsOptionName = "exporthtmlstats";

            
            EliminationSettings::EliminationSettings() : ModuleSettings(moduleName) {
                std::vector<std::string> orders = {"fw", "fwrev", "bw", "bwrev", "rand", "spen", "dpen", "regex"};
                this->addOption(storm::settings::OptionBuilder(moduleName, eliminationOrderOptionName, true, "The order that is to be used for the elimination techniques.").addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The name of the order in which states are chosen for elimination.").addValidatorString(ArgumentValidatorFactory::createMultipleChoiceValidator(orders)).setDefaultValueString("fwrev").build()).build());
                
                std::vector<std::string> methods = {"state", "hybrid"};
                this->addOption(storm::settings::OptionBuilder(moduleName, eliminationMethodOptionName, true, "The elimination technique to use.").addArgument(storm::settings::ArgumentBuilder::createStringArgument("name", "The name of the elimination technique to use.").addValidatorString(ArgumentValidatorFactory::createMultipleChoiceValidator(methods)).setDefaultValueString("state").build()).build());
                
                this->addOption(storm::settings::OptionBuilder(moduleName, entryStatesLastOptionName, true, "Sets whether the entry states are eliminated last.").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, maximalSccSizeOptionName, true, "Sets the maximal size of the SCCs for which state elimination is applied.")
                                .addArgument(storm::settings::ArgumentBuilder::createUnsignedIntegerArgument("maxsize", "The maximal size of an SCC on which state elimination is applied.").setDefaultValueUnsignedInteger(20).setIsOptional(true).build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, useDedicatedModelCheckerOptionName, true, "Sets whether to use the dedicated model elimination checker (only DTMCs).").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, topologicalOrderingOptionName, false, "Perform topological ordering").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, exportHTMLOptionName, "", "Export the various matrizes to HTML file")
                                .addArgument(storm::settings::ArgumentBuilder::createStringArgument("filename", "the name of the HTML file.").build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, exportHTMLComplexityOptionName, true, "Show complexity in HTML export").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, exportHTMLStatsOptionName, true, "Show statistics in HTML export").build());

            }
            
            std::string EliminationSettings::getEliminationMethodAsString() const {
                return this->getOption(eliminationMethodOptionName).getArgumentByName("name").getValueAsString();
            }

            EliminationSettings::EliminationMethod EliminationSettings::getEliminationMethod() const {
                std::string method = getEliminationMethodAsString();

                if (method == "state") {
                    return EliminationMethod::State;
                } else if (method == "hybrid") {
                    return EliminationMethod::Hybrid;
                } else {
                    STORM_LOG_THROW(false, storm::exceptions::IllegalArgumentValueException, "Illegal elimination method selected.");
                }
            }

            std::string EliminationSettings::getEliminationOrderAsString() const {
                return this->getOption(eliminationOrderOptionName).getArgumentByName("name").getValueAsString();
            }

            EliminationSettings::EliminationOrder EliminationSettings::getEliminationOrder() const {
                std::string order = getEliminationOrderAsString();

                if (order == "fw") {
                    return EliminationOrder::Forward;
                } else if (order == "fwrev") {
                    return EliminationOrder::ForwardReversed;
                } else if (order == "bw") {
                    return EliminationOrder::Backward;
                } else if (order == "bwrev") {
                    return EliminationOrder::BackwardReversed;
                } else if (order == "rand") {
                    return EliminationOrder::Random;
                } else if (order == "spen") {
                    return EliminationOrder::StaticPenalty;
                } else if (order == "dpen") {
                    return EliminationOrder::DynamicPenalty;
                } else if (order == "regex") {
                    return EliminationOrder::RegularExpression;
                } else {
                    STORM_LOG_THROW(false, storm::exceptions::IllegalArgumentValueException, "Illegal elimination order selected.");
                }
            }
            
            bool EliminationSettings::isEliminateEntryStatesLastSet() const {
                return this->getOption(entryStatesLastOptionName).getHasOptionBeenSet();
            }
            
            uint_fast64_t EliminationSettings::getMaximalSccSize() const {
                return this->getOption(maximalSccSizeOptionName).getArgumentByName("maxsize").getValueAsUnsignedInteger();
            }
            
            bool EliminationSettings::isUseDedicatedModelCheckerSet() const {
                return this->getOption(useDedicatedModelCheckerOptionName).getHasOptionBeenSet();
            }

            bool EliminationSettings::isTopologicalOrderingSet() const {
                return this->getOption(topologicalOrderingOptionName).getHasOptionBeenSet();
            }

            bool EliminationSettings::isExportHTMLSet() const {
                return this->getOption(exportHTMLOptionName).getHasOptionBeenSet();
            }

            std::string EliminationSettings::getExportHTMLFilename() const {
                return this->getOption(exportHTMLOptionName).getArgumentByName("filename").getValueAsString();
            }

            bool EliminationSettings::isExportHTMLComplexitySet() const {
                return this->getOption(exportHTMLComplexityOptionName).getHasOptionBeenSet();
            }

            bool EliminationSettings::isExportHTMLStatsSet() const {
                return this->getOption(exportHTMLStatsOptionName).getHasOptionBeenSet();
            }

        } // namespace modules
    } // namespace settings
} // namespace storm
