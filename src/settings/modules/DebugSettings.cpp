#include "src/settings/modules/DebugSettings.h"

#include "src/settings/SettingsManager.h"
#include "src/settings/Option.h"
#include "src/settings/OptionBuilder.h"
#include "src/settings/Argument.h"
#include "src/settings/ArgumentBuilder.h"

namespace storm {
    namespace settings {
        namespace modules {
            
            const std::string DebugSettings::moduleName = "debug";
            const std::string DebugSettings::debugOptionName = "debug";
            const std::string DebugSettings::traceOptionName = "trace";
            const std::string DebugSettings::logfileOptionName = "logfile";
            const std::string DebugSettings::logfileOptionShortName = "l";
            const std::string DebugSettings::testOptionName = "test";
 
            DebugSettings::DebugSettings() : ModuleSettings(moduleName) {
                this->addOption(storm::settings::OptionBuilder(moduleName, debugOptionName, false, "Print debug output.").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, traceOptionName, false, "Print even more debug output.").build());
                this->addOption(storm::settings::OptionBuilder(moduleName, logfileOptionName, false, "If specified, the log output will also be written to this file.").setShortName(logfileOptionShortName)
                                .addArgument(storm::settings::ArgumentBuilder::createStringArgument("filename", "The name of the file to write the log.").build()).build());
                this->addOption(storm::settings::OptionBuilder(moduleName, testOptionName, false, "Activate a test setting.").build());
            }
            
            bool DebugSettings::isDebugSet() const {
                return this->getOption(debugOptionName).getHasOptionBeenSet();
            }
            
            bool DebugSettings::isTraceSet() const {
                return this->getOption(traceOptionName).getHasOptionBeenSet();
            }
            
            bool DebugSettings::isLogfileSet() const {
                return this->getOption(logfileOptionName).getHasOptionBeenSet();
            }
            
            std::string DebugSettings::getLogfilename() const {
                return this->getOption(traceOptionName).getArgumentByName("filename").getValueAsString();
            }
            
            bool DebugSettings::isTestSet() const {
                return this->getOption(testOptionName).getHasOptionBeenSet();
            }
            
        } // namespace modules
    } // namespace settings
} // namespace storm