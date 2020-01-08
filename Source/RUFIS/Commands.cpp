#include "Commands.h"

void add_rufis_commands(std::map<std::string, std::function<int(int, char **)>> &commands) {
    commands["mupa"]     = &mupa_main;
    commands["rf-sim"]   = &rf_sim_main;
    commands["rufis-mt"] = &rufis_mt_main;
}