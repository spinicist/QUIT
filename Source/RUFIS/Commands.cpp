#include "Commands.h"

void add_rufis_commands(std::map<std::string, std::function<int(int, char **)>> &commands) {
    commands["rufis-mupa"] = &rufis_mupa_main;
    commands["rf-sim"]     = &rf_sim_main;
    commands["rufis-ss"]   = &rufis_ss_main;
}