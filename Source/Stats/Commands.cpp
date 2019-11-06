#include "Commands.h"

void add_stats_commands(std::map<std::string, std::function<int(int, char **)>> &commands) {
    commands["glm_contrasts"] = &glm_contrasts_main;
    commands["glm_setup"]     = &glm_setup_main;
    commands["rois"]          = &rois_main;
}