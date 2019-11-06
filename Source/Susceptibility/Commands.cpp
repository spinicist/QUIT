#include "Commands.h"

void add_suscep_commands(std::map<std::string, std::function<int(int, char **)>> &commands) {
    commands["fieldmap"]       = &fieldmap_main;
    commands["unwrap_laplace"] = &unwrap_laplace_main;
    commands["unwrap_path"]    = &unwrap_path_main;
}