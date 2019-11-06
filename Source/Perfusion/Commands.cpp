#include "Commands.h"

void add_perfusion_commands(std::map<std::string, std::function<int(int, char **)>> &commands) {
    commands["ase_oef"] = &ase_oef_main;
    commands["asl"]     = &asl_main;
    commands["zshim"]   = &zshim_main;
}