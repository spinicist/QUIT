#include "Commands.h"

void add_b1_commands(std::map<std::string, std::function<int(int, char **)>> &commands) {
    commands["afi"]     = &afi_main;
    commands["dream"]   = &dream_main;
    commands["b1_papp"] = &b1_papp_main;
}