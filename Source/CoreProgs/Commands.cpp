#include "Commands.h"

void add_core_commands(std::map<std::string, std::function<int(int, char **)>> &commands) {
    commands["diff"]     = &diff_main;
    commands["hdr"]      = &hdr_main;
    commands["newimage"] = &newimage_main;
}