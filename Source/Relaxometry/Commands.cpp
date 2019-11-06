#include "Commands.h"

void add_relax_commands(std::map<std::string, std::function<int(int, char **)>> &commands) {
    commands["jsr"]          = &jsr_main;
    commands["mpm_r2s"]      = &mpm_r2s_main;
    commands["planet"]       = &planet_main;
    commands["ssfp_ellipse"] = &ssfp_ellipse_main;
    commands["despot1"]      = &despot1_main;
    commands["despot1hifi"]  = &despot1hifi_main;
    commands["despot2"]      = &despot2_main;
    commands["despot2fm"]    = &despot2fm_main;
    // commands["mcdespot"]     = &mcdespot_main;
    commands["mp2rage"]   = &mp2rage_main;
    commands["multiecho"] = &multiecho_main;
}