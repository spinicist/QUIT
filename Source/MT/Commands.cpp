#include "Commands.h"

void add_mt_commands(std::map<std::string, std::function<int(int, char **)>> &commands) {
    commands["lineshape"]    = &lineshape_main;
    commands["lorentzian"]   = &lorentzian_main;
    commands["mtr"]          = &mtr_main;
    commands["mtsat"]        = &mtsat_main;
    commands["qmt"]          = &qmt_main;
    commands["ssfp_emt"]     = &ssfp_emt_main;
    commands["zspec_b1"]     = &zspec_b1_main;
    commands["zspec_interp"] = &zspec_interp_main;
}