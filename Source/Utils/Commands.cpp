#include "Commands.h"

void add_utils(std::map<std::string, std::function<int(int, char **)>> &commands) {
    commands["affine"]       = &affine_main;
    commands["affine-angle"] = &affine_angle_main;
    commands["coil-combine"] = &coil_combine_main;
    commands["gradient"]     = &gradient_main;
    commands["pca"]          = &pca_main;
    commands["rfprofile"]    = &rfprofile_main;
    commands["select"]       = &select_main;
    commands["ssfp_bands"]   = &ssfp_bands_main;
    commands["tvmask"]       = &tvmask_main;
    commands["complex"]      = &complex_main;
    commands["kfilter"]      = &kfilter_main;
    commands["mask"]         = &mask_main;
    commands["polyfit"]      = &polyfit_main;
    commands["polyimg"]      = &polyimg_main;
}