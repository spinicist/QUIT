#pragma once

#include "JSON.h"

struct SSSequence {
    double         TR, Tramp;
    Eigen::ArrayXd FA, Trf;
    Eigen::ArrayXi groups_per_seg;
    int            spokes_per_seg;
    double         prep_p1, prep_p2, prep_Trf;
    Eigen::ArrayXd prep_FA, prep_df;
    Eigen::Index   size() const;
};
void from_json(const json &j, SSSequence &s);