#pragma once

#include "JSON.h"
#include <Eigen/Dense>

struct RFPulse {
    Eigen::ArrayXd B1x, B1y, timestep;
};
void from_json(const json &j, RFPulse &p);
void to_json(json &j, RFPulse const &p);

struct MTPulse {
    Eigen::ArrayXd B1x, B1y;
    double         FA, width;
};
void from_json(const json &j, MTPulse &p);
void to_json(json &j, MTPulse const &p);

struct PrepPulse {
    double FAeff, int_b1_sq, T_long, T_trans;
};
void from_json(const json &j, PrepPulse &p);
void to_json(json &j, const PrepPulse &s);