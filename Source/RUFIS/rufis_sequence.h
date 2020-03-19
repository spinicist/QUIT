#pragma once

#include "SequenceBase.h"
#include "rufis_pulse.h"
#include <unordered_map>

struct RUFISSequence : QI::SequenceBase {
    double                                     TR, Tramp;
    Eigen::ArrayXd                             FA, Trf;
    Eigen::ArrayXi                             groups_per_seg;
    int                                        spokes_per_seg;
    std::unordered_map<std::string, PrepPulse> prep_pulses;
    std::vector<std::string>                   prep;
    QI_SEQUENCE_DECLARE(RUFIS);
    Eigen::Index size() const override { return prep.size(); };
};
void from_json(const json &j, RUFISSequence &s);
