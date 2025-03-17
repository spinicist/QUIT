#pragma once

#include "SequenceBase.h"
#include "rf_pulse.h"
#include <unordered_map>

struct PrepSequence : QI::SequenceBase {
    double                                     TR, Tramp, Trf, Tprep;
    Eigen::ArrayXd                             FA, FAprep;
    int                                        SPS, spoilers;
    QI_SEQUENCE_DECLARE(PrepSequence);
    Eigen::Index preps() const { return FAprep.size(); };
    Eigen::Index size() const override { return preps() * SPS; };
};
void from_json(const json &j, RUFISSequence &s);
