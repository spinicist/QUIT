#pragma once

#include "SequenceBase.h"
#include "rufis_pulse.h"
#include <unordered_map>

struct MTSequence : QI::SequenceBase {
    double         TR, Tramp, Trf, Tspoil, MT_pulsewidth;
    int            SPS;
    Eigen::ArrayXd RUFIS_FA, MT_FA, MT_offsets;

    MTPulse MT_pulse;
    QI_SEQUENCE_DECLARE(MT);
    Eigen::Index size() const override { return RUFIS_FA.rows(); };
};
void from_json(const json &j, MTSequence &s);
void to_json(json &j, MTSequence const &s);