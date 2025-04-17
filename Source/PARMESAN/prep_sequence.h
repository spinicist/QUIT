#pragma once

#include "SequenceBase.h"

struct PrepSequence : QI::SequenceBase {
    double          TR, Tramp, Trf, Dprep, Dseg;
    int             SPS, spoilers;
    Eigen::ArrayXd  FA, FAprep, Tprep, fprep;
    QI_SEQUENCE_DECLARE(Prep);
    auto preps() const -> Eigen::Index;
    Eigen::Index size() const override;
};
void from_json(const json &j, PrepSequence &s);
