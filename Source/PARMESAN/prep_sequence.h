#pragma once

#include "SequenceBase.h"

struct PrepSequence : QI::SequenceBase {
    double         TR, Tramp, Trf, Dprep, Dseg;
    int            SPS, spoilers;
    Eigen::ArrayXd FA, FAprep, Tprep, fprep;
    QI_SEQUENCE_DECLARE(Prep);
    auto         preps() const -> Eigen::Index;
    Eigen::Index size() const override;
    PrepSequence(double         TR,
                 double         Tramp,
                 double         Trf,
                 double         Dprep,
                 double         Dseg,
                 int            SPS,
                 int            spoilers,
                 Eigen::ArrayXd FA,
                 Eigen::ArrayXd FAprep,
                 Eigen::ArrayXd Tprep,
                 Eigen::ArrayXd fprep);
};
void from_json(const json &j, PrepSequence &s);
