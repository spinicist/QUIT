#pragma once

#include "SequenceBase.h"

struct PrepSequence : QI::SequenceBase {
    double         TR, Tramp, Trf;
    int            SPS, spoilers;
    Eigen::ArrayXd FA, FAprep, Tprep, Tpreseg, Tpostseg, fprep;
    QI_SEQUENCE_DECLARE(Prep);
    auto         preps() const -> Eigen::Index;
    Eigen::Index size() const override;
    PrepSequence(double         TR,
                 double         Tramp,
                 double         Trf,
                 int            SPS,
                 int            spoilers,
                 Eigen::ArrayXd FA,
                 Eigen::ArrayXd FAprep,
                 Eigen::ArrayXd Tprep,
                 Eigen::ArrayXd Tpreseg_,
                 Eigen::ArrayXd Tpostseg_,
                 Eigen::ArrayXd fprep);
};
void from_json(const json &j, PrepSequence &s);
void to_json(json &j, const PrepSequence &s);