#pragma once

#include "SequenceBase.h"

struct PrepZTESequence : QI::SequenceBase {
    double          TR, Tramp, Trf, Tprep, Dprep, Dseg;
    int             SPS, spoilers;
    Eigen::ArrayXd  FA, FAprep, fprep;
    Eigen::MatrixXd basis;
    QI_SEQUENCE_DECLARE(PrepZTE);
    Eigen::Index preps() const { return FAprep.size(); }
    Eigen::Index size() const override {
        if (basis.size()) {
            return basis.rows();
        } else {
            return SPS * FAprep.size();
        }
    };
};
void from_json(const json &j, PrepZTESequence &s);
