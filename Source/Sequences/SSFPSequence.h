/*
 *  SSFP.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_SSFP_H
#define SEQUENCES_SSFP_H

#include "SequenceBase.h"
#include "fmt/format.h"

namespace QI {

struct SSFPBase : SequenceBase {
    double         TR;
    Eigen::ArrayXd FA;
    Eigen::Index   size() const override;
};

struct SSFPSequence : SSFPBase {
    Eigen::ArrayXd PhaseInc;
    QI_SEQUENCE_DECLARE(SSFP);
    Eigen::ArrayXd weights(const double f0) const override;
};
void from_json(const json &j, SSFPSequence &s);
void to_json(json &j, const SSFPSequence &s);

// struct SSFPEchoSequence : SSFPSequence {
//     QI_SEQUENCE_DECLARE(SSFPEcho);
// };

// struct SSFPFiniteSequence : SSFPBase {
//     double         Trf;
//     Eigen::ArrayXd PhaseInc;

//     QI_SEQUENCE_DECLARE(SSFPFinite);
//     Eigen::ArrayXd weights(const double f0) const override;
// };

// struct SSFPGSSequence : SSFPBase {
//     QI_SEQUENCE_DECLARE(SSFPGS);
// };

struct SSFPMTSequence : SequenceBase {
    Eigen::ArrayXd FA, TR, Trf, intB1;
    QI_SEQUENCE_DECLARE(SSFPMT);
    Eigen::Index size() const override;
};
void from_json(const json &j, SSFPMTSequence &s);
void to_json(json &j, const SSFPMTSequence &s);

} // End namespace QI

namespace fmt {
template <> struct formatter<QI::SSFPSequence> {
    template <typename ParseContext> constexpr auto parse(ParseContext &ctx) { return ctx.begin(); }

    template <typename FormatContext> auto format(const QI::SSFPSequence &s, FormatContext &ctx) {
        return format_to(ctx.out(),
                         "SSFP:\n\tTR: {}\n\tFA: {}\n\tPhaseInc: {}",
                         s.TR,
                         (s.FA * 180. / M_PI).transpose(),
                         (s.PhaseInc * 180. / M_PI).transpose());
    }
};
} // namespace fmt

#endif // SEQUENCES_SSFP_H
