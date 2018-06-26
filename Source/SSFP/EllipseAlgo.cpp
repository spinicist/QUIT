/*
 *  EllipseAlgo.cpp
 *
 *  Copyright (c) 2016, 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "EllipseAlgo.h"
/*#include "EllipseHelpers.h"
#include "Banding.h"
#include "GoldenSection.h"*/

namespace QI {

EllipseAlgo::EllipseAlgo(const QI::SSFPEllipseSequence &seq, bool debug) :
    m_debug(debug), m_seq(seq)
{
    m_zero = TOutput(m_seq.size());
    m_zero.Fill(0.);
}

EllipseAlgo::TStatus EllipseAlgo::apply(const std::vector<TInput> &inputs,
                        const std::vector<TConst> &consts,
                        const TIndex &, // Unused
                        std::vector<TOutput> &outputs, TOutput &residual,
                        TInput & /* Unused */, TIterations & /* Unused */) const
{
    const int np = m_seq.PhaseInc.rows();
    const double B1 = consts[0];
    if (m_debug) std::cout << "Input: " << inputs[0] << std::endl;
    for (int f = 0; f < m_seq.FA.rows(); f++) {
        Eigen::ArrayXcf data(np);
        for (int i = 0; i < np; i++) {
            data[i] = inputs[0][f*np + i];
        }
        if (m_debug) {
            std::cout << "Flip: " << m_seq.FA[f] << " Data: " << data.transpose() << std::endl;
        }
        Eigen::ArrayXd tempOutputs = this->apply_internal(data, B1 * m_seq.FA[f], m_seq.TR, m_seq.PhaseInc, m_debug, residual[f]);
        if (m_debug) {
            std::cout << "Outputs: " << tempOutputs.transpose() << std::endl;
        }
        for (size_t o = 0; o < this->numOutputs(); o++) {
            outputs[o][f] = tempOutputs[o];
        }
    }
    return std::make_tuple(true, "");
}

} // End namespace QI