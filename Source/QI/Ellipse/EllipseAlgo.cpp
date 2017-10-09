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

#include "QI/Ellipse/EllipseAlgo.h"
/*#include "QI/Ellipse/EllipseHelpers.h"
#include "QI/Banding.h"
#include "QI/GoldenSection.h"*/

namespace QI {

EllipseAlgo::EllipseAlgo(std::shared_ptr<QI::SSFPEcho> &seq, bool debug) :
    m_sequence(seq), m_debug(debug)
{
    m_zero = TOutput(m_sequence->flip().rows());
    m_zero.Fill(0.);
}

bool EllipseAlgo::apply(const std::vector<TInput> &inputs,
                        const std::vector<TConst> &consts,
                        std::vector<TOutput> &outputs, TConst &residual,
                        TInput &resids, TIters &its) const
{
    const int np = m_sequence->phase_incs().rows();
    const double B1 = consts[0];
    for (int f = 0; f < m_sequence->flip().rows(); f++) {
        Eigen::ArrayXcf data(np);
        for (int i = 0; i < np; i++) {
            data[i] = inputs[0][f*np + i];
        }
        if (m_debug) {
            std::cout << "Flip: " << m_sequence->flip()[f] << " Data: " << data.transpose() << std::endl;
        }
        Eigen::ArrayXd tempOutputs = this->apply_internal(data, B1 * m_sequence->flip()[f], m_sequence->TR(), m_sequence->phase_incs(), m_debug, residual);
        if (m_debug) {
            std::cout << "Outputs: " << tempOutputs.transpose() << std::endl;
        }
        for (int o = 0; o < this->numOutputs(); o++) {
            outputs[o][f] = tempOutputs[o];
        }
    }
    return true;
}

} // End namespace QI