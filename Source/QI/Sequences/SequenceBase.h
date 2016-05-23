/*
 *  SequenceBase.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_BASE_H
#define SEQUENCES_BASE_H

#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <Eigen/Core>

#include "QI/Macro.h"
#include "QI/Util.h"
#include "QI/Signals/SignalEquations.h"
#include "QI/Models/Models.h"

namespace QI {

class SequenceBase {
    protected:
        double m_TR = 0.;
        Eigen::ArrayXd m_flip;

    public:
        SequenceBase();
        SequenceBase(const Eigen::ArrayXd &flip, const double TR);
    
        virtual Eigen::ArrayXcd signal(const std::shared_ptr<Model> m, const Eigen::VectorXd &p) const = 0;
        virtual size_t size() const = 0;
        virtual void write(std::ostream &os) const = 0;
        virtual std::string name() const = 0;
        virtual size_t count() const { return 1; }
        double TR() const { return m_TR; }
        void setTR(const double TR) { m_TR = TR; }
        const Eigen::ArrayXd & flip() const { return m_flip; }
        void setFlip(const Eigen::ArrayXd &f) { m_flip = f; }
        virtual Eigen::ArrayXd weights(double f0 = 0.0) const { return Eigen::ArrayXd::Ones(size()); }        
};
std::ostream& operator<<(std::ostream& os, const SequenceBase& s);


class SequenceGroup : public SequenceBase {
private:
    std::vector<std::shared_ptr<SequenceBase>> m_sequences;

public:
    SequenceGroup();
    void write(std::ostream &os) const override;
    std::string name() const override { return "Sequences"; }

    size_t count() const override;
    std::shared_ptr<SequenceBase> sequence(const size_t i) const;
    std::vector<std::shared_ptr<SequenceBase>> &sequences();

    size_t size() const override;
    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
    Eigen::ArrayXd weights(const double f0 = 0.0) const override;
    
    void addSequence(const std::shared_ptr<SequenceBase> &seq);
};

} // End namespace QI

#endif // SEQUENCES_BASE_H