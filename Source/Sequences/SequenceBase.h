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
#include <Eigen/Dense>

#include "Util.h"
#include "Signals/SignalEquations.h"
#include "Models/Models.h"

using namespace std;
using namespace Eigen;

class SequenceBase {
    protected:
        double m_TR = 0.;
        ArrayXd m_flip;

    public:
        virtual ArrayXcd signal(const shared_ptr<Model> m, const VectorXd &p) const = 0;
        virtual size_t size() const = 0;
        virtual void write(ostream &os) const = 0;
        virtual string name() const = 0;
        virtual size_t count() const { return 1; }
        double TR() const { return m_TR; }
        void setTR(const double TR) { m_TR = TR; }
        const ArrayXd & flip() const { return m_flip; }
        void setFlip(const ArrayXd &f) { m_flip = f; }
        virtual ArrayXd weights(double f0 = 0.0) const { return ArrayXd::Ones(size()); }        
};
ostream& operator<<(ostream& os, const SequenceBase& s);


class SequenceGroup : public SequenceBase {
private:
    vector<shared_ptr<SequenceBase>> m_sequences;

public:
    SequenceGroup();
    void write(ostream &os) const override;
    string name() const override { return "Sequences"; }

    size_t count() const override;
    shared_ptr<SequenceBase> sequence(const size_t i) const;
    vector<shared_ptr<SequenceBase>> &sequences();

    size_t size() const override;
    ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
    ArrayXd weights(const double f0 = 0.0) const override;
    
    void addSequence(const shared_ptr<SequenceBase> &seq);
};

#endif // SEQUENCES_BASE_H