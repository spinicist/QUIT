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
#include <vector>
#include <memory>
#include <Eigen/Core>
#include "Models.h"
#include <cereal/archives/json.hpp>
#include <cereal/types/polymorphic.hpp>
#include "EigenCereal.h"

namespace QI {

struct SequenceBase {
    virtual Eigen::ArrayXcd signal(const std::shared_ptr<Model> m, const Eigen::VectorXd &p) const = 0;
    virtual Eigen::ArrayXd  signal_magnitude(const std::shared_ptr<Model> m, const Eigen::VectorXd &p) const;
    virtual size_t size() const = 0;
    virtual size_t count() const;
    virtual Eigen::ArrayXd weights(double f0 = 0.0) const;
};

template<typename Sequence>
Sequence ReadSequence(std::istream &is, std::string name, bool verbose) {
    cereal::JSONInputArchive in_archive(is);
    Sequence sequence;
    in_archive(sequence);

    if (verbose) {
        std::cout << "Read " << name << ": " << std::endl;
        cereal::JSONOutputArchive archive(std::cout);
        archive(sequence);
    }
    
    return sequence;
}

} // End namespace QI

//CEREAL_FORCE_DYNAMIC_INIT(SequenceBase);

#endif // SEQUENCES_BASE_H
