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
#include <cereal/archives/json.hpp>
#include "Models.h"

namespace QI {

struct SequenceBase {
    virtual std::string &name() const = 0;
    virtual Eigen::Index size() const = 0;
    virtual Eigen::ArrayXcd signal(const std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const = 0;
    virtual void load(cereal::JSONInputArchive &ar) = 0;
    virtual void save(cereal::JSONOutputArchive &ar) const = 0;

    virtual size_t count() const;
    virtual Eigen::ArrayXd weights(double f0 = 0.0) const;
    virtual Eigen::ArrayXd  signal_magnitude(const std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const;
};

#define QI_SEQUENCE_DECLARE( N ) \
    std::string &name() const override { static std::string name = #N; return name; }\
    Eigen::ArrayXcd signal(std::shared_ptr<QI::Model::ModelBase> m, const Eigen::VectorXd &par) const override;\
    void load(cereal::JSONInputArchive &ar) override;\
    void save(cereal::JSONOutputArchive &ar) const override;

#define QI_SEQUENCE_DECLARE_NOSIG( N ) \
    std::string &name() const override { static std::string name = #N; return name; }\
    void load(cereal::JSONInputArchive &ar) override;\
    void save(cereal::JSONOutputArchive &ar) const override;

} // End namespace QI

#endif // SEQUENCES_BASE_H
