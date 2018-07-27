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
#include "JSON.h"

namespace QI {

struct SequenceBase {
    virtual std::string &name() const = 0;
    virtual Eigen::Index size() const = 0;
    virtual rapidjson::Value toJSON(rapidjson::Document::AllocatorType &) const = 0;

    virtual size_t count() const;
    virtual Eigen::ArrayXd weights(double f0 = 0.0) const;
};
SequenceBase *SequenceFromJSON(const rapidjson::Value &);
rapidjson::Value JSONFromSequence(const SequenceBase *, rapidjson::Document::AllocatorType &);

#define QI_SEQUENCE_DECLARE( NAME ) \
    std::string &name() const override { static std::string name = #NAME; return name; }\
    NAME ## Sequence (const rapidjson::Value &);\
    NAME ## Sequence () = default;\
    rapidjson::Value toJSON(rapidjson::Document::AllocatorType &) const override;

#define QI_SEQUENCE_DECLARE_NOSIG( NAME ) \
    std::string &name() const override { static std::string name = #NAME; return name; }\
    NAME ## Sequence (const rapidjson::Value &);\
    NAME ## Sequence () = default;\
    rapidjson::Value toJSON(rapidjson::Document::AllocatorType &) const override;

} // End namespace QI

#endif // SEQUENCES_BASE_H
