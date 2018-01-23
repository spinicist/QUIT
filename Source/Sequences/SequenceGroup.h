/*
 *  SequenceGroup.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef SEQUENCES_GROUP_H
#define SEQUENCES_GROUP_H

#include <cereal/types/vector.hpp>
#include "SequenceBase.h"
#include "SPGRSequence.h"

namespace QI {

struct SequenceWrapper {
    enum class Tag {SPGR, SPGREcho, SPGRFinite};
    Tag tag;
    union {
        SPGRSequence SPGR;
        SPGREchoSequence SPGREcho;
        SPGRFiniteSequence SPGRFinite;
    };

    SequenceWrapper() :
        tag(Tag::SPGR), SPGR()
    { }

    SequenceWrapper(const SequenceWrapper &sw) :
        tag(Tag::SPGR), SPGR()
    {
        tag = sw.tag;
        switch (tag) {
            case Tag::SPGR:       SPGR = sw.SPGR; break;
            case Tag::SPGREcho:   SPGREcho = sw.SPGREcho; break;
            case Tag::SPGRFinite: SPGRFinite = sw.SPGRFinite; break;
        }
    }

    ~SequenceWrapper() {
        switch (tag) {
            case Tag::SPGR: SPGR.~SPGRSequence(); break;
            case Tag::SPGREcho: SPGREcho.~SPGREchoSequence(); break;
            case Tag::SPGRFinite: SPGRFinite.~SPGRFiniteSequence(); break;
        }
    }

    SequenceWrapper(const SPGRSequence &s) : tag(Tag::SPGR), SPGR(s) { }
    SequenceWrapper(const SPGREchoSequence &s) : tag(Tag::SPGREcho), SPGREcho(s) { }
    SequenceWrapper(const SPGRFiniteSequence &s) : tag(Tag::SPGRFinite), SPGRFinite(s) { }

    SequenceBase *ptr() {
        SequenceBase *p;
        switch (tag) {
            case Tag::SPGR:       p = &SPGR; break;
            case Tag::SPGREcho:   p = &SPGREcho; break;
            case Tag::SPGRFinite: p = &SPGRFinite; break;
        }
        return p;
    }

    size_t size() const {
        switch (tag) {
            case Tag::SPGR:       return SPGR.size(); break;
            case Tag::SPGREcho:   return SPGREcho.size(); break;
            case Tag::SPGRFinite: return SPGRFinite.size(); break;
        }
    }

    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &p) const {
        switch (tag) {
            case Tag::SPGR:       return SPGR.signal(m, p); break;
            case Tag::SPGREcho:   return SPGREcho.signal(m, p); break;
            case Tag::SPGRFinite: return SPGRFinite.signal(m, p); break;
        }
    }

    Eigen::ArrayXd weights(const double f0 = 0.0) const {
        switch (tag) {
            case Tag::SPGR:       return SPGR.weights(f0); break;
            case Tag::SPGREcho:   return SPGREcho.weights(f0); break;
            case Tag::SPGRFinite: return SPGRFinite.weights(f0); break;
        }
    }

    template<typename Archive> void serialize(Archive &ar) {
        ar(cereal::make_nvp("sequence", tag));
        switch (tag) {
            case Tag::SPGR:       ar(CEREAL_NVP(SPGR)); break;
            case Tag::SPGREcho:   ar(CEREAL_NVP(SPGREcho)); break;
            case Tag::SPGRFinite: ar(CEREAL_NVP(SPGRFinite)); break;
        }
    }
};

struct SequenceGroup : SequenceBase {
    std::vector<SequenceWrapper> sequences;
    SequenceWrapper &operator[](const size_t i);

    size_t count() const override;
    size_t size() const override;
    Eigen::ArrayXcd signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const override;
    Eigen::ArrayXd weights(const double f0 = 0.0) const override;

    void addSequence(const SequenceWrapper &w);
    template<typename Archive>
    void serialize(Archive &archive) {
        archive(CEREAL_NVP(sequences));
    }
};

} // End namespace QI

#endif // SEQUENCES_BASE_H
