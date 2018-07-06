/*
 *  SequenceBase.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "SequenceBase.h"
#include "AFISequence.h"
#include "CASLSequence.h"
#include "MPRAGESequence.h"
#include "MTSatSequence.h"
#include "MultiEchoSequence.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"

using namespace std::string_literals;

namespace QI {

std::shared_ptr<SequenceBase> SequenceFromJSON(rapidjson::Value &json) {
    rapidjson::Value::MemberIterator seq = json.MemberBegin();
    #define CHECK_MAKE( NAME ) \
        (seq->name.GetString() == #NAME ## s ) { return std::make_shared<QI::NAME ## Sequence>(seq->value); }
    if CHECK_MAKE( SPGR )
    else if CHECK_MAKE( SPGREcho )
    else if CHECK_MAKE( SPGRFinite )
    else if CHECK_MAKE( MPRAGE )
    else if CHECK_MAKE( SSFP )
    else if CHECK_MAKE( SSFPEcho )
    else if CHECK_MAKE( SSFPFinite )
    else if CHECK_MAKE( SSFPGS )
    else if CHECK_MAKE( SSFPMT )
    else if CHECK_MAKE( MTSat )
    else if CHECK_MAKE( MultiEcho )
    else if CHECK_MAKE( MultiEchoFlex )
    else if CHECK_MAKE( CASL )
    else { QI_FAIL("Unimplemented JSON load for sequence type: " << seq->name.GetString()); }
    #undef CHECK_MAKE
}

rapidjson::Value JSONFromSequence(const std::shared_ptr<SequenceBase> &s, rapidjson::Document::AllocatorType &a) {
    rapidjson::Value json(rapidjson::kObjectType);
    #define MAKE_JSON( NAME ) \
        (s->name() == #NAME ) { json.AddMember(#NAME, std::static_pointer_cast< QI::NAME ## Sequence >(s)->toJSON(a), a); }
    if MAKE_JSON( SPGR )
    else if MAKE_JSON( SPGREcho )
    else if MAKE_JSON( SPGRFinite )
    else if MAKE_JSON( MPRAGE )
    else if MAKE_JSON( SSFP )
    else if MAKE_JSON( SSFPEcho )
    else if MAKE_JSON( SSFPFinite )
    else if MAKE_JSON( SSFPGS )
    else if MAKE_JSON( SSFPMT )
    else if MAKE_JSON( MTSat )
    else if MAKE_JSON( MultiEcho )
    else if MAKE_JSON( MultiEchoFlex )
    else if MAKE_JSON( CASL )
    else { QI_FAIL("Unimplemented save for sequence type: " << s->name()); }
    #undef MAKE_JSON
    return json;
}

size_t SequenceBase::count() const{
    return 1;
}

Eigen::ArrayXd SequenceBase::weights(const double /* Unused */) const {
    return Eigen::ArrayXd::Ones(size()); // Default weights are constant
}

Eigen::ArrayXd SequenceBase::signal_magnitude(const std::shared_ptr<Model::ModelBase> m, const Eigen::VectorXd &p) const {
    Eigen::ArrayXcd c_signal = this->signal(m, p);
    return c_signal.abs();
}

} // End namespace QI
