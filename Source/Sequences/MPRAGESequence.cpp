/*
 *  MPRAGESequence.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "MPRAGESequence.h"

namespace QI
{

Eigen::Index MPRAGESequence::size() const { return 1; }

MPRAGESequence::MPRAGESequence(const rapidjson::Value &json)
{
    if (json.IsNull())
        QI_FAIL("Could not read sequence: " << name());
    TR = GetMember(json, "TR").GetDouble();
    TI = GetMember(json, "TI").GetDouble();
    TD = GetMember(json, "TD").GetDouble();
    k0 = GetMember(json, "k0").GetInt();
    FA = GetMember(json, "FA").GetDouble() * M_PI / 180;
    ETL = GetMember(json, "ETL").GetInt();
    eta = GetMember(json, "eta").GetDouble();
}

rapidjson::Value
MPRAGESequence::toJSON(rapidjson::Document::AllocatorType &a) const
{
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("TI", TI, a);
    json.AddMember("TD", TD, a);
    json.AddMember("eta", eta, a);
    json.AddMember("FA", FA * 180 / M_PI, a);
    json.AddMember("ETL", ETL, a);
    json.AddMember("k0", k0, a);
    return json;
}

/*
 * MP2RAGE
 */

Eigen::Index MP2RAGESequence::size() const { return 2; }

MP2RAGESequence::MP2RAGESequence(const rapidjson::Value &json)
{
    if (json.IsNull())
        QI_FAIL("Could not read sequence: " << name());
    TR = GetMember(json, "TR").GetDouble();
    auto TI = ArrayFromJSON(json, "TI");
    auto TRPrep = GetMember(json, "TRPrep").GetDouble();
    FA = ArrayFromJSON(json, "FA", M_PI / 180);
    SegLength = GetMember(json, "SegLength").GetInt();
    k0 = GetMember(json, "k0").GetInt();
    TD[0] = TI[0] - k0*TR;
    TD[1] = TI[1] - (TI[0] + SegLength*TR);
    TD[2] = TRPrep - (TI[1] + (SegLength - k0)*TR);
}

rapidjson::Value
MP2RAGESequence::toJSON(rapidjson::Document::AllocatorType &a) const
{
    Eigen::Array2d TI{TD[0], TD[1] + SegLength*TR + TD[0]};
    double TRPrep = TI[1] + (SegLength - k0) * TR + TD[2];
    rapidjson::Value json;
    json.AddMember("TR", TR, a);
    json.AddMember("TRPrep", TRPrep, a);
    json.AddMember("TI", ArrayToJSON(TI, a), a);
    json.AddMember("FA", ArrayToJSON(FA, a, 180 / M_PI), a);
    json.AddMember("SegLength", SegLength, a);
    json.AddMember("k0", k0, a);
    return json;
}

} // End namespace QI
