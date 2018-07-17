/*
 *  qisequence.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <string>
#include <iostream>

#include "Util.h"
#include "Args.h"
#include "JSON.h"

#include "SequenceBase.h"
#include "AFISequence.h"
#include "CASLSequence.h"
#include "MPRAGESequence.h"
#include "MTSatSequence.h"
#include "MultiEchoSequence.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"
#include "SequenceGroup.h"

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Outputs skeleton sequence JSON objects to help users.\n"
                                "http://github.com/spinicist/QUIT");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     AFI(parser, "AFI", "AFI sequence", {'a',"afi"});
    args::Flag     CASL(parser, "CASL", "CASL sequence", {'c',"casl"});
    args::Flag     MPRAGE(parser, "MP-RAGE", "MP-RAGE sequence", {'m',"mprage"});
    args::Flag     MP2RAGE(parser, "MP2-RAGE", "MP2-RAGE sequence", {'2',"mp2rage"});
    args::Flag     MTSat(parser, "MT-Sat", "Gradient echo sequence with MT saturation pulse", {"mtsat"});
    args::Flag     MultiEcho(parser, "ME", "Multi-echo sequence with linear spacing", {"multiecho"});
    args::Flag     MultiEchoFlex(parser, "ME-FLEX", "Multi-echo sequence with flexible spacing", {"meflex"});
    args::Flag     SPGR(parser, "SPGR", "Spoiled gradient echo sequence", {'s',"spgr"});
    args::Flag     SPGREcho(parser, "SPGR-Echo", "Spoiled gradient echo sequence with echo-time", {"spgrecho"});
    args::Flag     SPGRFinite(parser, "SPGR-Finite", "Spoiled gradient echo sequence with echo-time and finite-pulse", {"spgrfinite"});
    args::Flag     SSFP(parser, "SSFP", "Balanced Steady-State Free Precession sequence", {"ssfp"});
    args::Flag     SSFPEcho(parser, "SSFP-Echo", "Balanced Steady-State Free Precession  sequence with echo-time", {"ssfpecho"});
    args::Flag     SSFPFinite(parser, "SSFP-Finite", "Balanced Steady-State Free Precession  sequence with echo-time and finite-pulse", {"ssfpfinite"});
    args::Flag     SSFPGS(parser, "SSFPGS", "Geometric solution to SSFP (no phase increment)", {"ssfpgs"});
    args::Flag     SSFPEllipse(parser, "SSFPEllipse", "Elliptic SSFP (tile phase increment and flip-angles)", {"ssfpellipse"});
    args::Flag     SSFPMT(parser, "SSFPMT", "SSFP with MT effects (specify pulse integrals etc.)", {"ssfpmt"});
    args::Flag     group(parser, "GROUP", "A group of sequences", {"group"});
    args::Flag     all(parser, "ALL", "Write all sequences", {"all"});
    QI::ParseArgs(parser, argc, argv, verbose);
    rapidjson::Document doc;
    doc.SetObject();
    QI_LOG(verbose, "Generating sequences");

    #define ADD_SEQUENCE( NAME ) \
    if ( NAME || all) {\
        QI::NAME ## Sequence sequence; \
        doc.AddMember( # NAME , sequence.toJSON(doc.GetAllocator()), doc.GetAllocator());\
        QI_LOG(verbose, "Made " #NAME " sequence");\
    }

    ADD_SEQUENCE( AFI )
    ADD_SEQUENCE( CASL )
    ADD_SEQUENCE( MPRAGE )
    ADD_SEQUENCE( MP2RAGE )
    ADD_SEQUENCE( MTSat )
    ADD_SEQUENCE( MultiEcho )
    ADD_SEQUENCE( MultiEchoFlex )
    ADD_SEQUENCE( SPGR )
    ADD_SEQUENCE( SPGREcho )
    ADD_SEQUENCE( SPGRFinite )
    ADD_SEQUENCE( SSFP )
    ADD_SEQUENCE( SSFPEcho )
    ADD_SEQUENCE( SSFPFinite )
    ADD_SEQUENCE( SSFPGS )
    ADD_SEQUENCE( SSFPEllipse )
    ADD_SEQUENCE( SSFPMT )
    if ( group || all ) { // SequenceGroup doesn't follow the naming convention
        QI::SequenceGroup group;
        group.addSequence(std::make_shared<QI::SPGRSequence>());
        group.addSequence(std::make_shared<QI::SSFPSequence>());
        doc.AddMember( "Sequences", group.toJSON(doc.GetAllocator()), doc.GetAllocator());
        QI_LOG(verbose, "Made sequence group");
    }

    QI_LOG(verbose, "Writing sequences...");
    QI::WriteJSON(std::cout, doc);
    QI_LOG(verbose, "Finished.");
    return EXIT_SUCCESS;
}

