/*
 *  Args.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Args.h"

#include <iostream>
#include <sstream>
#include <iomanip>

namespace QI {

std::ostream &operator<< (std::ostream &os, const TOption &o) {
    if (o.short_name != '\0') {
        os << "  -" << o.short_name << ",";
    } else {
        os << "     ";
    }
    os << "--"
       << std::setw(10) << std::left << o.long_name
       << std::setw(0) << " : " << o.usage;
    return os;
}

ArgParser::ArgParser(int argc, char **argv, const std::string &usage, 
                     const std::vector<TOption> &opts) :
    m_opts(opts)
{
    int optind = 1;
    while (optind < argc) {
        std::string thisopt(argv[optind]);
        auto it = m_opts.end();
        std::string arg = "";
        if (thisopt == "--") { // Premature end of TOptions
            optind++;
            while (optind < argc) {
                thisopt = argv[optind++];
                m_nonopts.push_back(thisopt);
            }
        } else if (thisopt[0] != '-') { // Non-TOption
            m_nonopts.push_back(thisopt);
        } else if (thisopt[1] == '-') { // Long TOption
            size_t pos_end = thisopt.find("=");
            if (pos_end != std::string::npos) {
                arg = thisopt.substr(pos_end + 1);
                thisopt = thisopt.substr(2, pos_end - 2);
            } else {
                thisopt.erase(0, 2);
            }
            it = std::find_if(m_opts.begin(), m_opts.end(), [&] (const TOption &o) { return o.long_name == thisopt; });
            if (it == m_opts.end()) {
                QI_EXCEPTION("Unknown long option '" << thisopt << "' given on command-line.");
            }
            m_args.push_back({it->long_name, arg});
        } else { // Short TOption
            char sopt = thisopt[1];
            it = std::find_if(m_opts.begin(), m_opts.end(), [&] (const TOption &o) { return o.short_name == sopt; });
            if (it == m_opts.end()) {
                QI_EXCEPTION("Unknown short option '" << std::string(thisopt) << "' given on command-line.");
            }
            if (it->has_arg) {
                if (thisopt.size() > 2) { // Handle ArgParser without spaces
                    arg = thisopt.substr(2);
                } else if ((optind+1) == argc) {
                    QI_EXCEPTION("Missing required argument '" << std::to_string(sopt) << "'");
                } else {
                    arg = argv[++optind];
                }
            }
            m_args.push_back({it->long_name, arg});
        }
        optind++;
    }
    if (option_present("help")) {
        std::cout << usage << std::endl << std::endl;
        for (auto &o : m_opts) {
            std::cout << o << std::endl;
        }
        exit(EXIT_SUCCESS);
    }
}

bool ArgParser::option_present(const std::string &name) {
    auto it = std::find_if(m_args.cbegin(), m_args.cend(), [&] (const std::pair<std::string, std::string> &p) { return p.first == name; });
    if (it == m_args.cend()) {
        return false;
    } else {
        const std::string value = it->second;
        return true;
    }
}

std::string ArgParser::string_value(const std::string &name,
                                    const std::string &def_value)
{
    auto it = std::find_if(m_args.cbegin(), m_args.cend(), [&] (const std::pair<std::string, std::string> &p) { return p.first == name; });
    if (it == m_args.cend()) {
        return def_value;
    } else {
        return it->second;
    }
}

}