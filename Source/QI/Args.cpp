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

ArgParser::ArgParser(int argc, char **argv, std::vector<const TOption> opts) :
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
                QI_EXCEPTION("Unhandled long TOption: " + thisopt);
            }
            m_found.push_back({it->long_name, arg});
        } else { // Short TOption
            char sopt = thisopt[1];
            it = std::find_if(m_opts.begin(), m_opts.end(), [&] (const TOption &o) { return o.short_name == sopt; });
            if (it == m_opts.end()) {
                QI_EXCEPTION("Unhandled short TOption: " + std::string(thisopt));
            }
            if (it->has_arg) {
                if (thisopt.size() > 2) { // Handle ArgParser without spaces
                    arg = thisopt.substr(2);
                } else if ((optind+1) == argc) {
                    QI_EXCEPTION("Missing required argument for TOption " + std::to_string(sopt));
                } else {
                    arg = argv[++optind];
                }
            }
            m_found.push_back({it->long_name, arg});
        }
        optind++;
    }
}

bool ArgParser::found(const std::string &name) {
    auto it = std::find_if(m_found.cbegin(), m_found.cend(), [&] (const std::pair<std::string, std::string> &p) { return p.first == name; });
    if (it == m_found.cend()) {
        return false;
    } else {
        const std::string value = it->second;
        return true;
    }
}

std::pair<bool, const std::string> ArgParser::consume(const std::string &name) {
    auto it = std::find_if(m_found.cbegin(), m_found.cend(), [&] (const std::pair<std::string, std::string> &p) { return p.first == name; });
    if (it == m_found.cend()) {
        return {false, ""};
    } else {
        const std::string value = it->second;
        m_found.erase(it);
        return {true, value};
    }
}

bool From_Switch(const std::string &name, ArgParser &a) {
    auto o = a.consume(name);
    if (o.first) {
        return true;
    } else {
        return false;
    }
}

void Help(ArgParser &a, const std::string &usage) {
    if (a.found("help")) {
        std::cout << usage << std::endl << std::endl;
        for (auto &o : a.options()) {
            std::cout << o << std::endl;
        }
        exit(EXIT_SUCCESS);
    }
}

}