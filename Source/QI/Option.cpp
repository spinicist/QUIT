/*
 *  Option.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Option.h"

namespace QI {

// Global default options
OptionList &DefaultOptions() {
    static OptionList s_options;
    return s_options;
}

// OptionBase
OptionBase::OptionBase(const char s, const std::string &l, const std::string &u,
                       OptionList &list) :
    m_short(s), m_long(l), m_usage(u)
{
    list.add(this);
}

std::ostream &operator<< (std::ostream &os, const OptionBase &o) {
    os << "  -" << o.shortOption() << ",--"
       << std::setw(10) << std::left << o.longOption()
       << std::setw(0) << " : " << o.usage();
    return os;
}

// OptionList
void OptionList::add(OptionBase *o) {
    m_options.push_front(o);
}

void OptionList::parse(int argc, char *const *argv, std::vector<std::string> &nonopts) {
    int optind = 1;
    while (optind < argc) {
        std::string thisopt(argv[optind]);
        auto it = m_options.end();
        std::string arg;
        if (thisopt == "--") { // Premature end of options
            optind++;
            while (optind < argc) {
                thisopt = argv[optind++];
                nonopts.push_back(thisopt);
            }
            return;
        } else if (thisopt[0] != '-') { // Non-option
            nonopts.push_back(thisopt);
        } else if (thisopt[1] == '-') { // Long option
            size_t epos = thisopt.find("=");
            std::string arg;
            if (epos != std::string::npos) {
                arg = thisopt.substr(epos + 1);
                thisopt = thisopt.substr(2, epos - 2);
            } else {
                thisopt.erase(0, 2);
            }
            it = std::find_if(m_options.begin(), m_options.end(), [&] (OptionBase *const &o) { return o->longOption() == thisopt; });
            if (it == m_options.end()) {
                QI_EXCEPTION("Unhandled long option: " + thisopt);
            }
            if ((*it)->hasArgument()) {
                (*it)->setArgument(arg);
            } else {
                (*it)->setValue();
            }
        } else { // Short option
            char sopt = thisopt[1];
            it = std::find_if(m_options.begin(), m_options.end(), [&] (OptionBase *const &o) { return o->shortOption() == sopt; });
            if (it == m_options.end()) {
                QI_EXCEPTION("Unhandled short option: " + std::string(thisopt));
            }
            if ((*it)->hasArgument()) {
                if (thisopt.size() > 2) { // Handle arguments without spaces
                    (*it)->setArgument(thisopt.substr(2));
                } else if ((optind+1) == argc) {
                    QI_EXCEPTION("Missing required argument for option " + std::to_string(sopt));
                } else {
                    (*it)->setArgument(argv[++optind]);
                }
            } else {
                (*it)->setValue();
            }
        }
        optind++;
    }
}

void OptionList::print(std::ostream &os) const {
    for (OptionBase *o : m_options) {
        os << *o << std::endl;
    }
}

std::ostream &operator<< (std::ostream &os, const OptionList &o) {
    o.print(os);
    return os;
}

} // End namespace QI