/*
 *  Args.h
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_ARGS_H
#define QI_ARGS_H

#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "QI/Macro.h"

namespace QI {

struct TOption {
    const std::string long_name;
    const char short_name;
    const std::string usage;
    const bool has_arg;
};
std::ostream &operator<< (std::ostream &os, const TOption &o);

class ArgParser {
protected:
    std::vector<const TOption> m_opts;
    std::deque<std::pair<std::string, std::string>> m_found;
    std::deque<const std::string> m_nonopts;

public:
    ArgParser(int argc, char **argv, std::vector<const TOption> opts);

    const std::vector<const TOption> &options() { return m_opts; }
    std::deque<const std::string> &nonoptions() { return m_nonopts; }
    bool found(const std::string &name);
    std::pair<bool, const std::string> consume(const std::string &name);
};

template<typename T>
T From_Option(const std::string &name, const T &def_value, ArgParser &a) {
    auto o = a.consume(name);
    if (o.first) {
        std::stringstream ss(o.second);
        T temp;
        ss >> temp;
        return temp;
    } else {
        return def_value;
    }
}
bool From_Switch(const std::string &name, ArgParser &a);

void Help(ArgParser &a, const std::string &usage);

} // End namespace QI

#endif // QI_ARGS_H