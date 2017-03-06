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
#include <algorithm>

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
    const std::vector<TOption> &m_opts;
    std::deque<std::pair<std::string, std::string>> m_args;
    std::deque< std::string> m_nonopts;

public:
    ArgParser(int argc, char **argv, const std::string &usage,
              const std::vector<TOption> &opts);

    const std::vector<TOption> &options() { return m_opts; }
    const std::deque<std::string> &nonoptions() { return m_nonopts; }
    bool option_present(const std::string &name);
    std::pair<bool, std::string> consume(const std::string &name);

    template<typename T = std::string>
    T option_value(const std::string &name, const T &def_value = "") {
        auto it = std::find_if(m_args.cbegin(), m_args.cend(), [&] (const std::pair<std::string, std::string> &p) { return p.first == name; });
        if (it == m_args.cend()) {
            return def_value;
        } else {
            std::stringstream ss(it->second);
            T temp;
            ss >> temp;
            m_args.erase(it);
            return temp;
        }
    }
};

} // End namespace QI

#endif // QI_ARGS_H