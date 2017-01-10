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

#include <string>
#include <list>
#include <map>

#include "QI/Macro.h"

namespace QI {

struct Option {
    const char short_name;
    const std::string long_name;
    const std::string usage;
    const bool hasValue;
};

class Arguments {
protected:
    std::list<const std::string> m_args;
    std::map<std::string, Option> m_found_options;

public:
    Arguments(int argc, char **argv) :
        m_args(argv + 1, argv + argc)
    {
        std::cout << "Raw args: " << std::endl;
        for (auto &a: m_args) {
            std::cout << a << std::endl;
        }
    }

    std::pair<bool, std::string> findOption(const Option &o) {
        auto iArg = m_args.cbegin();
        while (iArg != m_args.cend()) {
            if (*iArg == "--") { // Premature end of options, didn't find it
                return {false, ""};
            } else if ((*iArg)[0] != '-') {
                // Non-option, ignore
            } else if ((*iArg)[1] == '-') { // Long option
                size_t posEnd = iArg->find("=");
                std::string name;
                if (posEnd != std::string::npos) {
                    name = iArg->substr(2, posEnd - 2);
                } else {
                    name = iArg->substr(2, -1);
                }
                if (name == o.long_name) { // We have a match
                    std::string value = "";
                    if (posEnd != std::string::npos) {
                        value = iArg->substr(posEnd + 1);
                        if (!o.hasValue) {
                            QI_EXCEPTION("Long option: " << o.long_name << " given unexpected value: " << value);
                        }
                    } else {
                        if (o.hasValue) {
                            QI_EXCEPTION("Long option: " << o.long_name << " requires a value.");
                        }
                    }
                    m_args.erase(iArg);
                    return {true, value};
                }
            } else if (o.short_name == (*iArg)[1]) { // Matching short option
                std::string value;
                if (o.hasValue) {
                    if (iArg->size() > 2) { // Handle arguments without spaces
                        value = iArg->substr(2);
                        m_args.erase(iArg);
                    } else {
                        auto iVal = iArg;
                        ++iVal;
                        if (iVal == m_args.cend()) {
                            QI_EXCEPTION("Missing required value for option: " + std::to_string(o.short_name));
                        }
                        value = *(++iVal);
                        m_args.erase(iArg,++iVal);
                    }
                }
                return {true,value};
            }
            iArg++;
        }
        return {false,""};
    }

    std::list<const std::string> finishOptions() {
        if (m_args.front() == "--") {
            m_args.pop_front();
        } else { // Should have no more options now
            for (auto iArg = m_args.cbegin(); iArg != m_args.cend(); iArg++) {
                if ((*iArg)[0] == '-') {
                    QI_EXCEPTION("Unhandled option: " << *iArg);
                }
            }
        }
        return m_args;
    }
};

template<typename T>
T From_Option(const T &def_value, const Option &opt, Arguments &a);

} // End namespace QI

#endif // QI_ARGS_H