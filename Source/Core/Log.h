/*
 *  Log.h
 *  Part of the QUantitative Image Tools
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_LOG_H
#define QUIT_LOG_H

#include "fmt/color.h"
#include "fmt/ostream.h"
#include "fmt/time.h"
using namespace fmt::literals;

namespace QI {

template <typename S, typename... Args>
inline void Log(const bool verbose, const S &fmt_str, const Args &... args) {
    if (verbose) {
        fmt::print(stderr, fmt::fg(fmt::terminal_color::bright_white), fmt_str, args...);
        fmt::print(stderr, "\n");
    }
}

template <typename S, typename... Args>
inline void Info(const bool verbose, const S &fmt_str, const Args &... args) {
    if (verbose) {
        const std::time_t now = std::time(nullptr);
        fmt::print(stderr, fmt::fg(fmt::color::green), "{:%T} ", *std::localtime(&now));
        fmt::print(stderr, fmt::fg(fmt::terminal_color::bright_white), fmt_str, args...);
        fmt::print(stderr, "\n");
    }
}

template <typename S, typename... Args> inline void Warn(const S &fmt_str, const Args &... args) {
    const std::time_t now = std::time(nullptr);
    fmt::print(stderr, fmt::fg(fmt::color::orange), "{:%T} ", *std::localtime(&now));
    fmt::print(stderr, fmt::fg(fmt::terminal_color::bright_white), fmt_str, args...);
    fmt::print(stderr, "\n");
}

template <typename S, typename... Args>
[[noreturn]] inline void Fail(const S &fmt_str, const Args &... args) {
    fmt::print(stderr, fmt::fg(fmt::color::red), "Error");
    fmt::print(stderr, fmt_str, args...);
    exit(EXIT_FAILURE);
}

} // End namespace QI

#endif
