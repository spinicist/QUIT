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

#pragma once

#define FMT_DEPRECATED_OSTREAM

#include "fmt/chrono.h"
#include "fmt/color.h"
#include "fmt/ostream.h"
using namespace fmt::literals;

namespace QI {

template <typename... Args>
inline void Log(const bool verbose, fmt::format_string<Args...> fmt_str, Args &&...args) {
    if (verbose) {
        fmt::print(stderr, fmt_str, std::forward<Args>(args)...);
        fmt::print(stderr, "\n");
    }
}

template <typename... Args>
inline void Info(const bool verbose, fmt::format_string<Args...> fmt_str, Args &&...args) {
    if (verbose) {
        const std::time_t now = std::time(nullptr);
        fmt::print(
            stderr, fmt::fg(fmt::terminal_color::bright_green), "{:%T} ", *std::localtime(&now));
        fmt::print(stderr, fmt_str, std::forward<Args>(args)...);
        fmt::print(stderr, "\n");
    }
}

template <typename... Args>
inline void Warn(fmt::format_string<Args...> fmt_str, Args &...args) {
    const std::time_t now = std::time(nullptr);
    fmt::print(
        stderr, fmt::fg(fmt::terminal_color::bright_yellow), "{:%T} ", *std::localtime(&now));
    fmt::print(stderr, fmt_str, std::forward<Args>(args)...);
    fmt::print(stderr, "\n");
}

template <typename... Args>
__attribute__((noreturn)) inline void Fail(fmt::format_string<Args...> fmt_str, Args &&...args) {
    fmt::print(stderr, fmt::fg(fmt::terminal_color::bright_red), "Error ");
    fmt::print(stderr, fmt_str, std::forward<Args>(args)...);
    fmt::print(stderr, "\n");
    exit(EXIT_FAILURE);
}

} // End namespace QI
