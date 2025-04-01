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
#include "fmt/ranges.h"
using namespace fmt::literals;

namespace QI {

enum struct Level {
    Info, Warn, Fail
};

void SetVerbose(bool const vb);
void FormatLog(Level lvl, fmt::string_view fmt, fmt::format_args args);
__attribute__((noreturn)) void FormatFail(fmt::string_view fmt, fmt::format_args args);

template <typename... Args>
inline void Log(Level lvl, fmt::format_string<Args...> fstr, Args &&...args) {
        FormatLog(lvl, fstr, fmt::make_format_args(args...));
}

template <typename... Args>
inline void Info(fmt::format_string<Args...> fstr, Args &&...args) {
        FormatLog(Level::Info, fstr, fmt::make_format_args(args...));
}

template <typename... Args>
inline void Warn(fmt::format_string<Args...> fstr, Args &&...args) {
        FormatLog(Level::Warn, fstr, fmt::make_format_args(args...));
}

template <typename... Args> __attribute__((noreturn))
inline void Fail(fmt::format_string<Args...> fstr, Args &&...args) {
        FormatFail(fstr, fmt::make_format_args(args...));
}

} // End namespace QI
