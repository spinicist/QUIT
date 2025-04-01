#include "Log.h"

namespace QI {

namespace {
bool verbose = false;
}

void SetVerbose(bool const vb) {
    verbose = vb;
}

void FormatLog(Level lvl, fmt::string_view fstr, fmt::format_args args) {
    fmt::terminal_color color;
    switch (lvl) {
    default:
    case Level::Info:
        color = fmt::terminal_color::bright_green;
        break;
    case Level::Warn:
        color = fmt::terminal_color::bright_yellow;
        break;
    case Level::Fail:
        color = fmt::terminal_color::bright_red;
        break;
    }

    if (verbose || lvl != Level::Info) {
        fmt::print(stderr, fmt::fg(color), "{:%H:%M:%S} ", fmt::localtime(std::time(nullptr)));
        fmt::vprintln(stderr, fstr, args);
    }

    if (lvl == Level::Fail) {
        std::abort();
    }
}

void FormatFail(fmt::string_view fstr, fmt::format_args args) {
    fmt::print(stderr,
               fmt::fg(fmt::terminal_color::bright_red),
               "{:%H:%M:%S} ",
               fmt::localtime(std::time(nullptr)));
    fmt::vprintln(stderr, fstr, args);
    std::abort();
}

} // namespace QI