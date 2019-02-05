/*
 *  Args.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Args.h"
#include "Util.h"

/*
 * Moved these here to avoid duplicate symbol errors
 */
namespace args {
/** (INTERNAL) Count UTF-8 glyphs
 *
 * This is not reliable, and will fail for combinatory glyphs, but it's
 * good enough here for now.
 *
 * \param string The string to count glyphs from
 * \return The UTF-8 glyphs in the string
 */
std::string::size_type Glyphs(const std::string &string_) {
    std::string::size_type length = 0;
    for (const char c : string_) {
        if ((c & 0xc0) != 0x80) {
            ++length;
        }
    }
    return length;
}

/** (INTERNAL) Wrap a string into a vector of lines
 *
 * This is quick and hacky, but works well enough.  You can specify a
 * different width for the first line
 *
 * \param width The width of the body
 * \param the widtho f the first line, defaults to the width of the body
 * \return the vector of lines
 */
std::vector<std::string> Wrap(const std::string &in, const std::string::size_type width,
                              std::string::size_type firstlinewidth) {
    // Preserve existing line breaks
    const auto newlineloc = in.find('\n');
    if (newlineloc != in.npos) {
        auto first  = Wrap(std::string(in, 0, newlineloc), width);
        auto second = Wrap(std::string(in, newlineloc + 1), width);
        first.insert(std::end(first), std::make_move_iterator(std::begin(second)),
                     std::make_move_iterator(std::end(second)));
        return first;
    }
    if (firstlinewidth == 0) {
        firstlinewidth = width;
    }
    auto currentwidth = firstlinewidth;

    std::istringstream       stream(in);
    std::vector<std::string> output;
    std::ostringstream       line;
    std::string::size_type   linesize = 0;
    while (stream) {
        std::string item;
        stream >> item;
        auto itemsize = Glyphs(item);
        if ((linesize + 1 + itemsize) > currentwidth) {
            if (linesize > 0) {
                output.push_back(line.str());
                line.str(std::string());
                linesize     = 0;
                currentwidth = width;
            }
        }
        if (itemsize > 0) {
            if (linesize) {
                ++linesize;
                line << " ";
            }
            line << item;
            linesize += itemsize;
        }
    }
    if (linesize > 0) {
        output.push_back(line.str());
    }
    return output;
}
std::ostream &operator<<(std::ostream &os, const ArgumentParser &parser) {
    parser.Help(os);
    return os;
}

} // namespace args

namespace QI {

void ParseArgs(args::ArgumentParser &parser, int argc, char **argv, const args::Flag &verbose) {
    try {
        parser.ParseCLI(argc, argv);
        QI::Log(verbose, "Starting {} {}", argv[0], QI::GetVersion());
    } catch (args::Help) {
        fmt::print("{}\n", parser);
        exit(EXIT_SUCCESS);
    } catch (args::ParseError e) {
        QI::Fail("{}\n{}", parser, e.what());
    } catch (args::ValidationError e) {
        QI::Fail("{}\n{}", parser, e.what());
    }
}

void ParseArgs(args::ArgumentParser &parser, int argc, char **argv, const args::Flag &verbose,
               args::ValueFlag<int> &threads) {
    try {
        parser.ParseCLI(argc, argv);
        QI::Log(verbose, "Starting {} {}", argv[0], QI::GetVersion());
        QI::Log(verbose, "Max threads = {}", threads.Get());
        itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(threads.Get());
    } catch (args::Help) {
        fmt::print("{}\n", parser);
        exit(EXIT_SUCCESS);
    } catch (args::ParseError e) {
        QI::Fail("{}\n{}", parser, e.what());
    } catch (args::ValidationError e) {
        QI::Fail("{}\n{}", parser, e.what());
    }
}

} // End namespace QI
