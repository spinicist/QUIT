#include "Log.h"

#include "B1/Commands.h"
#include "CoreProgs/Commands.h"
#include "MT/Commands.h"
#include "Perfusion/Commands.h"
#include "RUFIS/Commands.h"
#include "Relaxometry/Commands.h"
#include "Stats/Commands.h"
#include "Susceptibility/Commands.h"
#include "Utils/Commands.h"

#include <functional>

using MainFunc = std::function<int(int, char **)>;

int main(int argc, char **argv) {

    std::map<std::string, MainFunc> commands;

    add_core_commands(commands);
#ifdef BUILD_B1
    add_b1_commands(commands);
#endif
#ifdef BUILD_MT
    add_mt_commands(commands);
#endif
#ifdef BUILD_PERFUSION
    add_perfusion_commands(commands);
#endif
#ifdef BUILD_RELAX
    add_relax_commands(commands);
#endif
#ifdef BUILD_RUFIS
    add_rufis_commands(commands);
#endif
#ifdef BUILD_STATS
    add_stats_commands(commands);
#endif
#ifdef BUILD_SUSCEP
    add_suscep_commands(commands);
#endif
#ifdef BUILD_UTILS
    add_utils(commands);
#endif

    auto print_command_list = [&commands]() {
        fmt::print("Available commands:\n");
        size_t max_width = 0;
        for (auto const &kv_pair : commands) {
            auto const &this_width = kv_pair.first.size();
            if (this_width > max_width) {
                max_width = this_width;
            }
        }
        auto num_printed = 0;
        fmt::print("  ");
        for (auto const &kv_pair : commands) {
            fmt::print("{:<{}}", kv_pair.first, max_width + 1);
            num_printed++;
            if (num_printed % 5 == 0) {
                fmt::print("\n  ");
            }
        }
        fmt::print("\n");
    };

    if (argc < 2) {
        print_command_list();
        QI::Fail("Must specify command");
    }

    std::string command(argv[1]);
    auto        find_command = commands.find(command);
    if (find_command == commands.end()) {
        print_command_list();
        QI::Fail("Unknown command {}", command);
    }
    return find_command->second(argc - 1, argv + 1);
}