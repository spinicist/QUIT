#include <functional>
#include <map>
#include <string>

int rufis_mupa_main(int argc, char **argv);
int rf_sim_main(int argc, char **argv);
int rufis_ss_main(int argc, char **argv);

void add_rufis_commands(std::map<std::string, std::function<int(int, char **)>> &commands);
