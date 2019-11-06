#include <functional>
#include <map>
#include <string>

int ase_oef_main(int argc, char **argv);
int asl_main(int argc, char **argv);
int zshim_main(int argc, char **argv);

void add_perfusion_commands(std::map<std::string, std::function<int(int, char **)>> &commands);
