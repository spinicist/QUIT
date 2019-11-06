#include <functional>
#include <map>
#include <string>

int lineshape_main(int argc, char **argv);
int lorentzian_main(int argc, char **argv);
int mtsat_main(int argc, char **argv);
int qmt_main(int argc, char **argv);
int ssfp_emt_main(int argc, char **argv);
int zspec_b1_main(int argc, char **argv);
int zspec_interp_main(int argc, char **argv);

void add_mt_commands(std::map<std::string, std::function<int(int, char **)>> &commands);
