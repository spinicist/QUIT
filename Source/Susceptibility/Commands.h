#include <functional>
#include <map>
#include <string>

int fieldmap_main(int argc, char **argv);
int unwrap_laplace_main(int argc, char **argv);
int unwrap_path_main(int argc, char **argv);

void add_suscep_commands(std::map<std::string, std::function<int(int, char **)>> &commands);
