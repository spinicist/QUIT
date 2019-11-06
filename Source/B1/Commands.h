#include <functional>
#include <map>
#include <string>

int afi_main(int argc, char **argv);
int dream_main(int argc, char **argv);
int b1_papp_main(int argc, char **argv);

void add_b1_commands(std::map<std::string, std::function<int(int, char **)>> &commands);
