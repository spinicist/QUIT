#include <functional>
#include <map>
#include <string>

int diff_main(int argc, char **argv);
int hdr_main(int argc, char **argv);
int newimage_main(int argc, char **argv);

void add_core_commands(std::map<std::string, std::function<int(int, char **)>> &commands);
