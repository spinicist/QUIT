#include <functional>
#include <map>
#include <string>

int glm_contrasts_main(int argc, char **argv);
int glm_setup_main(int argc, char **argv);
int rois_main(int argc, char **argv);

void add_stats_commands(std::map<std::string, std::function<int(int, char **)>> &commands);
