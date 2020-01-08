#include <functional>
#include <map>
#include <string>

int jsr_main(int argc, char **argv);
int mpm_r2s_main(int argc, char **argv);
int planet_main(int argc, char **argv);
int ssfp_ellipse_main(int argc, char **argv);
int despot1_main(int argc, char **argv);
int despot1hifi_main(int argc, char **argv);
int despot2_main(int argc, char **argv);
int despot2fm_main(int argc, char **argv);
int mcdespot_main(int argc, char **argv);
int mp2rage_main(int argc, char **argv);
int multiecho_main(int argc, char **argv);

void add_relax_commands(std::map<std::string, std::function<int(int, char **)>> &commands);
