#include <functional>
#include <map>
#include <string>

int coil_combine_main(int argc, char **argv);
int gradient_main(int argc, char **argv);
int pca_main(int argc, char **argv);
int rfprofile_main(int argc, char **argv);
int select_main(int argc, char **argv);
int ssfp_bands_main(int argc, char **argv);
int tvmask_main(int argc, char **argv);
int affine_main(int argc, char **argv);
int complex_main(int argc, char **argv);
int kfilter_main(int argc, char **argv);
int mask_main(int argc, char **argv);
int polyfit_main(int argc, char **argv);
int polyimg_main(int argc, char **argv);

void add_utils(std::map<std::string, std::function<int(int, char **)>> &commands);
