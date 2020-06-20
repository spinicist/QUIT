#include "Args.h"
#include "Commands.h"
#include "Util.h"
#include <iostream>

int main(int argc, char **argv) {
    args::ArgumentParser parser("QUIT http://github.com/spinicist/QUIT");
    args::Group          commands(parser, "COMMANDS");
    args::GlobalOptions  globals(parser, global_group);
    args::Flag           version(parser, "VERSION", "Print the version of QUIT", {"version"});

#define ADD(CMD, HELP) args::Command CMD(commands, #CMD, HELP, &CMD##_main);
    ADD(newimage, "Create a new image");
    ADD(diff, "Calcualte the difference between two images");
    ADD(hdr, "Print header information from an image");
#ifdef BUILD_B1
    ADD(afi, "Actual Flip-Angle Imaging");
    ADD(dream, "DREAM B1-Mapping");
    ADD(b1_papp, "B1- estimation via Papp method");
#endif
#ifdef BUILD_MT
    ADD(lineshape, "Calculate lineshapes");
    ADD(lorentzian, "Lorentzian fitting for CEST");
    ADD(mtr, "Calculate Magnetization Transfer Ratios");
    ADD(mtsat, "Calculate MTsat");
    ADD(qmt, "Fit qMT data with Ramani model");
    ADD(ssfp_emt, "Fit qMT SSFP data with Wood model");
    ADD(zspec_b1, "B1 correction for Z-spectra");
    ADD(zspec_interp, "Interpolate Z-spectra (B0 correction etc.)");
#endif
#ifdef BUILD_PERFUSION
    ADD(ase_oef, "Fit the Blockley model for Oxygen Extraction Fraction");
    ADD(asl, "Calculate perfusion values with standard ASL model");
    ADD(zshim, "Perform Z-shimming");
#endif
#ifdef BUILD_RELAX
    ADD(jsr, "Fit Joint System Relaxometry model (T1/T2)")
    ADD(mpm_r2s, "Fit for R2* across multiple contrasts (ECSTATICS)");
    ADD(ssfp_ellipse, "Calculuate the SSFP ellipse properties");
    ADD(despot1, "DESPOT1/VFA methods");
    ADD(despot1hifi, "DESPOT1-HIFI simultaneous T1/B1 mapping");
    ADD(despot2, "DESPOT2");
    ADD(despot2fm, "DESPOT2-FM simultaneous T2/B0 mapping");
    ADD(mp2rage, "MP2-RAGE estimation of T1");
    ADD(multiecho, "Multi-echo T2/T2*");
    ADD(planet, "PLANET method of T1/T2 mapping");
#endif
#ifdef BUILD_PARMESAN
    ADD(transient, "Transient type PARMESAN");
    ADD(ss, "Steady-state PARMESAN");
    ADD(rf_sim, "Calculate RF pulse properties for PARMESAN");
#endif
#ifdef BUILD_STATS
    ADD(glm_contrasts, "Calculate group means etc. for a GLM");
    ADD(glm_setup, "Setup files for input to FSL randomise etc.");
    ADD(rois, "Extract ROI values");
#endif
#ifdef BUILD_SUSCEP
    ADD(fieldmap, "Calculate a B0 map from multi-echo data");
    ADD(unwrap_laplace, "Laplacian phase unwrapping");
    ADD(unwrap_path, "Path-based phase unwrapping");
#endif
#ifdef BUILD_UTILS
    ADD(affine, "Edit the affine transformation in an image header");
    ADD(affine_angle, "Calculuate the angle from the Z-axis of the header affine transform");
    ADD(coil_combine, "Combine images from multi-channel coils");
    ADD(gradient, "Calculate the gradients of an image");
    ADD(pca, "Perform PCA noise reduction on multi-volume data");
    ADD(rfprofile, "Multiply a B1 map by a slab profile");
    ADD(select, "Choose volumes from a 4D image");
    ADD(ssfp_bands, "Remove banding artefacts from SSFP images");
    ADD(tvmask, "Calculate a mask from 4D data using Total-Variation");
    ADD(complex, "Convert real/imaginary/magnitude/phase/complex data");
    ADD(kfilter, "Filter an image via k-space");
    ADD(mask, "Calculate a mask using various threshold based methods");
    ADD(polyfit, "Fit a polynomial to an image");
    ADD(polyimg, "Create an image from a polynomial");
#endif
#undef ADD

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cerr << parser << '\n';
        exit(EXIT_SUCCESS);
    } catch (args::Error e) {
        std::cerr << parser << '\n' << e.what() << '\n';
        exit(EXIT_FAILURE);
    }
    if (version) {
        std::cout << QI::GetVersion() << '\n';
    }
    exit(EXIT_SUCCESS);
}