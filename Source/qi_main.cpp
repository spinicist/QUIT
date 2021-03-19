#include "Args.h"
#include "Commands.h"
#include "Util.h"
#include <iostream>

int main(int argc, char **argv) {
    args::ArgumentParser parser("http://github.com/spinicist/QUIT");
    args::GlobalOptions  globals(parser, global_group);

#define ADD(CMD, GROUP, HELP) args::Command CMD(GROUP, #CMD, HELP, &CMD##_main);

    args::Group core(parser, "CORE");
    args::Flag  version(core, "VERSION", "Print the version of QUIT", {"version"});
    ADD(newimage, core, "Create a new image");
    ADD(diff, core, "Calcualte the difference between two images");
    ADD(hdr, core, "Print header information from an image");
#ifdef BUILD_B1
    args::Group b1(parser, "B1");
    ADD(afi, b1, "Actual Flip-Angle Imaging");
    ADD(b1_papp, b1, "B1- estimation via Papp method");
    ADD(dream, b1, "DREAM B1-Mapping");
#endif
#ifdef BUILD_MT
    args::Group mt(parser, "MT");
    ADD(lineshape, mt, "Calculate lineshapes");
    ADD(lorentzian, mt, "Lorentzian fitting for CEST");
    ADD(mtr, mt, "Calculate Magnetization Transfer Ratios");
    ADD(mtsat, mt, "Calculate MTsat");
    ADD(qmt, mt, "Fit qMT data with Ramani model");
    ADD(ssfp_emt, mt, "Fit qMT SSFP data with Wood model");
    ADD(zspec_b1, mt, "B1 correction for Z-spectra");
    ADD(zspec_interp, mt, "Interpolate Z-spectra (B0 correction etc.)");
#endif
#ifdef BUILD_PERFUSION
    args::Group perfusion(parser, "PERFUSION");
    ADD(ase_oef, perfusion, "Fit the Blockley model for Oxygen Extraction Fraction");
    ADD(asl, perfusion, "Calculate perfusion values with standard ASL model");
    ADD(zshim, perfusion, "Perform Z-shimming");
#endif
#ifdef BUILD_RELAX
    args::Group relax(parser, "RELAXOMETRY");
    ADD(despot1, relax, "DESPOT1/VFA methods");
    ADD(despot1hifi, relax, "DESPOT1-HIFI simultaneous T1/B1 mapping");
    ADD(despot2, relax, "DESPOT2");
    ADD(despot2fm, relax, "DESPOT2-FM simultaneous T2/B0 mapping");
    ADD(jsr, relax, "Fit Joint System Relaxometry model (T1/T2)")
    ADD(mcdespot, relax, "mcDESPOT");
    ADD(mpm_r2s, relax, "Fit for R2* across multiple contrasts (ECSTATICS)");
    ADD(mp2rage, relax, "MP2-RAGE estimation of T1");
    ADD(multiecho, relax, "Multi-echo T2/T2*");
    ADD(planet, relax, "PLANET method of T1/T2 mapping");
    ADD(ssfp_ellipse, relax, "Calculuate the SSFP ellipse properties");
#endif
#ifdef BUILD_PARMESAN
    args::Group parmesan(parser, "PARMESAN");
    ADD(transient, parmesan, "Transient type PARMESAN");
    ADD(ss, parmesan, "Steady-state PARMESAN");
    ADD(rf_sim, parmesan, "Calculate RF pulse properties for PARMESAN");
#endif
#ifdef BUILD_STATS
    args::Group stats(parser, "STATS");
    ADD(glm_contrasts, stats, "Calculate group means etc. for a GLM");
    ADD(glm_setup, stats, "Setup files for input to FSL randomise etc.");
    ADD(rois, stats, "Extract ROI values");
#endif
#ifdef BUILD_SUSCEP
    args::Group suscep(parser, "SUSCEPTIBILITY");
    ADD(fieldmap, suscep, "Calculate a B0 map from multi-echo data");
    ADD(unwrap_laplace, suscep, "Laplacian phase unwrapping");
    ADD(unwrap_path, suscep, "Path-based phase unwrapping");
#endif
#ifdef BUILD_UTILS
    args::Group utils(parser, "UTILITIES");
    ADD(affine, utils, "Edit the affine transformation in an image header");
    ADD(affine_angle, utils, "Calculuate the angle from the Z-axis of the header affine transform");
    ADD(coil_combine, utils, "Combine images from multi-channel coils");
    ADD(complex, utils, "Convert real/imaginary/magnitude/phase/complex data");
    ADD(denoise, utils, "Apply TGV denoising");
    ADD(gradient, utils, "Calculate the gradients of an image");
    ADD(kfilter, utils, "Filter an image via k-space");
    ADD(mask, utils, "Calculate a mask using various threshold based methods");
    ADD(noise_est, utils, "Calculate noise values in a region/mask");
    ADD(pca, utils, "Perform PCA noise reduction on multi-volume data");
    ADD(polyfit, utils, "Fit a polynomial to an image");
    ADD(polyimg, utils, "Create an image from a polynomial");
    ADD(rfprofile, utils, "Multiply a B1 map by a slab profile");
    ADD(select, utils, "Choose volumes from a 4D image");
    ADD(ssfp_bands, utils, "Remove banding artefacts from SSFP images");
    ADD(tvmask, utils, "Calculate a mask from 4D data using Total-Variation");
#endif
#undef ADD

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cerr << parser << '\n';
        exit(EXIT_SUCCESS);
    } catch (args::Error e) {
        if (version) {
            std::cout << QI::GetVersion() << '\n';
            exit(EXIT_SUCCESS);
        }
        std::cerr << parser << '\n' << e.what() << '\n';
        exit(EXIT_FAILURE);
    }

    exit(EXIT_SUCCESS);
}