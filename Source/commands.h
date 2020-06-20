#pragma once

int diff_main(args::Subparser &parser);
int hdr_main(args::Subparser &parser);
int newimage_main(args::Subparser &parser);

#ifdef BUILD_B1
int afi_main(args::Subparser &parser);
int dream_main(args::Subparser &parser);
int b1_papp_main(args::Subparser &parser);
#endif
#ifdef BUILD_MT
int lineshape_main(args::Subparser &parser);
int lorentzian_main(args::Subparser &parser);
int mtr_main(args::Subparser &parser);
int mtsat_main(args::Subparser &parser);
int qmt_main(args::Subparser &parser);
int ssfp_emt_main(args::Subparser &parser);
int zspec_b1_main(args::Subparser &parser);
int zspec_interp_main(args::Subparser &parser);
#endif
#ifdef BUILD_PERFUSION
int ase_oef_main(args::Subparser &parser);
int asl_main(args::Subparser &parser);
int zshim_main(args::Subparser &parser);
#endif
#ifdef BUILD_RELAX
int jsr_main(args::Subparser &parser);
int mpm_r2s_main(args::Subparser &parser);
int planet_main(args::Subparser &parser);
int ssfp_ellipse_main(args::Subparser &parser);
int despot1_main(args::Subparser &parser);
int despot1hifi_main(args::Subparser &parser);
int despot2_main(args::Subparser &parser);
int despot2fm_main(args::Subparser &parser);
int mcdespot_main(args::Subparser &parser);
int mp2rage_main(args::Subparser &parser);
int multiecho_main(args::Subparser &parser);
int planet_main(args::Subparser &parser);
#endif
#ifdef BUILD_PARMESAN
int transient_main(args::Subparser &parser);
int rf_sim_main(args::Subparser &parser);
int ss_main(args::Subparser &parser);
#endif
#ifdef BUILD_STATS
int glm_contrasts_main(args::Subparser &parser);
int glm_setup_main(args::Subparser &parser);
int rois_main(args::Subparser &parser);
#endif
#ifdef BUILD_SUSCEP
int fieldmap_main(args::Subparser &parser);
int unwrap_laplace_main(args::Subparser &parser);
int unwrap_path_main(args::Subparser &parser);
#endif
#ifdef BUILD_UTILS
int affine_main(args::Subparser &parser);
int affine_angle_main(args::Subparser &parser);
int coil_combine_main(args::Subparser &parser);
int gradient_main(args::Subparser &parser);
int pca_main(args::Subparser &parser);
int rfprofile_main(args::Subparser &parser);
int select_main(args::Subparser &parser);
int ssfp_bands_main(args::Subparser &parser);
int tvmask_main(args::Subparser &parser);
int complex_main(args::Subparser &parser);
int kfilter_main(args::Subparser &parser);
int mask_main(args::Subparser &parser);
int polyfit_main(args::Subparser &parser);
int polyimg_main(args::Subparser &parser);
#endif