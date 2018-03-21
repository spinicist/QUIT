# SSFP

The Steady-State Free-Precession (SSFP), or more precisely balanced-SSFP (bSSFP), sequence is one of the oldest NMR sequences and can be used to give high SNR MR images with mixed T1/T2 contrast in very short scan time. However, it suffers from banding artefacts in areas of off-resonance which limit its clinical applicability. This module contains a tool for removing those banding artefacts, and then further tools for quantitative mapping using the ellipse signal model.

* [qissfpbands](#qissfpbands)
* [qi_ssfp_ellipse](#qi_ssfp_ellipse)
* [qi_ssfp_relax](#qi_ssfp_relax)
* [qi_ssfp_emt](#qi_ssfp_emt)

## qissfpbands

There are several different methods for removing SSFP bands in the literature. Most of them rely on acquiring multiple SSFP images with different phase-increments (also called phase-cycling or phase-cycling patterns). Changing the phase-increments moves the bands to a different location, after which the images can be combined to reduce the banding. The different approaches are discussed further below.

**Example Command Line**

```bash
qissfpbands
```

The Z-spectrum must be a 4D file with each volume acquired at a different offset frequency.

**Example Input File**

```json
{
    "freq" : [ -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
}
```

These are the offset frequencies for each volume in the Z-spectrum input.

**Outputs**

* `LTZ_f0.nii.gz`  - The center frequency of the fitted Lorentzian.
* `LTZ_w.nii.gz`   - The width of the fitted Lorentzian.
* `LTZ_sat.nii.gz` - The saturation ratio of the fitted Lorentzian.
* `LTZ_PD.nii.gz`  - The apparent Proton Density of the fitted Lorentzian.

## qi_mtasym

Calculates the MT asymmetry of a Z-spectrum.

**Example Command Line**

```bash
qi_mtasym zspectrum.nii.gz --f0=LTZ_f0.nii.gz < input.txt
```

The off-resonance map units must match the input frequencies (e.g. either PPM or Hertz)

**Example Input File**

```json
{
    "freq" : [ -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5],
    "asym_freq" : [-3, -2]
}
```

`freq` is the offset frequencies the Z-spectrum was acquired at. `asym_freq` are the frequencies you want the asymmetry calculated at.

**Outputs**

* `MT_asymmetry.nii.gz` The asymmetry value at each asymmetry frequency.

**Important Options**

* `--f0, -f`

    Specify an off-resonance map. Units must be the same as the input & asymmetry frequencies.

##qi_dipolar_mtr

Calculates dipolar/inhomogeneous Magnetization Transfer Ratios (MTRs). Dipolar/inhomogeneous MT is a new (see note) contrast mechanism that is present in highly structured materials such as myelin and tendon. By applying off-resonance saturation at both positive and negative frequencies (instead of only one side as in classic MTR) it is possible to decouple the dipolar pool and hence produce an enhanced Magnetization Transfer (eMT) effect. The different between eMT and normal MT is the dipolar/inhomogeneous MT and is potentially highly specific to myelin within the brain.

Although the majority of the existing literature refers to this effect as inhomogeneous MT, this name was chosen before the physical phenomena underlying the effect was well understood. Current theory does not rely on inhomogeneous effects at all, so the name is a misnomer.

*Note - The original ihMT abstracts are from around 2005. There was 10 years between the conference abstracts and the corresponding full papers. So the method is not that new*

**Example Command Line**

```bash
qi_dipolar_mtr dipolar_mt_volumes.nii.gz
```

The input must consist of 5 volumes: Dipolar +/-, Dipolar -/+, Unsaturated, MT+, MT-. This scheme is not flexible and will be improved in a future version.

**Outputs**

* `DMT_mtr.nii.gz` - The classic MTR, expressed as a percentage
* `DMT_emtr.nii.gz` - The enhanced MTR, expressed as a percentage
* `DMT_dmtr.nii.gz` - The dipolar MTR, expressed as a percentage. This is the difference between eMTR and MTR.
* `DMT_mta.nii.gz` - The first-order MT-asymmetry (MT- subtracted from MT+, relative to unsaturated, in percent).
