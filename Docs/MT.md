# MT

MR voxels often contain complex microstructure with multiple different components or pools, each with unique relaxation properties. It is possible for magnetization to be transferred between these pools via several mechanisms, such as exchange of individual protons or entire molecules, or simple dipolar coupling from molecules that are in close proximity. These mechanisms can be studied in the related fields of Magnetization Transfer (MT) and Chemical Exchange Saturation Transfer (CEST). QUIT currently contains some basic CEST analysis tools and one for calculating simple dipolar/inhomogeneous MT ratios.

* [qi_lorentzian](#qi_lorentzian)
* [qi_mtasym](#qi_mtasym)
* [qi_dipolar_mtr](#qi_dipolar_mtr)

## qi_lorentzian

Fits a single Lorentzian to a Z-spectrum for B0 correction. Currently hard-coded to only fit the spectrum between +/-2ppm to avoid background MT contamination.

**Example Command Line**

```bash
qi_lorentzian zspectrum.nii.gz < input.txt
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
