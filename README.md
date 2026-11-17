# Hubble_constant_estimation
# Estimating the Hubble Constant from Galaxy Spectra

## Project Overview

This project investigates the relationship between the distances of galaxies and their recession velocities to estimate the Hubble constant, which describes the expansion rate of the Universe.

Using observational spectral data, I extracted the central wavelengths of Hβ spectral lines through Gaussian curve fitting, calculated galaxy recession velocities using the relativistic Doppler relation, and applied linear regression to estimate the Hubble constant.

## Methods

### 1. Spectral Line Fitting

- Processed galaxy spectral data containing frequency and intensity measurements.
- Modelled each spectral feature using a Gaussian profile superimposed on a linear background.
- Used `scipy.optimize.curve_fit` to determine the best-fitting parameters, including the central frequency of each spectral line.
- Used the covariance matrix returned by the fitting procedure to estimate uncertainty in the fitted frequency.

The fitting model was:

\[
I(f)=a\exp\left[-\frac{(f-\mu)^2}{2\sigma^2}\right]+mf+c
\]

where the Gaussian component represents the spectral feature and the linear component accounts for the local background.

### 2. Recession Velocity Calculation

- Converted fitted frequencies into observed wavelengths.
- Used the rest wavelength of the Hβ line (486.1 nm) to calculate the wavelength shift.
- Applied the relativistic Doppler relation to estimate the recession velocity of each galaxy.
- Propagated the fitted frequency uncertainty into the velocity uncertainty.

### 3. Hubble Constant Estimation

- Matched the calculated recession velocities with their corresponding galaxy distance measurements.
- Excluded observations flagged with an instrument response of `E`.
- Applied weighted linear regression to the distance–velocity relationship.
- Estimated the Hubble constant from the gradient of the fitted line.

## Results

The analysis produced an approximately linear relationship between galaxy distance and recession velocity, consistent with Hubble's law.

**Estimated Hubble constant:** 

[Insert verified result and uncertainty after finalising the weighted regression.]

## Visualisations

### Spectral Line Fitting

Representative examples of Gaussian-plus-linear fitting applied to the observed spectral data.

### Hubble Diagram

Recession velocity plotted against galaxy distance, with the fitted linear relationship.

## Technologies

- Python
- NumPy
- Pandas
- SciPy (`curve_fit`, `polyfit`)
- Matplotlib

## Data and Reproducibility

The analysis was carried out using the spectral and distance datasets supplied for the project coursework. The original datasets are not included in this repository.

To run the analysis, install the required Python libraries:

```bash
pip install numpy pandas scipy matplotlib
```

The analysis notebook and figures are provided in this repository.
