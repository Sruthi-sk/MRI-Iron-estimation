# Project for Computational Modelling for Biomedical Imaging module

# Quantification of Hepatic Iron Overload using TV-Regularized MRI Optimization
 
This project implements a robust and non-invasive method to quantify hepatic iron overload in thalassemia patients using MRI T2* relaxometry. The approach combines Total Variation (TV)-regularized optimization with the Alternating Direction Method of Multipliers (ADMM) to estimate T2* maps from multi-echo MR images, validated on both synthetic phantoms and clinical data

## Motivation

Thalassemia patients often require repeated blood transfusions, leading to excessive iron accumulation in the liver. Accurate quantification of iron concentration is crucial for initiating and tailoring chelation therapy. While liver biopsy is the gold standard, it is invasive, expensive, and prone to sampling error. MRI-based T2* mapping provides a non-invasive alternative, and this project enhances its precision using advanced optimization techniques.

### Key Features

- T2* Estimation: Fits exponential decay curves to multi-echo signal intensities.

- TV-Regularized ADMM Solver: Enforces spatial smoothness while preserving anatomical boundaries.

- Validation on Synthetic and Clinical Data: Tested on simulated phantoms and patient MRI scans.

- Quantitative Evaluation: Achieves <5% error compared to console-based estimation methods.

## T2Star Relaxation Estimation in MRI

This repository contains MATLAB code for estimating T2* (T2-star) relaxation parameters from MRI data using total variation (TV)-regularized optimization and phantom validation. The project includes phantom creation, patient DICOM data processing, parameter estimation, and visualization of signal intensity and T2* maps.

## Contents

T2* Estimation Algorithm
- `relaxationEst.m`: Core function to estimate signal amplitude `a` and relaxation rate `r = 1/T2*` using an alternating minimization approach with TV regularization.
- `gradientDescentAK.m`: Gradient descent solver to optimize pixel-wise T2* parameters.
- `chambolle_prox_TV_stop.m`: Implementation of Chambolle’s total variation denoising algorithm.
Phantom Simulation
- `createPhantoms.m`: Generates synthetic MRI phantom data with exponential decay.varying T2* and S0 values for testing
- `main.m`: Tests the algorithm on phantoms and compares estimated vs. ground truth T2* and `a`.
Patient Data Loading & Processing
- `Patient_main.m`: Processes real DICOM MRI data from multiple patients and applies the estimator to a selected ROI.
- `analyzePatientdata.m`: Standalone tool to load a single patient’s MRI data and compute mean signal decay over ROI.
- `plotPatientdata.m`: Visualizes estimated `a` and T2* images and calculates mean T2* values in ROI.
- `Results.txt`, `ResultsGT.txt`: Benchmark T2* and `a` values from the scanner console and algorithm outputs.
Comparison to Levenberg–Marquardt
- `fitmodel_lm.m`: baseline model using the Levenberg–Marquardt (LM) algorithm for comparison

## Methodology

The estimation framework minimizes the following energy functional:

```math
min_{a, r} ||y - a * exp(-r * TE)||^2 + λ_A TV(a) + λ_R TV(r)

Where:

y: observed signals across TE (echo times)

a: signal amplitude map

r: relaxation rate map (r = 1/T2*)

TV(.): total variation for regularization

λ_A, λ_R: regularization weights

The optimization proceeds via alternating updates using:

Chambolle’s algorithm for TV proximal mapping

Pixel-wise gradient descent for r
```

### Sample Results
Estimated a and T2* maps are produced for each patient and phantom.

Mean T2* values are computed over a circular region of interest and compared to scanner-reported values.

### Requirements
- MATLAB 
- Image Processing Toolbox
- Access to DICOM-formatted MRI images

### Dataset
The scripts expect DICOM image folders for each patient
Each folder should contain multiple .IMA files corresponding to different echo times.

### Visualization
The code includes rich plotting of:
- Original MRI observations
- Signal amplitude (a)
- T2* relaxation time maps
- Mean signal intensity decay across TE

## Model Comparison: Levenberg–Marquardt vs ADMM-TV
We also implemented a baseline model using the Levenberg–Marquardt (LM) algorithm for voxel-wise T2* estimation via non-linear least squares fitting. Unlike our ADMM-based method, which incorporates spatial total variation (TV) regularization, LM fits each pixel independently and does not exploit spatial coherence. While LM performs well in high-SNR regions, our ADMM-TV approach produces smoother and more anatomically consistent T2* maps, especially in noisy or clinical data.
