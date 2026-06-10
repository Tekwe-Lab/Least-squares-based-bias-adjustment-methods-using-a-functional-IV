# Least squares-based bias adjustment methods for scalar-on-function linear errors-in-variables models using a functional instrumental variable

Xiwei Chen 1, Ufuk Beyaztas 2, Caihong Qin 1, Heyang Ji 1, Gilson Honvoh3, Roger S. Zoh 1, Lan Xue 4, and Carmen D. Tekwe*1

1. Department of Epidemiology and Biostatistics, Indiana University, School of Public Health, Indiana, USA
2. Department of Statistics, Marmara University, Turkey
3. Department of Pediatrics, Cincinnati Children’s Hospital Medical Center, Ohio, USA
4. Department of Statistics, Oregon State University, Oregon, USA

---

**Abstract** \
Instrumental variables (IVs) are widely used to adjust for measurement error (ME) bias when assessing associations of health outcomes with ME-prone independent variables. IV approaches addressing ME in longitudinal models are well established, but few methods exist for functional regression. We develop two methods to adjust for ME bias in scalar-on-function linear errors-in-variables (EIV) models. We regress a scalar outcome on an ME-prone functional variable using a functional IV for model identification and propose two least squares–based methods to adjust for ME bias. Our methods alleviate potential computational challenges encountered when applying classical regression calibration methods for bias adjustment in high-dimensional settings and adjust for potential serial correlations across time. Simulations demonstrate faster run times, lower bias, and lower AIMSE for the proposed methods when compared to existing approaches. The proposed methods were applied to investigate the association between body mass index and wearable device-based physical activity intensity among community-dwelling adults living in the United States.

---

This repository contains R code accompanying the paper:

**“Least squares-based bias adjustment methods for scalar-on-function linear errors-in-variables models using a functional instrumental variable.”**

We provide:
- Data simulation code for scalar-on-function linear errors-in-variables (EIV) models with functional covariates
- Implementations of measurement error correction methods:
  - **MULTI-2SLS** (proposed)
  - **PW-2SLS** (proposed)
  - **SIMEX** (existing comparator, previously described by [Tekwe et al.](https://doi.org/10.1093/biostatistics/kxac017))
  - **Oracle** (benchmark using the true latent covariate; simulation-only)
  - **Naive** (ignores measurement error)
- Performance metric computations for the five estimators:
  - **ABias<sup>2</sup>**, **AVar**, **AIMSE**, **MSPEE**
- Simulation example code to assess sample-size effects
- NHANES (2005–2006 cycle) data-processing code used in the real-data application  
  All raw data files were obtained from the NHANES 2005–2006 cycle through the [U.S. CDC website](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2005).  
  The following component files were downloaded:
  - Demographic Variables & Sample Weights: [DEMO_D](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DEMO_D.xpt)
  - Body Measures: [BMX_D](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/BMX_D.xpt)
  - Diabetes: [DIQ_D](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/DIQ_D.xpt)
  - Physical Activity Monitor: [PAXRAW_D](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/PAXRAW_D.zip)
  
  Running the code available in the accompanying GitHub repository on the publicly available NHANES files reproduces the analytic dataset used in this study.

---

## Methods implemented

### Proposed bias adjustment methods
- **MULTI-2SLS**: Multivariate two-stage least squares approach.
- **PW-2SLS**: Pointwise two-stage least squares approach.

### Baselines / comparators
- **SIMEX**: Simulation-extrapolation method for measurement error correction.
- **Naive**: Fits the model ignoring measurement error.
- **Oracle**: Uses the true latent functional covariate (simulation benchmark; not available for real data).

---

## Repository structure

This repo contains three main R scripts:

1. **Core functions**  
   - Data simulation function(s)  
   - Measurement error procedures: MULTI-2SLS, PW-2SLS, SIMEX, Oracle, Naive  
   - Performance metrics: ABias<sup>2</sup>, AVar, AIMSE, MSPEE  
   **File:** `R/<Simulation_LinearFuncIV>.R`

2. **Simulation example: sample size effect**  
   - Demonstrates how to call the core functions to run simulation studies across varying sample sizes  
   - Produces summary tables used to evaluate estimator behavior  
   **File:** `R/<Example_DataSimulation>.R`

3. **NHANES (2005–2006) data preparation**  
   - Code to reproduce the analytic dataset used in this study from the publicly available NHANES 2005–2006 data, including merging, cleaning, and variable derivation procedures  
   **File:** `R/<Application_NHANES2005_DataPrep>.R`

---

## Requirements

- R (>= 4.0 recommended)

---
