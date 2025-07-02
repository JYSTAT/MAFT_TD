# MAFT_TD

Gaussian Mixture Accelerated Failure Time (AFT) Model with Time-dependent Covariates

## Description

This repository contains an R implementation of the `MAFT_TD()` function, which fits an Accelerated Failure Time model using a **nonparametric Gaussian scale mixture** approach. It accommodates **time-dependent covariates** via the start–stop counting process formulation.

## Main Features

- Support for **time-dependent covariates**
- **Nonparametric Gaussian scale mixtures** for error distribution
- Efficient optimization via **directional derivative-based support update**
- **EM-type algorithm** for coefficient estimation
- **Bootstrapped standard errors** available

## Files

- `MAFT_TD.R` – Main function including all supporting routines
- `MAFT_TD.Rproj` – RStudio project file

## Installation & Example Usage

```r
# Install required packages
install.packages(c("dplyr", "survival", "nloptr", "nnls"))

# Load the function after cloning the repository
source("MAFT_TD.R")

# Fit the model with your start–stop format survival data
result <- MAFT_TD(dat = your_data, X = c("Z1", "X1"), process = TRUE)

# Extract estimated regression coefficients
result$Final_Est

## License & Author

```r
This software is licensed for academic and research use only.  
Copyright (c) 2025  
Department of Applied Statistics, Yonsei University  
Department of Statistics, Yeungnam University  
All rights reserved. Redistribution or commercial use is prohibited without explicit permission.  

**Ju-young Park**  
Assistant Professor, Department of Statistics, Yeungnam University  
GitHub: [JYSTAT](https://github.com/JYSTAT)
