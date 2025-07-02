# MAFT_TD

Gaussian Mixture Accelerated Failure Time (AFT) Model with Time-dependent Covariates

## Description

This repository contains an R implementation of the `MAFT_TD()` function, which fits an Accelerated Failure Time model using a **nonparametric Gaussian scale mixture** approach. It accommodates **time-dependent covariates** via the start–stop counting process formulation.

## Main Features

- Support for **time-dependent covariates**
- Assumes **a nonparametric Gaussian scale mixture distribution** for the log-transformed failure time (log(T₀))
- Efficient optimization via **directional derivative-based support update**
- **EM-type algorithm** for coefficient estimation
- **Bootstrapped standard errors** available

## Files

- `MAFT_TD.R` – Main function including all supporting routines
- `MAFT_TD.Rproj` – RStudio project file
- `Normal_b111.csv` – Example data simulated under log-normal model with true coefficients (intercept = 1, Z1 = 1, X = 1)

## Installation & Example Usage

```r
# Install required packages
install.packages(c("dplyr", "survival", "nloptr", "nnls"))

# Load the main function
source("MAFT_TD.R")

# Load example data
dat <- read.csv("Normal_b111.csv")

# Fit the model
result <- MAFT_TD(dat = dat, X = c("Z1", "X"), subject_id = "ID", process = TRUE)

# View results
result$Final_Est
```

## License & Author

This software is licensed for academic and research use only.  
Copyright (c) 2025  
Department of Applied Statistics, Yonsei University  
Department of Statistics, Yeungnam University  
All rights reserved. Redistribution or commercial use is prohibited without explicit permission.  

**Ju-young Park**  
Assistant Professor, Department of Statistics, Yeungnam University  
GitHub: [JYSTAT](https://github.com/JYSTAT)
