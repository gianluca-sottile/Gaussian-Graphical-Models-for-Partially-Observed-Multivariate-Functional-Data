# Gaussian Graphical Models for Partially Observed Multivariate Functional Data

This repository contains the implementation of the methodology described in **"Gaussian Graphical Models for Partially Observed Multivariate Functional Data"**.

The project provides R and Rcpp code for estimating functional Gaussian graphical models in the presence of partially observed multivariate functional data. The estimation algorithm is based on an EM approach with a joint graphical lasso at the M-step.

---

## General structure and requirements

The general implementation is written in **R**, while all computationally intensive routines are implemented in **C++** using the **Rcpp** and **RcppArmadillo** frameworks.  

All required libraries are loaded in the preamble of each R script.  

Before running the code, ensure that your R environment is properly configured to compile C++ code through Rcpp (a C++ compiler compatible with your R installation is required).

---

## Repository organization

The repository is organized into the following directories:

### **1. Code/**
Contains the main source code for the proposed methodology: 

- **`pofggm.R`**
- **`pofggm.cpp`**.

---

### **2. Simulations/**
Contains scripts and results for the Monte Carlo simulation study described in the manuscript.

- **`figures/`**  
Includes graphical summaries of the simulation results (e.g., \(\text{MSE}_X\), \(\text{MSE}_\Theta\), and AUC plots).  

---

### **3. Analysis/**
Includes the scripts used to process simulation outputs and reproduce the figures reported in the manuscript.

- **`DATA/`**  
  Folder containing the original datasets from the paper Liebl (2019).  

- **`figs/`**  
  Includes figures shown in Section 5 of the paper.  

---

## How to run the code
