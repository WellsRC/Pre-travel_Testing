# Testing for COVID-19 Is Much More Effective when Performed Immediately Prior to Social Mixing: MATLAB Code

## Wells C. R., Gokcebel S., Pandey A., Galvani A.P., Townsend J.P.
 

Copyright (C) <2022>, Chad R. Wells et. al. All rights reserved. Released under the GNU General Public License (GPL)

This repository contains codes and data used to simulate and analyze COVID- pre-arrival testing strategies in the scenarios of different variants

The model code is written in MATLAB and results are saved as MATLAB data files (extension .mat), with plots also being constructed in MATLAB. 

## OS System requirements
The codes developed here are tested on Windows operating system (Windows 10 Home: 64-bit). 

## Installation guide
### MATLAB
Installation instruction for MATLAB can be found at https://www.mathworks.com/help/install/install-products.html. Typical install time for MATLAB on a "normal" desktop is around 30-40 minutes. The current codes were developed and tested on MATLAB R2019b.

## Analysis functions
TX0 - Runs the MLE and uncertainty analysis for RT-PCR

NT- Runs the analysis for no tests

TX1 - Runs the MLE analysis for RA tests

TX1Fa - Runs a subset of uncertainty analysis for the RA tests

TX1Fb - Runs a subset of uncertainty analysis for the RA tests

## Result functions
FigureS4 - Generates Figure S4

Table_Output - Generates tables for supplement

Table_Beta_Coefficients- Generates table for PPA logistic regression coefficients

Plot_PPA_RA_Sensitivity- Produces plots for the different RA tests in the supplement

Plot_PCR_Fit - Generates the plot for the RT-PCR diagnostic sensitivity

Plot_NB_Probability- Produces the plot for the sensitivty analysis for the ngative binomial 

