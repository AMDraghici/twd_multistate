# New abundance and survival estimates for the critically endangered Taiwanese white dolphin indicate no signs of recovery
## Authors: Claryana Araújo-Wang, John Y. Wang, Alexandru Marian Draghici, Peter S. Ross, Simon J. Bonner

The `R` code in this repository was used to conduct the statistical analysis for the article "New abundance and survival estimates for the critically endangered Taiwanese white dolphin indicate no signs of recovery". 

### 00_twd_functions.R

This file contains functions used to process, model, and draw inference from the Taiwanese White Dolphin data (not included in this repo) studied in this manuscript. There are also some functions used to preliminary studies that informed our model development and investigation. 

### 01_twd_jags_model.R/02_twd_jags_model.R/03_twd_jags_model.R

These three files contain JAGS (software which utilizes MCMC in order to train Bayesian models, see: https://mcmc-jags.sourceforge.io/) code which we used to model the data in this study. 

### 10_twd_run_jags.R

Script to call the jags model files using the functions in 00. 

### 20_twd_mcmc_diag.R/21_twd_generate_summary.R

Scripts used to investigate convergence of MCMC samplers and summarize estimates of the demographic parameters in our study. 

Article Link: https://doi.org/10.1002/aqc.3831

Reference: Araújo-Wang, C., Wang, J.Y., Draghici, A.M., Ross, P.S. & Bonner, S.J. (2022). New abundance and survival estimates for the critically endangered Taiwanese white dolphin indicate no signs of recovery. Aquatic Conservation: Marine and Freshwater Ecosystems, 1– 10. https://doi.org/10.1002/aqc.3831
