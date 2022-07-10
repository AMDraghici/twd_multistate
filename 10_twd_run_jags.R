########################
#1. Setup R Environment#
########################

## Clear memory and apply settings
rm(list = ls()) 
options(scipen = 999)

## Define Paths and read-in custom functions
`%+%` <- function(a,b) paste0(a,b)
data_dir <- getwd() %+% "/Data"
script_dir <- getwd() %+% "/Scripts"
out_dir <- getwd() %+% "/Output" 
source(script_dir %+% "/00_twd_functions.R")

## Load Libraries
libs <- c("tidyverse", "openxlsx", "dclone", "lubridate", "readxl",
          "abind", "coda", "rjags","gridExtra" ,"parallel")
load_packages(libs)

#############################
#2. Extract and process data#
#############################

## Read in and process data
filters <- list(years = 2010:2018,
                months = c(6,7,8,9),
                blocks = 2:15,
                days = Inf)

filenames <- c("Individual_sightings_block_UPDATED.csv",
               "Block_location.csv",
               "Survey_information_UPDATED.csv",
               "old_females_UPDATED.csv")

jags_data <- twd_create_jags_data(data_dir, filenames, filters, Nmax = 200) 

## Save/read processed full data 
saveRDS(jags_data, out_dir %+% "/twd_data_full_200_b1_2010-18.rds") #write data

#########################
# 3. Run JAGS (Parallel)#
#########################

## MCMC parameters  
par_settings <- list('n.iter' = 3e5, 
                     'n.thin' = 200,
                     'n.burn' = 1e5,
                     'n.chains' = 4,
                     'n.adapt' = 1e5,
                     'psi_setting' = "random") 

## Jags parameters and model script

# Run Full Model + No Groups
jags_params1 <- c("xi","beta","phi","p","psi","abundance","recruited","phi.xi","phi.tau","p.xi","p.tau")
jags_model1 <- script_dir %+% "/01_twd_jags_model.R"

## Run jags in parallel and save results
jags_samples1 <- twd_run_jags_parallel(jags_data, 
                                       jags_model1,
                                       jags_params1, 
                                       par_settings,
                                       out_dir,
                                       outname = "Full_mod1_N200_b1_2010-18")


#Run with Old Females
jags_params2 <- c("xi","beta","theta",
                 "phi","phi.of","p","p.of",
                 "phi.xi","phi.tau","p.xi","p.tau",
                 "phi.xi.of","phi.tau.of","p.xi.of","p.tau.of",
                 "psi","psi.of",
                 "abundance","abundance.of","abundance.other",
                 "recruited","recruited.of","recruited.other",
                 "OF")

jags_model2 <- script_dir %+% "/02_twd_jags_model2.R"

## Run jags in parallel and save results
jags_samples2 <- twd_run_jags_parallel(jags_data, 
                                       jags_model2, 
                                       jags_params2,
                                       par_settings, 
                                       out_dir, 
                                       outname = "Full_mod2_N200_b1_2010-18")

#Run Full Model 3 + No Groups
jags_params3 <- c("xi","beta","phi","p","psi","abundance","recruited")
jags_model3  <- script_dir %+% "/03_twd_jags_model3.R"

# ## Run jags in parallel and save results
jags_samples <- twd_run_jags_parallel(jags_data, 
                                      jags_model3,
                                      jags_params3, 
                                      par_settings,
                                      out_dir,
                                      outname = "Full_mod1_N200_b1_2010-18")