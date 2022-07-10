########################
#1. Setup R Environment#
########################
## Clear memory
rm(list = ls()) 
options(scipen = 999)

## Set Working Directory and reference paths
`%+%` <- function(a, b) paste0(a, b)
data_dir <- getwd() %+% "/Data/"
out_dir <- getwd() %+% "/Output/"
script_dir <- getwd() %+% "/Scripts/"

#Pull in Custom Functions
source(script_dir %+% "00_twd_functions.R")

## Load Libraries
libs <- c("tidyverse", "coda", "ggmcmc", "gridExtra")
load_packages(libs)

######################################
#2. Model #1 Preparation and Summary #
######################################

## Read in MCMC Results#
coda_file <-  out_dir %+% "jags_samples_Full_mod1_N200_b1_2010-18.rds" 

jags_samples <- readRDS(coda_file)

# Uncondition Beta
gg_sample <- twd_construct_ggs(jags_samples,"beta") %>% twd_uncondition_beta()
thick_ci <- c(0.25, 0.75)
thin_ci <-  c(0.025, 0.975)
beta_ci <- ci(gg_sample,thin_ci = thin_ci, thick_ci = thick_ci) 
colnames(beta_ci) <- c("Parameter","2.5%","25%","50%","75%","97.5%")
beta_mean <- gg_sample %>% 
  group_by(Parameter) %>%
  summarize(Mean = mean(value), SD = sd(value)) %>% 
  ungroup()

beta_stats <- inner_join(beta_mean,beta_ci)

# Reparameterize Recruited (by year not cumulative)
rec_chain <- do.call(rbind, jags_samples)[,"recruited[" %+% 1:9 %+% "]"]
rec_new <- matrix(nrow = nrow(rec_chain), ncol = ncol(rec_chain))

rec_new[,1] <- rec_chain[,1] #First year is the same

for(i in 2:9){ # Transform recruited results accordingly (difference by year)
  rec_new[,i] <- rec_chain[,i] - rec_chain[,(i-1)]
}

rec_mean <- tibble(Parameter = colnames(rec_chain), 
       Mean = colMeans(rec_new),
       SD = apply(rec_new,2,sd))

rec_ci <- ci(ggs(mcmc(rec_new)),thin_ci = thin_ci, thick_ci = thick_ci) %>% 
  mutate(Parameter = "recruited[" %+% 1:9 %+% "]")

colnames(rec_ci) <- c("Parameter","2.5%","25%","50%","75%","97.5%")

rec_stats <- inner_join(rec_mean,rec_ci)

# Produce Tables
summ <- summary(jags_samples)

param_mean <- as_tibble(summ$statistics) %>% 
  mutate(Parameter = rownames(summ$statistics)) %>% 
  select(Parameter, Mean, SD)

param_ci <- as_tibble(summ$quantiles) %>% 
  mutate(Parameter = rownames(summ$quantiles)) %>% 
  select(Parameter, "2.5%","25%","50%","75%", "97.5%")

stats <- inner_join(param_mean,param_ci) %>% 
  filter(substr(Parameter,1,5) != "beta[" & 
         substr(Parameter,1,10) != "recruited[") 

stats <- rbind(stats,rec_stats,beta_stats)

# Write XLSX files 
write.csv(stats, getwd() %+% "/Mod1_Summary.csv")

######################################
#2. Model #1 Preparation and Summary #
######################################

## Read in MCMC Results#
coda_file <-  out_dir %+% "jags_samples_Full_mod2_N200_b1_2010-18.rds" 

jags_samples <- readRDS(coda_file)

# Uncondition Beta
gg_sample <- twd_construct_ggs(jags_samples,"beta") %>% twd_uncondition_beta()
thick_ci <- c(0.25, 0.75)
thin_ci <-  c(0.025, 0.975)
beta_ci <- ci(gg_sample,thin_ci = thin_ci, thick_ci = thick_ci) 
colnames(beta_ci) <- c("Parameter","2.5%","25%","50%","75%","97.5%")
beta_mean <- gg_sample %>% 
  group_by(Parameter) %>%
  summarize(Mean = mean(value), SD = sd(value)) %>% 
  ungroup()

beta_stats <- inner_join(beta_mean,beta_ci)

# Reparameterize Recruited (by year not cumulative)
rec_names <- "recruited" %+% c("",".of",".other")
rec_stat_list <- list()
rec_chain_list <- list()

for(j in 1:3){
  rec_name <- rec_names[j]
  rec_chain <- do.call(rbind, jags_samples)[,rec_name %+% "[" %+% 1:9 %+% "]"]
  rec_new <- matrix(nrow = nrow(rec_chain), ncol = ncol(rec_chain))
  
  rec_new[,1] <- rec_chain[,1] #First year is the same
  
  for(i in 2:9){ # Transform recruited results accordingly (difference by year)
    rec_new[,i] <- rec_chain[,i] - rec_chain[,(i-1)]
  }
  
  rec_mean <- tibble(Parameter = colnames(rec_chain), 
                     Mean = colMeans(rec_new),
                     SD = apply(rec_new,2,sd))
  
  rec_ci <- ci(ggs(mcmc(rec_new)),thin_ci = thin_ci, thick_ci = thick_ci) %>% 
    mutate(Parameter = rec_name %+% "[" %+% 1:9 %+% "]")
  
  colnames(rec_ci) <- c("Parameter","2.5%","25%","50%","75%","97.5%")
  rec_stats <- inner_join(rec_mean,rec_ci)
  
  rec_chain_list[[j]] <- rec_new
  rec_stat_list[[j]] <- rec_stats
}

# Produce Tables
summ <- summary(jags_samples)

param_mean <- as_tibble(summ$statistics) %>% 
  mutate(Parameter = rownames(summ$statistics)) %>% 
  select(Parameter, Mean, SD)

param_ci <- as_tibble(summ$quantiles) %>% 
  mutate(Parameter = rownames(summ$quantiles)) %>% 
  select(Parameter, "2.5%","25%","50%","75%", "97.5%")

stats <- inner_join(param_mean,param_ci) %>% 
  filter(substr(Parameter,1,4) != "beta" & 
           substr(Parameter,1,9) != "recruited") 

stats <- rbind(stats,rec_stat_list[[1]],rec_stat_list[[2]],rec_stat_list[[3]],beta_stats)

# Write XLSX files 
write.csv(stats, getwd() %+% "/Mod2_Summary.csv")

