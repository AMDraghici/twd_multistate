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

## Read in MCMC Results#
coda_file <-  out_dir %+% "jags_samples_Full_mod1_N200_b1_2010-18.rds"

jags_samples <- readRDS(coda_file)
jags_params <- c("xi","beta","phi","p","psi","abundance")

## Plot Effective Sample Sizes
twd_effective_size(jags_samples)

# Process MCMC Results with respect to parameter
gg_sample1 <- twd_construct_ggs(jags_samples,"phi.xi") %>% 
  mutate(year = as.integer(substr(Parameter,8,str_length(Parameter)-1))) %>% 
  select(-Parameter) %>% 
  rename("phi.xi" = value)
gg_sample2 <- twd_construct_ggs(jags_samples,"phi.tau") %>% 
  mutate(year = as.integer(substr(Parameter,9,str_length(Parameter)-1))) %>% 
  select(-Parameter) %>% 
  rename("phi.tau" = value)

############################
# Diagnose MCMC Convergence#
############################

## Summary Statistics
mask <- as.character(unique(gg_sample$Parameter))
summary(jags_samples[,mask])

## Compute/Plot Effective Sampling Rate
ggs_effective(gg_sample)

## Trace and Density
ggs_traceplot(gg_sample) #+ ylim(0,0.25)
ggs_density(gg_sample) + xlim(-20,20)
ggs_pairs(gg_sample,lower = list(continuous = "density")) #SLOW!


## Barplots/Histograms for Discrete Variables
gg_sample %>% ggplot() +
  geom_bar(mapping = aes(x = value,fill=factor(Chain)),position="dodge") + 
  facet_grid(Parameter~.)

gg_sample %>% 
  ggplot(mapping = aes(x=value,fill=factor(Chain))) +
  geom_bar() + facet_grid(Parameter~Chain) +
  labs(fill = "Chain") +
  theme(legend.position = "bottom")

#ggs_histogram(gg_sample) ugly - make own 

## Autocorrelation 
ggs_autocorrelation(gg_sample)

## Running Means
ggs_running(gg_sample)

## Partial Density
#ggs_compare_partial(gg_sample)

## Crosscorrelation
ggs_crosscorrelation(gg_sample)

## Check R-hat for between chain variation 
ggs_Rhat(gg_sample) + 
  xlab("R_hat")

# Check geweke for comparing initial and final chain means 
ggs_geweke(gg_sample)

## Inspect Caterpillar
gg_sample %>% 
  ggs_caterpillar(horizontal = T) + 
  xlim(0,0.3)

# ## Investigate Probability Spread
gg_sample %>%
  ggplot(aes(x=Iteration,
             y = value,
             color=Parameter)) +
  geom_line() #+ ylim(0,0.5)

## Inspect Caterpillar #2
X <- data.frame(Parameter = unique(pull(gg_sample,Parameter))) %>% 
  rownames_to_column() %>% 
  rename(x = rowname) %>% 
  select(Parameter,x)

jags_samples %>%
  twd_caterpillar(x=X, horiz= T) + 
  ylim(0,0.3)

## Export Results to a PDF
dev.copy2pdf(file="Figures/eff_red1.pdf")
