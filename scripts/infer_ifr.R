# Preamble -------------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(rstan)
library(foreach)

options(mc.cores = 5)
rm(list = ls())

source("scripts/utils.R")

remove_EMS <- F     # account for deaths in nursing homes separately
last_survey <- F    # use only last survey
redo_fit <- T

suffix <- ifelse(last_survey, "_lastsurvey", "")  # suffix for result file names

# Epi data ---------------------------------------------------------------------
age_epidata <- read_csv("data/stratified_dgs_data.csv") 

if (remove_EMS) {
  # Account for EMS cases and deaths separately
  age_epidata <- age_epidata %>% 
    mutate(age_class = case_when(EMS == 1 ~ "EMS", T ~ age_class)) %>% 
    select(-EMS) 
} 

age_epidata <- age_epidata %>% 
  group_by(date, age_class, var) %>% 
  summarise(value = sum(value)) %>% 
  ungroup() %>% 
  arrange(age_class, date) %>% 
  pivot_wider(id_cols = c("age_class", "date"),
              names_from = "var",values_from = "value")

# Add whole population
age_epidata <- age_epidata %>%
  rbind(age_epidata %>%
          select(-age_class) %>%
          group_by(date) %>%
          summarise_all(sum) %>%
          mutate(age_class = "all"))


# Population data --------------------------------------------------------------
age_popdata <- read_csv("data/stratified_pop_data.csv") # data from: http://www.ge.ch/statistique/tel/domaines/01/01_01/T_01_01_8_01.xls [Accesed June 3 2020]

if (remove_EMS) {
  EMS_pop <- 4065 # from https://www.ge.ch/statistique/tel/publications/2019/informations_statistiques/autres_themes/is_etablissements_sante_01_2019.pdf [Accesed June 3 2020]
  # Remove the population in the EMS
  age_popdata$pop[age_popdata$age_class == "[65, Inf)"] <- age_popdata$pop[age_popdata$age_class == "[65, Inf)"] - EMS_pop
  # Add EMS poplation
  age_popdata <- rbind(age_popdata, tibble(age_class = "EMS", pop = EMS_pop))
}

#  Seroprevalence estimates ----------------------------------------------------
prev_est_age <- read_csv("data/age-sex-week-est.csv") %>% 
  select(-X1) %>% 
  group_by(week, age_cat, sim) %>% 
  # Aggregat over age classes
  summarise(seropos = weighted.mean(seropos, pop),
            pop = sum(pop)) %>%
  ungroup() %>% 
  mutate(age_cat = str_replace_all(age_cat, "105", " Inf"))

# Add whole population seroprev estimates
prev_est_age <- rbind(prev_est_age,
                      prev_est_age %>% 
                        group_by(week, sim) %>% 
                        summarise(seropos = weighted.mean(seropos, pop),
                                  pop = sum(pop)) %>% 
                        ungroup() %>% 
                        mutate(age_cat = "all")
)

if (remove_EMS) {
  # Assume EMS seroprevalence is as high as the 20-50 age class
  prev_est_age <- rbind(
    mutate(prev_est_age, pop = case_when(age_cat == "[65, Inf)" ~ pop - EMS_pop, T ~ pop)), 
    filter(prev_est_age, age_cat == "[20,50)") %>% 
      mutate(age_cat = "EMS", 
             pop = EMS_pop))
}

if (last_survey) {
  prev_est_age <- filter(prev_est_age, week == 5)
}

age_classes <- unique(prev_est_age$age_cat)    # unique age classes
n_sero_weeks <- length(unique(prev_est_age$week))  # number of seroseruvey weeks
n_post <- length(unique(prev_est_age$sim))   # number of posterior samples

# This is only used for serosurvey dates
serodata <- tribble(
  ~date, ~median, ~q025, ~q975,
  "2020-04-06", 2.9, 0.5, 6.5,
  "2020-04-16", 6.7, 3.1, 11,
  "2020-04-21", 9.4,  4.7, 15.2,
  "2020-04-29", 9.4,  4.7, 15.2,
  "2020-05-06", 10.8, 8.2, 13.9
) %>% 
  mutate(date = as.Date(date))

if (last_survey) {
  serodata <- slice(serodata, 5)
}

# Delay distributions ----------------------------------------------------------

delay_params <- tribble(
  ~delay, ~logmu, ~logmu.sd, ~logmu.low, ~logmu.high, ~logsigma, ~logsigma.sd, ~logsigma.low, ~logsigma.high,
  "inc", 1.57, NA, 1.44, 1.69, 0.65, NA, 0.56, 0.73,     # Bi et al. 2020
  "report", 1.49, 0.065, NA, NA, 0.756, 0.0458, NA, NA,  # Scire et al. 2020
  "symp_sero", 2.34, 0.114, NA, NA, 0.38, 0.26, NA, NA,  # Stringhini et al. 2020
  "report_death", 2.1, 0.055, NA, NA, 0.87,  0.039, NA, NA # DGS
) %>% 
  group_by(delay) %>% 
  mutate(
    logmu.sd = case_when(is.na(logmu.sd) ~ getSD(logmu, logmu.low, logmu.high), T ~ logmu.sd),
    logsigma.sd = case_when(is.na(logsigma.sd) ~ getSD(logsigma, logsigma.low, logsigma.high),  T ~ logsigma.sd)
  ) %>% 
  ungroup() %>% 
  mutate(mean = lnormMean(logmu, logsigma),   # Mean and SD on natural scale
         sd = lnormSD(logmu, logsigma))

if (!remove_EMS) {
  save(list = ls(), file = "data/data_setup.rda")
} else {
  save(list = ls(), file = "data/data_setup_with_EMS.rda")
}

# Sample delay PDFs 
n_dist_samples <- n_post # number of PDF samples to draw from
nmax <- 100  # max number of days in PMFs

# Initialize PDF matrices
pdf_inc <- matrix(rep(0, nmax * n_dist_samples), nrow = n_dist_samples)
pdf_report <- pdf_inc
pdf_symp_sero <- pdf_inc
pdf_sero <- pdf_inc
pdf_report_death <- pdf_inc
pdf_symp_death <- pdf_inc
pdf_death <- pdf_inc

# Sample distributions
for (i in 1:n_dist_samples) {
  pdf_inc[i, ] <- setDelayPMF(delay_params, "inc", nmax, rnd = T)
  pdf_report[i, ] <- setDelayPMF(delay_params, "report", nmax, rnd = T)
  pdf_symp_sero[i, ] <- setDelayPMF(delay_params, "symp_sero", nmax, rnd = T)
  pdf_sero[i, ] <- convolve(pdf_inc[i, ], rev(pdf_symp_sero[i, ]), type = "o")[1:nmax]  # convolve with incubation period
  pdf_report_death[i, ] <- setDelayPMF(delay_params, "report_death", nmax, rnd = T)
  pdf_symp_death[i, ] <- convolve(pdf_report[i, ], rev(pdf_report_death[i, ]), type = "o")[1:nmax]  # convolve with delay from symptoms to hospitalization
  pdf_death[i, ] <- convolve(pdf_inc[i, ], rev(pdf_symp_death[i, ]), type = "o")[1:nmax]  # convolve with incubation period
}

# Inference --------------------------------------------------------------------

# Settings for stan
control <- list(adapt_delta = .9, max_treedepth = 12, metric = "dense_e")
nwarmup <- 4000
niter <- nwarmup + 1000
n_chains <- 5

IFRs <- foreach(agec = age_classes[age_classes != "all"],
                .combine = rbind) %do% 
  {
    if (agec == "[65, Inf)" & remove_EMS) {
      res_file <- paste0("results/age_stratified_IFRs_with_EMS_", agec, suffix, ".rds") 
      stanfit_file <- paste0("results/stanfi_age_stratified_IFRs_with_EMS_distsample_", agec, suffix, ".rds") 
    } else {
      res_file <- paste0("results/age_stratified_IFRs_distsample_", agec, suffix, ".rds") 
      stanfit_file <- paste0("results/stanfit_age_stratified_IFRs_distsample_", agec, suffix, ".rds") 
    }
    
    if (!file.exists(res_file) | redo_fit) { 
      
      # Compile stan model
      ifr_stan <- stan_model("scripts/ifr_betabinomial.stan")
      
      epidata <- filter(age_epidata, age_class == agec)
      pop_age <- age_popdata$pop[age_popdata$age_class == agec]  # population of the age class
      
      # Dates on which serosurveys were done
      ind_date <- which(epidata$date %in% serodata$date)
      
      # Deaths for model
      deaths <- epidata$death_cumul[ind_date]
      n_data <- nrow(epidata)
      
      # Seroprevalence estimates 
      # Posterior draws as matrix with each row corresponding to a given date, 
      # and each column to a sample. 
      thetas <- prev_est_age %>% 
        filter(age_cat == agec) %>% 
        select(week, seropos, sim) %>% 
        spread(sim, seropos) %>% 
        select(-week) %>% 
        as.matrix()
      
      # Posterior draws of seroprevalent population
      Isero <- pop_age * thetas
      
      # Initialize Istar and phi matrices [n data points x n samples]
      Istar <- matrix(rep(0, n_data * n_dist_samples), ncol = n_dist_samples)
      phi <- Istar
      
      # Loop over samples of delay distributions
      for (i in 1:n_dist_samples) {
        
        # Deconvolve state of cumulative infections up to proportionality (Istar)
        Istar[, i] <- computeIstar(epidata$case_cumul, pdf_inc[i, ], pdf_report[i, ])
        
        # Compute phi(t)
        phi[, i] <- (convolve(Istar[, i], rev(pdf_death[i, ]), type = "o")/
                       convolve(Istar[, i], rev(pdf_sero[i, ]), type = "o"))[1:n_data]
        
      }
      
      # Get phis for dates of scerosurveys
      phis <- phi[ind_date, ]
      
      # Compute number of Infected at risk of dying
      I <- round(Isero * phis) 
      
      # Set min of I to observed deaths to avoid errors in stan
      for (i in 1:nrow(I)) {
        I[i, I[i,] < deaths[i]] <- deaths[i]
      }
      
      # Data for stan model
      data <- list(N = n_sero_weeks, 
                   M = n_post, 
                   deaths = as.array(deaths), 
                   I = I)
      
      ifr_stanfit <- sampling(ifr_stan,
                              data = data,
                              init = rndInit, 
                              chains = n_chains,
                              warmup = nwarmup, 
                              iter = niter,
                              control = control,
                              save_warmup = FALSE)
      
      saveRDS(ifr_stanfit, file = stanfit_file)
      
      # Extract and save
      ifrs <- extract(ifr_stanfit, pars = "IFR")$IFR
      saveRDS(ifrs, res_file)
      
    } else {
      # Load results
      ifrs <- readRDS(res_file)
    }
    tibble(ifr = ifrs, age_class = agec)
  }

if(remove_EMS) {
  saveRDS(IFRs, paste0("results/age_stratified_IFRs_with_EMS", suffix, ".rds"))
} else {
  saveRDS(IFRs, paste0("results/age_stratified_IFRs", suffix, ".rds"))
}
