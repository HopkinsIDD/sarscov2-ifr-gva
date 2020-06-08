# Preamble ---------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(gridExtra)

source("scripts/utils.R")

remove_EMS <- F     # account for deaths in nursing homes separately
last_survey <- F    # use only last survey
redo_fit <- T

suffix <- ifelse(last_survey, "_lastsurvey", "")  # suffix for result file names

# Load data setup (this is saved in the main script)
if (remove_EMS) {
  load("data/data_setup_with_EMS.rda")
} else {
  load("data/data_setup.rda")
}

# Load restults
if(remove_EMS) {
  IFRs <- readRDS(paste0("results/age_stratified_IFRs_with_EMS", suffix, ".rds"))
} else {
  IFRs <- readRDS(paste0("results/age_stratified_IFRs", suffix, ".rds"))
}

# Post-stratification of IFR estimates -----------------------------------------

IFRs <- IFRs %>%
  group_by(age_class) %>% 
  mutate(sim = row_number()) %>% 
  ungroup() %>% 
  inner_join(prev_est_age %>% 
               filter(week == 5),
             by = c("sim", "age_class" = "age_cat")) %>% 
  group_by(sim) %>% 
  summarise(ifr = weighted.mean(ifr, seropos * pop)) %>% 
  mutate(age_class = "all") %>% 
  select(ifr, age_class) %>% 
  rbind(IFRs)

IFRs <- IFRs %>% 
  mutate(age_class = factor(age_class, 
                            levels = age_classes[c(3, 1, 2, 4:length(age_classes))]))
  
# Compute Bayesian p-values ----------------------------------------------------

ref_age <- "[20,50)"    # reference age class
ref_ifrs <- IFRs$ifr[IFRs$age_class == ref_age] # reference IFR posterior draws

IFR_pvals <- IFRs %>% 
  group_by(age_class) %>% 
  summarise(pval = computePval(ifr, ref_ifrs)) %>% 
  mutate(pval = case_when(age_class == ref_age ~ as.double(NA), T ~ pval),
         pval_string = case_when(pval == 0 ~ "<0.001",
                                 is.na(pval) ~ "-",
                                 T ~ format(pval, digits = 2)))

# Plots ------------------------------------------------------------------------
Sys.setlocale("LC_TIME", "C")


# - - - - 
# Plot epi data

var_dict <- c(
  "case_cumul" = "Cumulative incidence",
  "death_cumul" = "Cumulative deaths",
  "case_incid" = "Case incidence",
  "death_incid" = "Death incidence"
)

stratified_data <- read_csv("data/stratified_dgs_data.csv")

p_data <- ggplot(stratified_data, aes(x = date, y = value, color = age_class)) +
  geom_point(size = .6) +
  geom_line() +
  facet_grid(var~EMS, scales = "free", labeller =  labeller(var = var_dict, EMS = label_both)) +
  theme_bw() +
  labs(x = "", y = "") +
  guides(color = guide_legend("Age class"))

ggsave(p_data, filename = "results/fig/epidata.png", width = 8, height = 7)


# - - - - 
# Plot delay distributions

# Delays for plotting
nmax <- 50  # max number of days in PDFS
dt <- .1
ntot <- nmax/dt
n_dist_samples <- 1000
# Initialize PDF matrices
pdf_inc <- matrix(rep(0, ntot * n_dist_samples), nrow = n_dist_samples)
pdf_report <- pdf_inc
pdf_case <- pdf_inc
pdf_symp_sero <- pdf_inc
pdf_sero <- pdf_inc
pdf_report_death <- pdf_inc
pdf_symp_death <- pdf_inc
pdf_death <- pdf_inc

# Sample distributions
for (i in 1:n_dist_samples) {
  pdf_inc[i, ] <- getDelay(delay_params, "inc", nmax, dt = dt, rnd = T)*dt
  pdf_report[i, ] <- getDelay(delay_params, "report", nmax, dt = dt, rnd = T)*dt
  pdf_case[i, ] <- convolve( pdf_inc[i, ], rev(pdf_report[i, ]), type = "o")[1:ntot]
  pdf_symp_sero[i, ] <- getDelay(delay_params, "symp_sero", nmax, dt = dt, rnd = T)*dt
  pdf_sero[i, ] <- convolve(pdf_inc[i, ], rev(pdf_symp_sero[i, ]), type = "o")[1:ntot]  # convolve with incubation period
  pdf_report_death[i, ] <- getDelay(delay_params, "report_death", nmax, dt = dt, rnd = T)*dt
  pdf_symp_death[i, ] <- convolve(pdf_report[i, ], rev(pdf_report_death[i, ]), type = "o")[1:ntot]  # convolve with delay from symptoms to hospitalization
  pdf_death[i, ] <- convolve(pdf_inc[i, ], rev(pdf_symp_death[i, ]), type = "o")[1:ntot]  # convolve with incubation period
}

delay_dict <- c(
  "i2s" = "Infection to symptom onset",
  "s2c" = "Onset to reporting",
  "s2sero" = "Onset to seroconversion",
  "c2d" = "Confirmation to death",
  "i2c" = "Infection to reporting",
  "i2sero" = "Infection to seroconversion",
  "i2d" = "Infection to death"
)

legend_title <- "Delay distribution"
delays_2plot <- map2_dfr(
  list(pdf_inc, pdf_report, pdf_symp_sero, pdf_report_death, pdf_case, pdf_sero, pdf_death),
  c("i2s", "s2c", "s2sero", "c2d", "i2c", "i2sero", "i2d"), 
  ~getQuantiles(.x) %>%
    mutate(days = seq(dt/2, nmax, by = dt),
           var = .y)) %>% 
  mutate(type = case_when(
    var %in% c("i2s", "s2c", "s2sero", "c2d") ~ "literature",
    T ~ "derived"),
    var_name = delay_dict[var]) 

p_lit <- delays_2plot %>% 
  filter(type == "literature") %>% 
  ggplot(aes(x = days, y = cdf.median, ymin = cdf.q025, ymax = cdf.q975, fill = var_name, lty = var_name)) + 
  geom_ribbon(alpha = .25) +
  geom_line(aes(color = var_name)) +
  guides(color = guide_legend(legend_title),
         fill = guide_legend(legend_title),
         lty = guide_legend(legend_title)) +
  theme_bw() +
  labs(x = "time [days]", y = "CDF") +
  ggthemes::scale_color_few() +
  ggthemes::scale_fill_few() +
  theme(legend.position = c(.71, .2))

p_conv <- delays_2plot %>% 
  filter(type == "derived") %>% 
  ggplot(aes(x = days, y = cdf.median, ymin = cdf.q025, ymax = cdf.q975, fill = var_name, lty = var_name)) + 
  geom_ribbon(alpha = .25) +
  geom_line(aes(color = var_name)) +
  guides(color = guide_legend(legend_title),
         fill = guide_legend(legend_title),
         lty = guide_legend(legend_title)) +
  theme_bw() +
  labs(x = "time [days]", y = "CDF") +
  theme(legend.position = c(.71, .2))

p_comb <- arrangeGrob(grobs = list(p_lit, p_conv), nrow = 1)
ggsave(plot = p_comb, filename = "results/fig/delay_dist.png", width = 9, height = 4.5)


# - - - - 
# Plot IFR posterior draws

if (!remove_EMS) {
  p_ifrs <- IFRs %>%
    ggplot(aes(x = ifr*100)) +
    geom_histogram() +
    facet_wrap(~age_class, scales = "free") +
    xlab("IFR [%]") +
    theme_bw()
  
  ggsave(p_ifrs, filename = "results/fig/ifr_posteriors.png", width = 8, height = 4)
}

# Table ------------------------------------------------------------------------

IFR_stats <- 
  # Epi data
  age_epidata %>%
  group_by(age_class) %>% 
  arrange(desc(case_cumul)) %>% 
  select(case_cumul, death_cumul) %>% 
  slice(1) %>% ungroup() %>% 
  # Populations by age
  inner_join(age_popdata %>% select(age_class, pop)) %>% 
  # Seroprevalence estimates
  inner_join(prev_est_age %>% 
               rename(age_class = age_cat) %>% 
               filter(week == 5) %>% 
               group_by(age_class) %>% 
               summarise(sero.mean = mean(seropos),
                         sero.025 = quantile(seropos, 0.025),
                         sero.975 = quantile(seropos, 0.975))) %>% 
  mutate(seropop.mean = sero.mean * pop,
         seropop.025 = sero.025 * pop,
         seropop.975 = sero.975 * pop) %>% 
  # IFR estimates
  inner_join(
    IFRs %>% 
      group_by(age_class) %>% 
      mutate(ifr = ifr*1e2) %>% 
      summarise(mean = mean(ifr),
                q025 = quantile(ifr, 0.025),
                q975 = quantile(ifr, 0.975))
  ) %>% 
  inner_join(IFR_pvals) %>% 
  group_by(age_class) %>% 
  mutate_at(vars(mean, q025, q975), 
            function(x) format(x, digits = 2)) %>% 
  mutate_at(vars(contains("seropop")), 
            function(x) 100*round(x/100)) %>% 
  ungroup() %>% 
  mutate(ifr = paste0(mean, " (", q025, "-", q975, ")"),
         seropop = paste0(format(seropop.mean, big.mark = "'"),
                          " (", format(seropop.025, big.mark = "'"), 
                          "-", 
                          format(seropop.975, big.mark = "'"), ")"),
         age_class = factor(age_class, 
                            levels = age_classes[c(3, 1, 2, 4:length(age_classes))])) %>% 
  arrange(age_class) %>% 
  mutate(age_class = as.character(age_class),
         age_class = str_replace_all(age_class, c("\\[" = "", "\\)" = "", "," = "-")),
         pop = format(pop, big.mark = "'")) %>% 
  select(age_class, pop, seropop, death_cumul, ifr) %>% 
  rename(`Age class` = age_class,
         Population = pop,
         `Estimated infected` = seropop,
         `Deaths` = death_cumul,
         IFR = ifr)

if(remove_EMS) {
  knitr::kable(IFR_stats,
               format = "latex",
               booktabs = T) %>% 
    write(file = paste0("results/age_stratified_IFRs_with_EMS_stats", suffix, ".tex"))
    
} else {
  knitr::kable(IFR_stats,
               format = "latex",
               booktabs = T) %>% 
    write(file = paste0("results/age_stratified_IFRs_stats", suffix, ".tex"))
}
