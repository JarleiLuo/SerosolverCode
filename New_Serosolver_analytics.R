
library(serosolver)

library(rngtools)
library(doRNG)

library(plyr)
library(data.table)

library(foreach)
library(parallel)

library(magrittr)
library(iterators)
library(doParallel)

library(viridisLite)
library(viridis)

library(reshape2)
library(bayesplot)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyverse)

library(readxl)


set.seed(123)

filename <- "small_new_ancvic_lastyear_(1,1)"
setwd(paste0("/Users/luojialei/Desktop/BSP/Flu B/new/",filename))
resolution <- 1
serosolver::describe_priors()
prior_version <- 2

antigenic_map <- read_excel("/Users/luojialei/Desktop/BSP/Flu B/antigenic_maps/lastyear_ancvic_antigenic_map.xlsx")
antibody_data <- read_excel("/Users/luojialei/Desktop/BSP/Flu B/new/new_titre_dat_ancvic.xlsx")
antibody_data <- as.data.frame(antibody_data)
antigenic_map <- as.data.frame(antigenic_map)


names(antibody_data)[names(antibody_data) == "DOB"] <- "birth"
names(antibody_data)[names(antibody_data) == "group"] <- "population_group"
names(antibody_data)[names(antibody_data) == "run"] <- "repeat_number"
names(antibody_data)[names(antibody_data) == "samples"] <- "sample_time"
names(antibody_data)[names(antibody_data) == "titre"] <- "measurement"
names(antibody_data)[names(antibody_data) == "virus"] <- "biomarker_id"


strains_isolation_times <- unique(antigenic_map$inf_times)

par_tab <- read.csv("/Users/luojialei/Downloads/par_tab(in).csv")
par_tab <- as.data.frame(par_tab)


par_tab[par_tab$names == "infection_model_prior_shape1","values"] <- 1
par_tab[par_tab$names == "infection_model_prior_shape2","values"] <- 1

antibody_data <- antibody_data %>% left_join(antibody_data %>% select(individual, population_group) %>% distinct() %>% mutate(new_id = 1:n())) %>% select(-individual) %>% rename(individual = new_id)


unique_indiv <- antibody_data[!duplicated(antibody_data$individual),]
ageMask <- create_age_mask(unique_indiv$DOB, strains_isolation_times)



start_inf <- setup_infection_histories(antibody_data, strains_isolation_times)

mcmc_pars <- c("iterations"=5000,"adaptive_frequency"=100,"thin"=1,
               "adaptive_iterations"=1000)
mcmc_pars <- c("iterations"=500000,"adaptive_frequency"=1000,"thin"=1000,
               "adaptive_iterations"=100000)
res <- serosolver(par_tab = par_tab, antibody_data = antibody_data, antigenic_map = antigenic_map, 
                              filename = filename, n_chains = 5, parallel = TRUE, 
                              mcmc_pars = mcmc_pars, start_inf_hist = NULL, verbose = TRUE)


load(paste0("/Users/luojialei/Desktop/BSP/Flu B/new/",filename,"/", filename, "_serosolver_settings.RData"))
res$all_diagnostics$p_thetas[[1]]
chains <- load_mcmc_chains(paste0("/Users/luojialei/Desktop/BSP/Flu B/new/", filename) ,par_tab,burnin=50000)
antibody_predictions <- plot_antibody_predictions(chains$theta_chain,chains$inf_chain,settings=res$settings)
print(antibody_predictions$proportion_correct)
antibody_predictions$p_hist_draws
antibody_predictions$p_pointrange

plot1 <- plot_model_fits(chains$theta_chain,chains$inf_chain,
                #individuals = indiv_start:(indiv_start + per_plot),
                individuals=1:9,
                known_infection_history=NULL, ## Set this to NULL for real data
                settings=res$settings,orientation="cross-sectional",expand_to_all_times = FALSE)
plot1 + facet_wrap(~individual, ncol=3)

chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no)) %>%
  filter(name %in% par_tab[par_tab$fixed == 0, "names"]) %>%
  rename(names=name,est=value) %>%
  left_join(par_tab %>% select(names,values)) %>%
  ggplot() + geom_density(aes(x=est,fill="Posterior"),alpha=0.5) + 
  geom_vline(aes(xintercept=values,linetype="True value"),col="black") +
  scale_color_viridis_d(name="") +
  scale_fill_viridis_d(name="") +
  scale_linetype_manual(name="",values=c("True value"="dashed")) +
  facet_wrap(~names,scales="free") +
  xlab("Value") +
  ylab("Density") +
  theme_minimal()
plot_attack_rates_pointrange(chains$inf_chain,settings = res$settings,
                             true_ar=NULL, ## Set this to NULL for real data
                             by_group=TRUE, 
                             plot_den = FALSE) +
  coord_cartesian(ylim=c(0,0.25)) + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(breaks=seq(8080,2022.5*4),labels=seq(8080,2022.5*4)/4)



