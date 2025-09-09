devtools::install_github("seroanalytics/serosolver")
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

filename <- "carlock_group_vl(1,1)"
setwd(paste0("/Users/luojialei/Desktop/BSP/Flu B/two_datasets/carlock/",filename))
resolution <- 1
serosolver::describe_priors()
prior_version <- 2

antigenic_map <- read_excel("/Users/luojialei/Desktop/BSP/Flu B/antigenic_maps/lastyear_ancvic_antigenic_map.xlsx")
antibody_data <- read_excel("/Users/luojialei/Desktop/BSP/Flu B/two_datasets/antibody_carlock_ancvic.xlsx")
antibody_data <- as.data.frame(antibody_data)
antigenic_map <- as.data.frame(antigenic_map)

strains_isolation_times <- unique(antigenic_map$inf_times)

par_tab <- read.csv("/Users/luojialei/Downloads/par_tab(in).csv")
par_tab <- as.data.frame(par_tab)
par_tab[par_tab$names == "max_measurement", "upper_bound"] <- 9

par_tab[par_tab$names == "infection_model_prior_shape1","values"] <- 1
par_tab[par_tab$names == "infection_model_prior_shape2","values"] <- 1

par_tab[par_tab$names=="infection_history_model_par1","stratification"] <-"Location"
par_tab[par_tab$names=="infection_history_model_par2","stratification"] <-"Location"

par_tab[par_tab$names=="infection_history_model_par1","stratification"] <-"Group"
par_tab[par_tab$names=="infection_history_model_par2","stratification"] <-"Group"

unique_indiv <- antibody_data[!duplicated(antibody_data$individual),]
ageMask <- create_age_mask(unique_indiv$birth, strains_isolation_times)

start_inf <- setup_infection_histories(antibody_data, strains_isolation_times)

mcmc_pars <- c("iterations"=5000,"adaptive_frequency"=100,"thin"=1,
               "adaptive_iterations"=1000)
mcmc_pars <- c("iterations"=500000,"adaptive_frequency"=1000,"thin"=1000,
               "adaptive_iterations"=100000)

ageMask <- create_age_mask(antibody_data %>% dplyr::select(individual,birth) %>% distinct() %>% pull(birth), strains_isolation_times)
sample_mask <- create_sample_mask(antibody_data,strains_isolation_times)
cbind(ageMask, sample_mask)

res <- serosolver(par_tab = par_tab, antibody_data = antibody_data, antigenic_map = antigenic_map,
                  filename = filename, n_chains = 5, parallel = TRUE,
                  mcmc_pars = mcmc_pars,
                  start_inf_hist = NULL,
                  verbose = TRUE)

load("/Users/luojialei/Desktop/BSP/Flu B/two_datasets/merged/merged_location_yl(1,1)/merged_location_yl(1,1)_serosolver_settings.RData")
#res$all_diagnostics$p_thetas[[1]]
par_tab <- serosolver_settings$par_tab
chains <- load_mcmc_chains("/Users/luojialei/Desktop/BSP/Flu B/two_datasets/merged/merged_location_yl(1,1)",par_tab,burnin=50000)
antibody_predictions <- plot_antibody_predictions(chains$theta_chain,chains$inf_chain, settings=serosolver_settings)
print(antibody_predictions$proportion_correct)
antibody_predictions$p_hist_draws
antibody_predictions$p_pointrange



p_ar <- plot_attack_rates(infection_histories=chains$inf_chain, antibody_data=serosolver_settings$antibody_data, 
                          demographics=NULL, par_tab=par_tab, possible_exposure_times=serosolver_settings$possible_exposure_times,
                          n_alive=NULL, pad_chain=TRUE, prior_pars=list(prior_version=serosolver_settings$prior_version,
                                                                        infection_model_prior_shape1=par_tab[par_tab$names=="infection_model_prior_shape1","values"],
                                                                        infection_model_prior_shape2=par_tab[par_tab$names=="infection_model_prior_shape2","values"]))

print(p_ar)

#n_alive <- get_n_alive(serosolver_settings$antibody_data, serosolver_settings$possible_exposure_times)

#ps_infhist <- plot_posteriors_infhist(inf_chain=chains$inf_chain, 
#                                      possible_exposure_times=serosolver_settings$possible_exposure_times, 
#                                      samples = 50,
#                                      n_alive=n_alive)

theta_chain <- as.data.frame(chains$theta_chain)
chain1 <- theta_chain[theta_chain$chain_no == 1,]
inf_chain <- chains$inf_chain
inf_chain1 <- inf_chain[inf_chain$chain_no == 1,]

antibody_preds <- get_antibody_level_predictions(chain = chain1, 
                                     infection_histories = inf_chain1, 
                                     antibody_data = serosolver_settings$antibody_data, 
                                     individuals = unique(serosolver_settings$antibody_data$individual),
                                     antigenic_map = serosolver_settings$antigenic_map, 
                                     nsamp = 50,
                                     par_tab = serosolver_settings$par_tab)

to_use <- antibody_preds$predictions

antibody_pred_p <- ggplot(to_use[to_use$individual %in% 1:9,])+
  geom_ribbon(aes(x=biomarker_id,ymin=lower, ymax=upper),fill="gray90")+
  geom_ribbon(aes(x=biomarker_id,ymin=lower_50, ymax=upper_50),fill="gray70")+
  geom_line(aes(x=biomarker_id, y=median))+
  geom_point(aes(x=biomarker_id, y=measurement))+
  coord_cartesian(ylim=c(0,8))+
  ylab("log titre") +
  xlab("Time of virus circulation") +
  theme_classic() +
  facet_wrap(~individual)
antibody_pred_p

infection_summary <- inf_chain[, .(total_x = sum(x)), by = .(i, chain_no, samp_no)]
setorder(infection_summary, i, chain_no)

compute_age <- function(df) {
  df <- df %>%
    group_by(individual) %>%
    mutate(age = max(sample_time, na.rm = TRUE) - birth) %>%
    ungroup()
  return(df)
}
age <- compute_age(serosolver_settings$antibody_data)
unique_individual_ages <- age %>%
  distinct(individual, age)
infection_summary <- infection_summary %>%
  left_join(unique_individual_ages, by = c("i" = "individual"))
infection_summary <- infection_summary %>%
  mutate(avg_infections_per_year_alive = total_x / age)
antibody_data <- serosolver_settings$antibody_data
latest_data <- infection_summary %>%
  group_by(i, chain_no) %>%
  filter(samp_no == max(samp_no))
#ggplot(latest_data, aes(x = avg_infections_per_year_alive, fill = factor(chain_no))) +
#  geom_density(alpha = 0.5) +
#  facet_wrap(~ chain_no, ncol = 1) +  # Separate plots for each chain
#  labs(
    #    title = "latest_data Distribution of avg_infections_per_year_alive Across All Individuals (5 chains)",
#    x = "Average Infections per Year Alive",
#    y = "Density"
#  ) +
#  theme_minimal()
ggplot(latest_data, aes(x = avg_infections_per_year_alive
                        #                        , fill = factor(chain_no)
)) +
  geom_density(alpha = 0.5,fill="gray") +
  #  facet_wrap(~ chain_no, ncol = 1) +  # Separate plots for each chain
  labs(
    #    title = "latest_data Distribution of avg_infections_per_year_alive Across All Individuals (5 chains combined)",
    x = "Average Infections per Year Alive",
    y = "Density"
  ) +
  theme_minimal()
setDT(latest_data)
latest_stats <- latest_data[, .(
  mean_infections = mean(total_x, na.rm = TRUE),
  lower_cri = quantile(total_x, 0.025, na.rm = TRUE),
  upper_cri = quantile(total_x, 0.975, na.rm = TRUE)
), by = .(i, age)][order(age)]
ggplot(latest_data, aes(x = total_x)) +
  geom_histogram(binwidth = 1,alpha = 1, fill = "gray", color = "black") +
  #  coord_flip()
  #  facet_wrap(~ chain_no, ncol = 1) +  # Separate plots for each chain
  labs(
    #    title = "total infection frequency",
    x = "Total Infections",
    y = "Frequency"
  ) +
  theme_minimal()
latest_data[, age_group := cut(age, 
                               breaks = c(0, 5, 10, 20, 30, 40, 50, 60, 70, 80, Inf),
                               labels = c("<5", "(5,10]", "(10,20]", "(20,30]", "(30,40]", "(40,50]", 
                                          "(50,60]", "(60,70]", "(70,80]", "80+"),
                               right = TRUE)]
latest_data[, infections_per_decade := total_x / (age / 10)]
ggplot(latest_data, aes(x = age_group, y = infections_per_decade)) +
  geom_boxplot(fill = "gray", color = "black") +
  labs(
    #    title = "Posterior Median Infections per 10-Year Period",
    #    subtitle = "Stratified by Age Group at Time of Infection",
    x = "Age Group at Time of Infection",
    y = "Infections per Decade"
  ) +
  theme_minimal()


plot1 <- plot_model_fits(chains$theta_chain,chains$inf_chain,
                #individuals = indiv_start:(indiv_start + per_plot),
                individuals=1:9,
                known_infection_history=NULL, ## Set this to NULL for real data
                settings=serosolver_settings,orientation="cross-sectional",expand_to_all_times = FALSE)
plot1 + facet_wrap(~individual, ncol=3)

chains$theta_chain %>% pivot_longer(-c(samp_no,chain_no)) %>%
  filter(name %in% serosolver_settings$par_tab[serosolver_settings$par_tab$fixed == 0, "names"]) %>%
  rename(names=name,est=value) %>%
  left_join(serosolver_settings$par_tab %>% select(names,values)) %>%
  ggplot() + geom_density(aes(x=est,fill="Posterior"),alpha=0.5) + 
  geom_vline(aes(xintercept=values,linetype="True value"),col="black") +
  scale_color_viridis_d(name="") +
  scale_fill_viridis_d(name="") +
  scale_linetype_manual(name="",values=c("True value"="dashed")) +
  facet_wrap(~names,scales="free") +
  xlab("Value") +
  ylab("Density") +
  theme_minimal()
#plot_attack_rates(chains$inf_chain,settings = serosolver_settings,
#                             true_ar=NULL, ## Set this to NULL for real data
#                             by_group=TRUE, 
#                             plot_den = FALSE) +
#  coord_cartesian(ylim=c(0,0.25)) + scale_y_continuous(expand=c(0,0)) + scale_x_continuous(breaks=seq(8080,2022.5*4),labels=seq(8080,2022.5*4)/4)

#plot_posteriors_infhist(chains$inf_chain, serosolver_settings$possible_exposure_times, n_alive)
