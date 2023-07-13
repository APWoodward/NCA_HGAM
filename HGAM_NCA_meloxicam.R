# Andrew P. Woodward, 2022-23.
# Noncompartmental analysis of the meloxicam dataset (https://doi.org/10.1007/bf01059644) as described in Warner et al. 2020 (https://doi.org/10.3389/fvets.2020.00548).
# Here we utilize a hierarchical generalized additive model for the change in meloxicam concentrations with time, with separate population splines for stage and route.
#     The hierarchical GAM are expressed following Pedersen et al. 2019 (https://doi.org/10.7717/peerj.6876).
#     Spline-based NCA was previously reported in (https://doi.org/10.1002/pst.336 and https://doi.org/10.1177/1471082X17706018); in this work we use the popular 'mgcv' implementation of regression splines via 'brms'.

# Load the required packages.
library(reshape2)
library(brms)
library(ggplot2)
library(viridis)
library(ncappc)
library(emmeans)
library(boot)
library(ggpubr)
library(future.apply)
library(rlist)
library(ggdist)

# The time will be transformed as log(time+1); the time-concentration function is presumed to be multiplicative.
#     Note that without the +1 adjustment the time 0 observations will be dropped.
#     For terminal phase observations, the log-time scale prevents the uncertainty from growing in-between the observations ('ballooning' effect), as the observations are similarly spaced on log-time scale.
#     The authors (https://doi.org/10.3389/fvets.2020.00548) for this data note the lower limit of quantification as 0.002.
meloxicam_IV_long <- melt(meloxicam_cow_data_IV, id.vars = c('ID','STAGE','ROUTE'))
meloxicam_PO_long <- melt(meloxicam_cow_data_PO, id.vars = c('ID','STAGE','ROUTE'))
meloxicam_cow_data <- rbind(meloxicam_IV_long, meloxicam_PO_long)
colnames(meloxicam_cow_data)[4] <- 'TIME'
colnames(meloxicam_cow_data)[5] <- 'CONC'
meloxicam_cow_data$CENS <- 0
meloxicam_cow_data$CENS[meloxicam_cow_data$CONC == 0.002] <- -1
meloxicam_cow_data$TIME <- as.character(meloxicam_cow_data$TIME)
meloxicam_cow_data$TIME <- as.numeric(meloxicam_cow_data$TIME)
meloxicam_cow_data$ID <- factor(meloxicam_cow_data$ID)
meloxicam_cow_data$TIME_log1 <- log(meloxicam_cow_data$TIME+1)
meloxicam_cow_data$TIME_hr <- meloxicam_cow_data$TIME/60
meloxicam_cow_data$TIME_hr_log1 <- log(meloxicam_cow_data$TIME_hr+1)

# Generalized additive model for the meloxicam data, using a covariate-dependent smoother, and subject-level smoothers with shared wiggliness (https://peerj.com/articles/6876/).
# Weak priors are set on the categorical terms ('intercepts') but the covariates on the smooths are flat.
# Then generate predicted concentrations and their uncertainty for visualization.
meloxicam_priors <- set_prior('normal(0,0.5)', class = 'sigma') + set_prior('normal(0,5)', coef = 'Intercept') + set_prior('normal(0,5)', coef = 'ROUTEPO') + set_prior('normal(0,5)', coef = 'ROUTEPO:STAGEPostMpartum') + set_prior('normal(0,5)', coef = 'STAGEPostMpartum')
meloxicam_dis_mod_m2 <- brm(bf(log(CONC)|cens(CENS) ~  1 + ROUTE * STAGE + s(TIME_hr_log1, by = interaction(ROUTE, STAGE)) + s(TIME_hr_log1, ID, bs = 'fs', m = 2), center = FALSE),  data = meloxicam_cow_data, cores = 4, chains = 4, iter = 2000, control = list(adapt_delta = 0.99, max_treedepth = 12), prior = meloxicam_priors)
summary(meloxicam_dis_mod_m2)
meloxicam_residuals <- data.frame(TIME_hr = meloxicam_cow_data$TIME_hr[meloxicam_cow_data$CENS == 0], CONC = meloxicam_cow_data$CONC[meloxicam_cow_data$CENS == 0], RESIDUAL = residuals(meloxicam_dis_mod_m2)[meloxicam_cow_data$CENS == 0,1], FITTED = fitted(meloxicam_dis_mod_m2)[meloxicam_cow_data$CENS == 0,1])
meloxicam_residual_time <- ggplot(data = meloxicam_residuals, aes(x = TIME_hr, y = RESIDUAL)) + geom_point(alpha = 0.3) + theme_bw() + ylim(c(-2,2)) + geom_hline(yintercept = 0, alpha = 0.4) + ylab('Residual (log scale)') + xlab('Time since administration (h)')
meloxicam_residual_pred <- ggplot(data = meloxicam_residuals, aes(x = FITTED, y = RESIDUAL)) + geom_point(alpha = 0.3) + theme_bw() + ylim(c(-2,2)) + geom_hline(yintercept = 0, alpha = 0.4)  + ylab('Residual (log scale)') + xlab('Posterior mean predicted concentration (log scale)')
ggarrange(meloxicam_residual_time, meloxicam_residual_pred, ncol = 1)

# Generate population predictions (global smoother) by covariate, ignoring the subject-level effects, and visualize.
new_meloxicam_data_all <- expand.grid(TIME_hr_log1 = seq(0,log(180),0.1), ROUTE = c('IV','PO'), STAGE = c('Mid-lactation','Post-partum'), ID = NA)
new_meloxicam_pred_all <- cbind(new_meloxicam_data_all, fitted(meloxicam_dis_mod_m2, newdata = new_meloxicam_data_all, re_formula = NA, probs = c(0.05,0.25,0.75,0.95)))
new_meloxicam_pred_all$TIME_hr <- exp(new_meloxicam_pred_all$TIME_hr_log1)-1
meloxicam_pop_plot <- ggplot(data = new_meloxicam_pred_all, aes(x = TIME_hr, y = exp(Estimate))) + geom_line() + geom_ribbon(aes(ymin = exp(Q5), ymax = exp(Q95)), alpha = 0.2) + geom_ribbon(aes(ymin = exp(Q25), ymax = exp(Q75)), alpha = 0.3) + facet_grid(cols = vars(STAGE), rows = vars(ROUTE)) + theme_bw() + theme(legend.position = 'none') + coord_cartesian(ylim = c(0.001,10), xlim = c(0,156)) + geom_point(data = meloxicam_cow_data, aes(x = TIME_hr, y = CONC, shape = factor(CENS)), alpha = 0.2, inherit.aes = FALSE, fill = 'white') + geom_line(data = meloxicam_cow_data, aes(x = TIME_hr, y = CONC, group = ID), alpha = 0.1) + scale_shape_manual(values = c(21,19)) + ylab('Concentration (mg/L)') + xlab('Time since administration (hours)') + scale_x_continuous(breaks = c(0,24,48,72,96,120,144), labels = c('0','24','48','72','96','120','144')) + scale_y_log10(breaks = c(0.001,0.01,0.1,1,10), labels = c(0.001,0.01,0.1,1,10))+ theme(panel.grid.minor = element_blank())
ggsave(meloxicam_pop_plot, file = 'meloxicam_pop_plot.svg', width = 200, height = 200, units = 'mm')

# Define a function that returns the predicted concentration for a specific time, subject and covariates, which will be used during integration.
single_AUC_int <- function(time_int, route_in, stage_in, subject, sample){
  time_int_data <- data.frame(TIME_hr_log1 = log(time_int+1), ID = subject, ROUTE = route_in, STAGE = stage_in)
  if(is.na(subject)){
    if(is.na(sample)){
      concentration_estimate <- (exp(fitted(meloxicam_dis_mod_m2, newdata = time_int_data, allow_new_levels = TRUE, re_formula = NA)[,1]))
    }
    else{
      concentration_estimate <- (exp(fitted(meloxicam_dis_mod_m2, newdata = time_int_data, allow_new_levels = TRUE, re_formula = NA, draw_ids = sample)[,1]))  
    }
  }
  else{
    if(is.na(sample)){
      concentration_estimate <- (exp(fitted(meloxicam_dis_mod_m2, newdata = time_int_data)[,1]))
    }
    else{
      concentration_estimate <- (exp(fitted(meloxicam_dis_mod_m2, newdata = time_int_data, draw_ids = sample)[,1]))   
    }
  }
  return(concentration_estimate)
}

# Generate population predictions of the AUC, Tmax and Cmax, by values of the covariates, using the HGAM, and conventional NCA.
meloxicam_pop_AUC <- expand.grid(ROUTE = c('IV','PO'), STAGE = c('Mid-lactation','Post-partum'))
meloxicam_pop_AUC$AUC_GAM <- 0
for (i in 1:length(meloxicam_pop_AUC$AUC_GAM)){
  meloxicam_pop_AUC$AUC_GAM[i] <- as.numeric(integrate(single_AUC_int,0,144,route_in = meloxicam_pop_AUC$ROUTE[i], stage_in = meloxicam_pop_AUC$STAGE[i], subject = NA, sample = NA)[1])
}
meloxicam_sub_AUC <- data.frame(ID = numeric(length(unique(meloxicam_cow_data$ID))), ROUTE = numeric(length(unique(meloxicam_cow_data$ID))), STAGE = numeric(length(unique(meloxicam_cow_data$ID))))
meloxicam_sub_AUC$AUCinf_trap <- 0
for (i in 1:length(unique(meloxicam_cow_data$ID))){
  subject_ind <- unique(meloxicam_cow_data$ID)[i]
  meloxicam_sub_AUC$ID[i] <- subject_ind
  meloxicam_sub_AUC$ROUTE[subject_ind] <- unique(meloxicam_cow_data$ROUTE[meloxicam_cow_data$ID == subject_ind])
  meloxicam_sub_AUC$STAGE[subject_ind] <- unique(meloxicam_cow_data$STAGE[meloxicam_cow_data$ID == subject_ind])
  meloxicam_sub_AUC$AUC_GAM[i] <- as.numeric(integrate(single_AUC_int,0,144,route_in = meloxicam_sub_AUC$ROUTE[i], stage_in = meloxicam_sub_AUC$STAGE[i], subject = subject_ind, sample = NA)[1])
  meloxicam_sub_AUC$AUC_trap[i] <- (ncappc(obsFile = data.frame(ID = meloxicam_cow_data$ID[((meloxicam_cow_data$ID == subject_ind) & (meloxicam_cow_data$CENS == 0))], TIME = meloxicam_cow_data$TIME_hr[(meloxicam_cow_data$ID == subject_ind & (meloxicam_cow_data$CENS == 0))], DV = meloxicam_cow_data$CONC[(meloxicam_cow_data$ID == subject_ind & (meloxicam_cow_data$CENS == 0))]), onlyNCA = TRUE, printOut = FALSE, extrapolate = TRUE, evid = FALSE))$ncaOutput$AUClast
  meloxicam_sub_AUC$AUCinf_trap[i] <- (ncappc(obsFile = data.frame(ID = meloxicam_cow_data$ID[((meloxicam_cow_data$ID == subject_ind) & (meloxicam_cow_data$CENS == 0))], TIME = meloxicam_cow_data$TIME_hr[(meloxicam_cow_data$ID == subject_ind & (meloxicam_cow_data$CENS == 0))], DV = meloxicam_cow_data$CONC[(meloxicam_cow_data$ID == subject_ind & (meloxicam_cow_data$CENS == 0))]), onlyNCA = TRUE, printOut = FALSE, extrapolate = TRUE, evid = FALSE))$ncaOutput$AUCINF_obs
  meloxicam_sub_AUC$Cmax_trap[i] <- (ncappc(obsFile = data.frame(ID = meloxicam_cow_data$ID[((meloxicam_cow_data$ID == subject_ind) & (meloxicam_cow_data$CENS == 0))], TIME = meloxicam_cow_data$TIME_hr[(meloxicam_cow_data$ID == subject_ind & (meloxicam_cow_data$CENS == 0))], DV = meloxicam_cow_data$CONC[(meloxicam_cow_data$ID == subject_ind & (meloxicam_cow_data$CENS == 0))]), onlyNCA = TRUE, printOut = FALSE, extrapolate = TRUE, evid = FALSE))$ncaOutput$Cmax
  meloxicam_sub_AUC$Tmax_trap[i] <- (ncappc(obsFile = data.frame(ID = meloxicam_cow_data$ID[((meloxicam_cow_data$ID == subject_ind) & (meloxicam_cow_data$CENS == 0))], TIME = meloxicam_cow_data$TIME_hr[(meloxicam_cow_data$ID == subject_ind & (meloxicam_cow_data$CENS == 0))], DV = meloxicam_cow_data$CONC[(meloxicam_cow_data$ID == subject_ind & (meloxicam_cow_data$CENS == 0))]), onlyNCA = TRUE, printOut = FALSE, extrapolate = TRUE, evid = FALSE))$ncaOutput$Tmax
  meloxicam_max_int <- optimize(single_AUC_int,c(0,144), maximum = TRUE, route_in = meloxicam_sub_AUC$ROUTE[i], stage_in = meloxicam_sub_AUC$STAGE[i], subject = subject_ind, sample = NA)
  meloxicam_sub_AUC$Cmax_GAM[i] <- as.numeric(meloxicam_max_int$objective)
  meloxicam_sub_AUC$Tmax_GAM[i] <- (meloxicam_max_int$maximum)
}
meloxicam_sub_AUC$DOSE[meloxicam_sub_AUC$ROUTE == 'IV'] <- 0.2
meloxicam_sub_AUC$DOSE[meloxicam_sub_AUC$ROUTE == 'PO'] <- 1

# Generate population predictions across the covariates.
meloxicam_pop_AUC$AUC_trap[1] <- exp(mean(log(meloxicam_sub_AUC$AUC_trap[((meloxicam_sub_AUC$ROUTE == meloxicam_pop_AUC$ROUTE[1]) & (meloxicam_sub_AUC$STAGE == meloxicam_pop_AUC$STAGE[1]))])))
meloxicam_pop_AUC$AUC_trap[2] <- exp(mean(log(meloxicam_sub_AUC$AUC_trap[((meloxicam_sub_AUC$ROUTE == meloxicam_pop_AUC$ROUTE[2]) & (meloxicam_sub_AUC$STAGE == meloxicam_pop_AUC$STAGE[2]))])))
meloxicam_pop_AUC$AUC_trap[3] <- exp(mean(log(meloxicam_sub_AUC$AUC_trap[((meloxicam_sub_AUC$ROUTE == meloxicam_pop_AUC$ROUTE[3]) & (meloxicam_sub_AUC$STAGE == meloxicam_pop_AUC$STAGE[3]))])))
meloxicam_pop_AUC$AUC_trap[4] <- exp(mean(log(meloxicam_sub_AUC$AUC_trap[((meloxicam_sub_AUC$ROUTE == meloxicam_pop_AUC$ROUTE[4]) & (meloxicam_sub_AUC$STAGE == meloxicam_pop_AUC$STAGE[4]))])))


# Generate the subject-level predictions, using the subject-level covariates, and visualize the results.
new_meloxicam_data_sub <- expand.grid(TIME_hr_log1 = seq(0,log(180),0.1), ID = unique(meloxicam_cow_data$ID))
new_meloxicam_data_sub$ROUTE <- 0
new_meloxicam_data_sub$STAGE <- 0
for(i in 1:length(unique(meloxicam_cow_data$ID))){
  data_sub_index <- which(new_meloxicam_data_sub$ID == (unique(meloxicam_cow_data$ID))[i])
  new_meloxicam_data_sub$ROUTE[data_sub_index] <- unique(meloxicam_cow_data$ROUTE[meloxicam_cow_data$ID == ((unique(meloxicam_cow_data$ID))[i])])
  new_meloxicam_data_sub$STAGE[data_sub_index] <- unique(meloxicam_cow_data$STAGE[meloxicam_cow_data$ID == ((unique(meloxicam_cow_data$ID))[i])])
}
new_meloxicam_pred_sub <- cbind(new_meloxicam_data_sub, fitted(meloxicam_dis_mod_m2, newdata = new_meloxicam_data_sub, probs = c(0.05,0.25,0.75,0.95)))
new_meloxicam_pred_sub$TIME_hr <- exp(new_meloxicam_pred_sub$TIME_hr_log1)-1
meloxicam_sub_plot <- ggplot(data = new_meloxicam_pred_sub, aes(x = TIME_hr, y = exp(Estimate))) + geom_line(aes(color = interaction(ROUTE,STAGE))) + geom_ribbon(aes(fill = interaction(ROUTE,STAGE), ymin = exp(Q25), ymax = exp(Q75)), alpha = 0.4) + geom_ribbon(aes(fill = interaction(ROUTE,STAGE), ymin = exp(Q5), ymax = exp(Q95)), alpha = 0.2) + theme_bw() + scale_y_log10(breaks = c(0.001,0.01,0.1,1,10), labels = c(0.001,0.01,0.1,1,10)) + facet_wrap(~ID, ncol = 3) + geom_point(data = meloxicam_cow_data, aes(x = TIME_hr, y = CONC, color = interaction(ROUTE,STAGE), shape = factor(CENS)), inherit.aes = FALSE, fill = 'white', alpha = 0.6) + coord_cartesian(ylim = c(0.001,10), xlim = c(0,150)) + scale_shape_manual(values = c(21,19), guide = 'none') + ylab('Concentration (mg/L)') + xlab('Time since administration (hours)') + scale_x_continuous(breaks = c(0,24,48,72,96,120,144), labels = c('0','24','48','72','96','120','144')) + theme(legend.title = element_blank(), legend.position = 'top', panel.grid.minor = element_blank()) + scale_color_viridis(discrete = TRUE, end = 0.8) + scale_fill_viridis(discrete = TRUE, end = 0.8)
ggsave(meloxicam_sub_plot, file = 'meloxicam_sub_plot.svg', width = 200, height = 200, units = 'mm')

# Define a bootstrap function for AUC.
boot_fun_AUC_raw <- function(input_data, index){
  boot_iter_data <- input_data[index,]
  AUC_trap_mod  <- lm(log(AUC_trap)~1+ROUTE*STAGE, data = boot_iter_data, offset = log(DOSE))
  # The order of outputs here are IV:ML, PO:ML, IV:PP, PO:PP; reported on log-scale.
  boot_iter_est <- predict(AUC_trap_mod, newdata = expand.grid(ROUTE = c('IV','PO'), STAGE = c('Mid-lactation','Post-partum'), DOSE = 1))
  # The order of outputs here are IV:PP/IV:ML, PO:PP/PO:ML, PO:ML/IV:ML, PO:PP/IV:PP; reported on log-scale.
  boot_iter_par <- c((boot_iter_est[3]-boot_iter_est[1]), (boot_iter_est[4]-boot_iter_est[2]), (boot_iter_est[2]-boot_iter_est[1]), (boot_iter_est[4]-boot_iter_est[3]))
  return(c(boot_iter_est, boot_iter_par))
}

# Generate bootstrap confidence statements for the AUC and AUC ratios.
meloxicam_trap_boot <- boot(data = meloxicam_sub_AUC, boot_fun_AUC_raw, 1000, strata = c(rep(1,5),rep(2,6),rep(3,6),rep(4,6)))
meloxicam_trap_conf <- data.frame(parameter = character(8), L90CI = numeric(8), L50CI = numeric(8), point = numeric(8), U50CI = numeric(8), U90CI = numeric(8))
meloxicam_trap_conf$parameter <- c('IV:ML','PO:ML','IV:PP','PO:PP','IV:PP/IV:ML','PO:PP/PO:ML','PO:ML/IV:ML','PO:PP/IV:PP')
for(i in 1:8){
  iter_boot <- boot(data = meloxicam_sub_AUC, boot_fun_AUC_raw, 1000, strata = c(rep(1,5),rep(2,6),rep(3,6),rep(4,6)))
  iter_boot_ci <- (boot.ci(iter_boot, conf = 0.9, type = 'perc', index = i))
  meloxicam_trap_conf[i,4] <- exp(iter_boot_ci$t0)
  meloxicam_trap_conf[i,2] <- exp(iter_boot_ci$perc[4])
  meloxicam_trap_conf[i,6] <- exp(iter_boot_ci$perc[5])
  iter_boot_ci <- (boot.ci(iter_boot, conf = 0.5, type = 'perc', index = i))
  meloxicam_trap_conf[i,3] <- exp(iter_boot_ci$perc[4])
  meloxicam_trap_conf[i,5] <- exp(iter_boot_ci$perc[5])
}
meloxicam_trap_conf

# Generate posterior samples for the population AUC (global smoothers ignoring the subject terms) by covariates.
#     Then generate posterior distributions for the ratios by condition (contrasts), and visualize them.
meloxicam_pop_cov <- expand.grid(unique(STAGE = meloxicam_cow_data$STAGE, ROUTE = meloxicam_cow_data$ROUTE))
plan(multisession, workers = 4)
AUC_meloxicam_list <- future_lapply((1:length(meloxicam_pop_cov$STAGE)), function(condition_ind){
  AUC_est_meloxicam_cov <- ndraws(meloxicam_dis_mod_m2)
  for (h in 1:(ndraws(meloxicam_dis_mod_m2))){
    AUC_est_meloxicam_cov[h] <- as.numeric(integrate(single_AUC_int, 0, 144, rel.tol = 0.001, abs.tol = 0.001, sample = h, subject = NA, route_in = meloxicam_pop_cov$ROUTE[condition_ind], stage_in = meloxicam_pop_cov$STAGE[condition_ind])[1])
  }
  AUC_est_meloxicam_cov
})
meloxicam_AUC_samples <- as.data.frame(c(AUC_meloxicam_list[[1]],AUC_meloxicam_list[[2]],AUC_meloxicam_list[[3]],AUC_meloxicam_list[[4]]))
colnames(meloxicam_AUC_samples)[1] <- 'AUC'
meloxicam_AUC_samples$STAGE <- c(rep(meloxicam_pop_cov$STAGE[1],ndraws(meloxicam_dis_mod_m2)),rep(meloxicam_pop_cov$STAGE[2],ndraws(meloxicam_dis_mod_m2)),rep(meloxicam_pop_cov$STAGE[3],ndraws(meloxicam_dis_mod_m2)),rep(meloxicam_pop_cov$STAGE[4],ndraws(meloxicam_dis_mod_m2)))
meloxicam_AUC_samples$ROUTE <- c(rep(meloxicam_pop_cov$ROUTE[1],ndraws(meloxicam_dis_mod_m2)),rep(meloxicam_pop_cov$ROUTE[2],ndraws(meloxicam_dis_mod_m2)),rep(meloxicam_pop_cov$ROUTE[3],ndraws(meloxicam_dis_mod_m2)),rep(meloxicam_pop_cov$ROUTE[4],ndraws(meloxicam_dis_mod_m2)))
meloxicam_AUC_samples$AUC[meloxicam_AUC_samples$ROUTE == 'IV'] <- (meloxicam_AUC_samples$AUC[meloxicam_AUC_samples$ROUTE == 'IV'])/0.2
posterior_AUCcov_plot <- ggplot(data = meloxicam_AUC_samples, aes(x = AUC, y = ROUTE:STAGE)) + stat_dotsinterval(.width = c(0.5,0.9), quantiles = 100, interval_alpha = 0.8, point_alpha = 0.8) + theme_bw() + coord_cartesian(xlim = c(0,200)) + xlab(expression(Dose-corrected~AUC[last])) + ylab('') + scale_y_discrete(limits = rev)
meloxicam_AUC_ratio     <- data.frame(AUC_ratio = c(meloxicam_AUC_samples$AUC[((meloxicam_AUC_samples$ROUTE == 'IV') & (meloxicam_AUC_samples$STAGE == 'Post-partum'))]/  meloxicam_AUC_samples$AUC[((meloxicam_AUC_samples$ROUTE == 'IV') & (meloxicam_AUC_samples$STAGE == 'Mid-lactation'))],
                                                    meloxicam_AUC_samples$AUC[((meloxicam_AUC_samples$ROUTE == 'PO') & (meloxicam_AUC_samples$STAGE == 'Post-partum'))]/  meloxicam_AUC_samples$AUC[((meloxicam_AUC_samples$ROUTE == 'PO') & (meloxicam_AUC_samples$STAGE == 'Mid-lactation'))],
                                                    meloxicam_AUC_samples$AUC[((meloxicam_AUC_samples$ROUTE == 'PO') & (meloxicam_AUC_samples$STAGE == 'Mid-lactation'))]/meloxicam_AUC_samples$AUC[((meloxicam_AUC_samples$ROUTE == 'IV') & (meloxicam_AUC_samples$STAGE == 'Mid-lactation'))],
                                                    meloxicam_AUC_samples$AUC[((meloxicam_AUC_samples$ROUTE == 'PO') & (meloxicam_AUC_samples$STAGE == 'Post-partum'))]/  meloxicam_AUC_samples$AUC[((meloxicam_AUC_samples$ROUTE == 'IV') & (meloxicam_AUC_samples$STAGE == 'Post-partum'))]))
meloxicam_AUC_ratio$comparison <- c(rep('IV:PP/IV:ML', ndraws(meloxicam_dis_mod_m2)),rep('PO:PP/PO:ML', ndraws(meloxicam_dis_mod_m2)),rep('PO:ML/IV:ML', ndraws(meloxicam_dis_mod_m2)),rep('PO:PP/IV:PP', ndraws(meloxicam_dis_mod_m2)))                                  
pop_meloxicam_table <- data.frame(parameter = character(length = 8), L90CI = numeric(length = 8), L50CI = numeric(length = 8), point = numeric(length = 8), U50CI = numeric(length = 8), U90CI = numeric(length = 8))
pop_meloxicam_table$parameter <- meloxicam_trap_conf$parameter
pop_meloxicam_table[(pop_meloxicam_table$parameter == 'IV:ML'),2:6] <- quantile((meloxicam_AUC_samples$AUC[(meloxicam_AUC_samples$STAGE == 'Mid-lactation') & (meloxicam_AUC_samples$ROUTE == 'IV')]), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
pop_meloxicam_table[(pop_meloxicam_table$parameter == 'PO:ML'),2:6] <- quantile((meloxicam_AUC_samples$AUC[(meloxicam_AUC_samples$STAGE == 'Mid-lactation') & (meloxicam_AUC_samples$ROUTE == 'PO')]), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
pop_meloxicam_table[(pop_meloxicam_table$parameter == 'IV:PP'),2:6] <- quantile((meloxicam_AUC_samples$AUC[(meloxicam_AUC_samples$STAGE == 'Post-partum') & (meloxicam_AUC_samples$ROUTE == 'IV')]), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
pop_meloxicam_table[(pop_meloxicam_table$parameter == 'PO:PP'),2:6] <- quantile((meloxicam_AUC_samples$AUC[(meloxicam_AUC_samples$STAGE == 'Post-partum') & (meloxicam_AUC_samples$ROUTE == 'PO')]), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
pop_meloxicam_table[(pop_meloxicam_table$parameter == 'IV:PP/IV:ML'),2:6] <- quantile((meloxicam_AUC_ratio$AUC_ratio[meloxicam_AUC_ratio$comparison == 'IV:PP/IV:ML']), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
pop_meloxicam_table[(pop_meloxicam_table$parameter == 'PO:PP/PO:ML'),2:6] <- quantile((meloxicam_AUC_ratio$AUC_ratio[meloxicam_AUC_ratio$comparison == 'PO:PP/PO:ML']), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
pop_meloxicam_table[(pop_meloxicam_table$parameter == 'PO:ML/IV:ML'),2:6] <- quantile((meloxicam_AUC_ratio$AUC_ratio[meloxicam_AUC_ratio$comparison == 'PO:ML/IV:ML']), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
pop_meloxicam_table[(pop_meloxicam_table$parameter == 'PO:PP/IV:PP'),2:6] <- quantile((meloxicam_AUC_ratio$AUC_ratio[meloxicam_AUC_ratio$comparison == 'PO:PP/IV:PP']), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
posterior_AUCratio_plot <- ggplot(data = meloxicam_AUC_ratio, aes(x = AUC_ratio, y = comparison)) + stat_dotsinterval(.width = c(0.5,0.9), quantiles = 100, interval_alpha = 0.8, point_alpha = 0.8) + theme_bw() + scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10)) + coord_cartesian(xlim = c(0.1,10)) + theme(panel.grid.minor = element_blank()) + geom_vline(xintercept = 1) + xlab(expression(Dose-corrected~AUC[last]~ratio)) + ylab('') + scale_y_discrete(limits = rev)
meloxicam_pop_plots <- ggarrange(posterior_AUCcov_plot,posterior_AUCratio_plot, ncol = 1)
ggsave(meloxicam_pop_plots, file = 'meloxicam_pop_plots.svg', units = 'mm', width = 200, height = 200)

# Generate posterior samples for the subject AUC, including the covariates (fully conditional on the subject parameters).
#     Then generate posterior distributions for the ratios by condition (contrasts), and visualize them.
#     Compared to the estimates for an average subject, these are for the average of the observed subjects.
#     Also visualize the individual subject AUC estimates and compare them to the trapezoidal estimates.
meloxicam_sub_cov <- meloxicam_sub_AUC[,1:3]
plan(multisession, workers = 12)
sub_meloxicam_list <- future_lapply((1:length(meloxicam_sub_cov$ID)), function(subject_ind){
  AUC_est_meloxicam_sub <-  numeric(length = ndraws(meloxicam_dis_mod_m2))
  for (h in 1:(ndraws(meloxicam_dis_mod_m2))){
    AUC_est_meloxicam_sub[h] <- as.numeric(integrate(single_AUC_int, 0, 144, rel.tol = 0.001, abs.tol = 0.001, sample = h, subject = meloxicam_sub_cov$ID[subject_ind], route_in = meloxicam_sub_cov$ROUTE[subject_ind], stage_in = meloxicam_sub_cov$STAGE[subject_ind])[1])
  }
  AUC_est_meloxicam_sub
})
sub_meloxicam_pred <- lapply((1:length(meloxicam_sub_cov$ID)), function(subject_ind){
  sub_meloxicam_iter <- data.frame(AUC = sub_meloxicam_list[[subject_ind]], ID = meloxicam_sub_cov$ID[subject_ind], ROUTE = meloxicam_sub_cov$ROUTE[subject_ind], STAGE = meloxicam_sub_cov$STAGE[subject_ind])
  sub_meloxicam_iter
})
sub_meloxicam_pred <- list.rbind(sub_meloxicam_pred)
sub_meloxicam_data <- sub_meloxicam_pred
sub_meloxicam_data$AUC[sub_meloxicam_data$ROUTE == 'IV'] <- sub_meloxicam_data$AUC[sub_meloxicam_data$ROUTE == 'IV']/0.2
meloxicam_AUC_plot <- ggplot(sub_meloxicam_data, aes(x = AUC, fill = interaction(ROUTE,STAGE))) + stat_histinterval(alpha = 0.7, size = 2) + facet_wrap(~ID) + theme_bw() + scale_fill_viridis(discrete = TRUE, end = 0.8) + theme(legend.position = 'top', legend.title = element_blank()) + ylab('') + scale_y_continuous(breaks = NULL) + coord_cartesian(xlim = c(0,200)) + geom_vline(data = meloxicam_sub_AUC, aes(xintercept = AUC_dose), alpha = 0.6, linetype = 'dashed') + xlab(expression(Dose-corrected~AUC[last]~(mg.hr/L)))
ggsave(meloxicam_AUC_plot, file = 'meloxicam_AUC_plot.svg', width = 200, height = 200, units = 'mm') 
emm_meloxicam_samp <- cbind(meloxicam_sub_cov, as.data.frame(list.rbind(sub_meloxicam_list)))
emm_meloxicam_pred <- as.data.frame(matrix(nrow = 4, ncol = (dim(emm_meloxicam_samp)[2])-1, 0))
emm_meloxicam_pred[,1:2] <- meloxicam_pop_cov
colnames(emm_meloxicam_pred)[1:2] <- colnames(meloxicam_pop_cov)[1:2]
for(i in 1:(dim(emm_meloxicam_pred)[2]-dim(meloxicam_sub_cov)[2])){
  emm_meloxicam_pred[1,i+(dim(meloxicam_sub_cov)[2]-1)] <- exp(mean(log(emm_meloxicam_samp[((emm_meloxicam_samp$ROUTE == meloxicam_pop_cov$ROUTE[1]) & (emm_meloxicam_samp$STAGE == meloxicam_pop_cov$STAGE[1])),i+(dim(meloxicam_sub_cov)[2])]/0.2)))
  emm_meloxicam_pred[2,i+(dim(meloxicam_sub_cov)[2]-1)] <- exp(mean(log(emm_meloxicam_samp[((emm_meloxicam_samp$ROUTE == meloxicam_pop_cov$ROUTE[2]) & (emm_meloxicam_samp$STAGE == meloxicam_pop_cov$STAGE[2])),i+(dim(meloxicam_sub_cov)[2])]/0.2)))
  emm_meloxicam_pred[3,i+(dim(meloxicam_sub_cov)[2]-1)] <- exp(mean(log(emm_meloxicam_samp[((emm_meloxicam_samp$ROUTE == meloxicam_pop_cov$ROUTE[3]) & (emm_meloxicam_samp$STAGE == meloxicam_pop_cov$STAGE[3])),i+(dim(meloxicam_sub_cov)[2])])))
  emm_meloxicam_pred[4,i+(dim(meloxicam_sub_cov)[2]-1)] <- exp(mean(log(emm_meloxicam_samp[((emm_meloxicam_samp$ROUTE == meloxicam_pop_cov$ROUTE[4]) & (emm_meloxicam_samp$STAGE == meloxicam_pop_cov$STAGE[4])),i+(dim(meloxicam_sub_cov)[2])])))
}
emm_meloxicam_ratio <- data.frame(AUC_ratio = as.numeric(t(emm_meloxicam_pred[2,3:4002]/emm_meloxicam_pred[1,3:4002])), comparison = 'IV:PP/IV:ML')
emm_meloxicam_ratio <- rbind(emm_meloxicam_ratio, data.frame(AUC_ratio = as.numeric(t(emm_meloxicam_pred[3,3:4002]/emm_meloxicam_pred[1,3:4002])), comparison = 'PO:ML/IV:ML'))
emm_meloxicam_ratio <- rbind(emm_meloxicam_ratio, data.frame(AUC_ratio = as.numeric(t(emm_meloxicam_pred[4,3:4002]/emm_meloxicam_pred[2,3:4002])), comparison = 'PO:PP/IV:PP'))
emm_meloxicam_ratio <- rbind(emm_meloxicam_ratio, data.frame(AUC_ratio = as.numeric(t(emm_meloxicam_pred[4,3:4002]/emm_meloxicam_pred[3,3:4002])), comparison = 'PO:PP/PO:ML'))
emm_meloxicam_pred  <- melt(emm_meloxicam_pred)
emm_meloxicam_table <- data.frame(parameter = character(length = 8), L90CI = numeric(length = 8), L50CI = numeric(length = 8), point = numeric(length = 8), U50CI = numeric(length = 8), U90CI = numeric(length = 8))
emm_meloxicam_table[(emm_meloxicam_table$parameter == 'IV:ML'),2:6] <- quantile((emm_meloxicam_pred$value[(emm_meloxicam_pred$STAGE == 'Mid-lactation') & (emm_meloxicam_pred$ROUTE == 'IV')]), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
emm_meloxicam_table[(emm_meloxicam_table$parameter == 'PO:ML'),2:6] <- quantile((emm_meloxicam_pred$value[(emm_meloxicam_pred$STAGE == 'Mid-lactation') & (emm_meloxicam_pred$ROUTE == 'PO')]), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
emm_meloxicam_table[(emm_meloxicam_table$parameter == 'IV:PP'),2:6] <- quantile((emm_meloxicam_pred$value[(emm_meloxicam_pred$STAGE == 'Post-partum') & (emm_meloxicam_pred$ROUTE == 'IV')]), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
emm_meloxicam_table[(emm_meloxicam_table$parameter == 'PO:PP'),2:6] <- quantile((emm_meloxicam_pred$value[(emm_meloxicam_pred$STAGE == 'Post-partum') & (emm_meloxicam_pred$ROUTE == 'PO')]), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
emm_meloxicam_table[(emm_meloxicam_table$parameter == 'IV:PP/IV:ML'),2:6] <- quantile((emm_meloxicam_ratio$AUC_ratio[emm_meloxicam_ratio$comparison == 'IV:PP/IV:ML']), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
emm_meloxicam_table[(emm_meloxicam_table$parameter == 'PO:PP/PO:ML'),2:6] <- quantile((emm_meloxicam_ratio$AUC_ratio[emm_meloxicam_ratio$comparison == 'PO:PP/PO:ML']), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
emm_meloxicam_table[(emm_meloxicam_table$parameter == 'PO:ML/IV:ML'),2:6] <- quantile((emm_meloxicam_ratio$AUC_ratio[emm_meloxicam_ratio$comparison == 'PO:ML/IV:ML']), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
emm_meloxicam_table[(emm_meloxicam_table$parameter == 'PO:PP/IV:PP'),2:6] <- quantile((emm_meloxicam_ratio$AUC_ratio[emm_meloxicam_ratio$comparison == 'PO:PP/IV:PP']), c(0.05,0.25,0.50,0.75,0.95), na.rm = TRUE)
emm_AUCcov_plot     <- ggplot(data = emm_meloxicam_pred, aes(x = value, y = ROUTE:STAGE)) + stat_dotsinterval(quantiles = 100) + theme_bw() + xlab(expression(Dose-corrected~AUC[last])) + ylab('') + scale_y_discrete(limits = rev) + coord_cartesian(xlim = c(0,200))
emm_AUCratio_plot   <- ggplot(data = emm_meloxicam_ratio, aes(x = AUC_ratio, y = comparison)) + stat_dotsinterval(.width = c(0.5,0.9), quantiles = 100, interval_alpha = 0.8, point_alpha = 0.8) + theme_bw() + scale_x_log10(breaks = c(0.1,0.2,0.5,1,2,5,10)) + coord_cartesian(xlim = c(0.1,10)) + theme(panel.grid.minor = element_blank()) + geom_vline(xintercept = 1) + xlab(expression(Dose-corrected~AUC[last]~ratio)) + ylab('') + scale_y_discrete(limits = rev)
emm_pop_plots <- ggarrange(emm_AUCcov_plot,emm_AUCratio_plot, ncol = 1)
ggsave(emm_pop_plots, file = 'emm_pop_plots.svg', units = 'mm', width = 200, height = 200)

# Aggregate all of the population-level estimates and save them.
all_meloxicam_table <- cbind(rbind(pop_meloxicam_table, emm_meloxicam_table, meloxicam_trap_conf), data.frame(method = c(rep('GAM_pop',length(pop_meloxicam_table$parameter)), rep('GAM_mean',length(emm_meloxicam_table$parameter)), rep('trap_mean',length(meloxicam_trap_conf$parameter)))))
all_meloxicam_table[,2:6] <- round((all_meloxicam_table[,2:6]),2)
write.csv(all_meloxicam_table, file = 'all_meloxicam_table.csv', row.names = FALSE)
