# Andrew P. Woodward, 2022-23.
# Noncompartmental analysis of the reference theophylline dataset (https://doi.org/10.1007/bf01059644) as reported by Upton and utilized by Beal et al. in NONMEM (1992-2009) and Pinheiro and Bates (http://dx.doi.org/10.1080/10618600.1995.10474663).
# Here we utilize a hierarchical generalized additive model for the change in theophylline concentrations with time.
#     The hierarchical GAM are expressed following Pedersen et al. 2019 (https://doi.org/10.7717/peerj.6876).
#     Spline-based NCA was previously reported in (https://doi.org/10.1002/pst.336 and https://doi.org/10.1177/1471082X17706018); in this work we use the popular 'mgcv' implementation of regression splines via 'brms'.

# Load the required packages.
library(brms)
library(ggplot2)
library(ggdist)
library(ncappc)
library(future.apply)
library(rlist)
library(reshape2)
library(viridis)
library(ggpubr)
library(boot)
library(bootstrap)

# The time will be transformed as log(time+1); the time-concentration function is presumed to be multiplicative.
#     Note that without the +1 adjustment the time 0 observations will be dropped (https://doi.org/10.1002/pst.336).
#     For terminal phase observations, the log-time scale prevents the uncertainty from growing in-between the observations ('ballooning' effect), as the observations are similarly spaced on log-time scale.
#     Note that this makes the time scale somewhat nonlinear for small values of time, so the units of time must be chosen with care.
theophylline_data$Time_log1 <- log(theophylline_data$Time+1)
theophylline_data$Conc_log  <- log(theophylline_data$conc)
Time_nominal <- c(0,0.25,0.5,1,2,4,5,7,9,12,24)
for (i in 1:length(theophylline_data$Time)){
  theophylline_data$Time_nominal[i] <- Closest(Time_nominal, theophylline_data$Time[i])
}
theophylline_pop_summary <- data.frame(Time = numeric(length = length(unique(theophylline_data$Time_nominal))), conc_mean = numeric(length = length(unique(theophylline_data$Time_nominal))), conc_L95 = numeric(length = length(unique(theophylline_data$Time_nominal))), conc_U95 = numeric(length = length(unique(theophylline_data$Time_nominal))))
for (i in 1:length(unique(theophylline_data$Time_nominal))){
  iter_time <- (unique(theophylline_data$Time_nominal))[i]
  theophylline_pop_summary$Time[i] <- iter_time
  theophylline_pop_summary$conc_mean[i] <- mean(theophylline_data$conc[theophylline_data$Time_nominal == iter_time])
  iter_jack <- as.numeric((jackknife(theophylline_data$conc[theophylline_data$Time_nominal == iter_time], mean))$jack.se)
  theophylline_pop_summary$conc_L95[i] <- theophylline_pop_summary$conc_mean[i] - (1.96*iter_jack)
  theophylline_pop_summary$conc_U95[i] <- theophylline_pop_summary$conc_mean[i] + (1.96*iter_jack)
}

# Generalized additive model for the theophylline data, using a common smoother and subject-level smoothers with shared wiggliness (https://peerj.com/articles/6876/).
# This model uses the M1-penalty on the subject terms (penalizing the squared first-derivative).
# Then generate predicted concentrations and their uncertainty for visualization.
theoph_all_mod_m1 <- brm(Conc_log ~ 1 + s(Time_log1, k = 6) + s(Time_log1, Subject, bs = 'fs', m = 1, k = 6), data = theophylline_data, cores = 4, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 12), prior = set_prior('normal(0,0.5)', class = 'sigma') + set_prior('normal(0,5)', class = 'Intercept'))
summary(theoph_all_mod_m1)
bayes_R2(theoph_all_mod_m1)
conditional_effects(theoph_all_mod_m1)
new_theoph_all <- expand.grid(Time_log1 = seq(0,round(log(25+1),1),0.01), Subject = unique(theophylline_data$Subject))
prd_theoph_all_m1 <- cbind(new_theoph_all, fitted(theoph_all_mod_m1,newdata = new_theoph_all, probs = c(0.05,0.95)))
prd_theoph_all_m1$Time <- exp(prd_theoph_all_m1$Time_log1)-1
prd_theoph_all_m1$Estimate <- exp(prd_theoph_all_m1$Estimate)
prd_theoph_all_m1$Q5 <- exp(prd_theoph_all_m1$Q5)
prd_theoph_all_m1$Q95 <- exp(prd_theoph_all_m1$Q95)
m1_mod_plot_lin <- ggplot(theophylline_data, aes(x = Time, y = conc)) + geom_point() + facet_wrap(~Subject) + geom_line(data = prd_theoph_all_m1, aes(x = Time, y = Estimate), inherit.aes = FALSE) + geom_ribbon(data = prd_theoph_all_m1, aes(x = Time, y = Estimate, ymin = Q5, ymax = Q95), inherit.aes = FALSE, alpha = 0.3) + theme_bw() + ggtitle('M1 Penalty, constant-error') + ylab('Theophylline concentration (mg/L)') + xlab('Time (h)')
m1_mod_plot_log <- ggplot(theophylline_data, aes(x = Time, y = conc)) + geom_point() + facet_wrap(~Subject) + geom_line(data = prd_theoph_all_m1, aes(x = Time, y = Estimate), inherit.aes = FALSE) + geom_ribbon(data = prd_theoph_all_m1, aes(x = Time, y = Estimate, ymin = Q5, ymax = Q95), inherit.aes = FALSE, alpha = 0.3) + theme_bw() + ggtitle('M1 Penalty, constant-error') + ylab('Theophylline concentration (mg/L)') + xlab('Time (h)') + scale_x_log10()
ggarrange(m1_mod_plot_lin, m1_mod_plot_log, nrow = 2)
pp_check(theoph_all_mod_m1, ndraws = 30)

# The same model, but using the M2-penalty on the subject terms (penalizing the squared second-derivative).
#     This appears to have good performance for this case. Note that the default from 'mgcv' is (2m > d+1), so is 2 for the single-predictor smooth.
theoph_all_mod_m2 <- brm(Conc_log ~ 1 + s(Time_log1, k = 6) + s(Time_log1, Subject, bs = 'fs', m = 2, k = 6), data = theophylline_data, cores = 4, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 12), prior = set_prior('normal(0,0.5)', class = 'sigma') + set_prior('normal(0,5)', class = 'Intercept'))
summary(theoph_all_mod_m2)
bayes_R2(theoph_all_mod_m2, probs = c(0.05,0.95))
conditional_effects(theoph_all_mod_m2)
new_theoph_all <- expand.grid(Time_log1 = seq(0,round(log(25+1),1),0.01), Subject = unique(theophylline_data$Subject))
prd_theoph_all_m2 <- cbind(new_theoph_all, fitted(theoph_all_mod_m2,newdata = new_theoph_all, probs = c(0.05,0.25,0.75,0.95)))
prd_theoph_all_m2$Time <- exp(prd_theoph_all_m2$Time_log1)-1
prd_theoph_all_m2$Estimate <- exp(prd_theoph_all_m2$Estimate)
prd_theoph_all_m2$Q5  <- exp(prd_theoph_all_m2$Q5)
prd_theoph_all_m2$Q95 <- exp(prd_theoph_all_m2$Q95)
prd_theoph_all_m2$Q25 <- exp(prd_theoph_all_m2$Q25)
prd_theoph_all_m2$Q75 <- exp(prd_theoph_all_m2$Q75)
m2_mod_plot_lin <- ggplot(theophylline_data, aes(x = Time, y = conc)) + geom_point(alpha = 0.5) + facet_wrap(~Subject) + geom_line(data = prd_theoph_all_m2, aes(x = Time, y = Estimate), inherit.aes = FALSE) + geom_ribbon(data = prd_theoph_all_m2, aes(x = Time, y = Estimate, ymin = Q5, ymax = Q95), inherit.aes = FALSE, alpha = 0.2) + geom_ribbon(data = prd_theoph_all_m2, aes(x = Time, y = Estimate, ymin = Q25, ymax = Q75), inherit.aes = FALSE, alpha = 0.4) + theme_bw() + ggtitle('M2 Penalty, constant-error') + ylab('Theophylline concentration (mg/L)') + xlab('Time (h)') + coord_cartesian(ylim = c(0,15), xlim = c(0,24.5))
m2_mod_plot_log <- ggplot(theophylline_data, aes(x = Time, y = conc)) + geom_point(alpha = 0.5) + facet_wrap(~Subject) + geom_line(data = prd_theoph_all_m2, aes(x = Time, y = Estimate), inherit.aes = FALSE) + geom_ribbon(data = prd_theoph_all_m2, aes(x = Time, y = Estimate, ymin = Q5, ymax = Q95), inherit.aes = FALSE, alpha = 0.2) + geom_ribbon(data = prd_theoph_all_m2, aes(x = Time, y = Estimate, ymin = Q25, ymax = Q75), inherit.aes = FALSE, alpha = 0.4) + theme_bw() + ggtitle('M2 Penalty, constant-error') + ylab('Theophylline concentration (mg/L)') + xlab('Time (h)') + scale_x_log10()
sub_plots_theoph_m2 <- ggarrange(m2_mod_plot_lin, m2_mod_plot_log, nrow = 2)
ggsave(sub_plots_theoph_m2, file = 'sub_plots_theoph_m2.svg', height = 220, width = 200, units = 'mm')
pp_check(theoph_all_mod_m2, ndraws = 30)
theoph_residuals <- data.frame(time = theophylline_data$Time, time_log1 = theophylline_data$Time_log1, residual = residuals(theoph_all_mod_m2)[,1], model = 'm2')

# The residuals from the previous model are a bit further spread at the early time points, so try a time-varying error model.
#     This appears to have good performance for this case. Note that the default from 'mgcv' is (2m > d+1), so is 2 for the single-predictor smooth.
theoph_all_mod_sd <- brm(bf(Conc_log ~ 1 + s(Time_log1, k = 6) + s(Time_log1, Subject, bs = 'fs', m = 2, k = 6), sigma ~ 1 + Time_log1), data = theophylline_data, cores = 4, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 15), prior = set_prior('normal(-1,1)', class = 'Intercept', dpar = 'sigma') + set_prior('normal(0,2)', class = 'b', dpar = 'sigma') + set_prior('normal(0,5)', class = 'Intercept'))
summary(theoph_all_mod_sd)
bayes_R2(theoph_all_mod_sd, probs = c(0.05,0.95))
conditional_effects(theoph_all_mod_sd)
conditional_effects(theoph_all_mod_sd, dpar = 'sigma')
new_theoph_all <- expand.grid(Time_log1 = seq(0,round(log(25+1),1),0.01), Subject = unique(theophylline_data$Subject))
prd_theoph_all_sd <- cbind(new_theoph_all, fitted(theoph_all_mod_sd,newdata = new_theoph_all, probs = c(0.05,0.25,0.75,0.95)))
prd_theoph_all_sd$Time <- exp(prd_theoph_all_sd$Time_log1)-1
prd_theoph_all_sd$Estimate <- exp(prd_theoph_all_sd$Estimate)
prd_theoph_all_sd$Q5  <- exp(prd_theoph_all_sd$Q5)
prd_theoph_all_sd$Q95 <- exp(prd_theoph_all_sd$Q95)
prd_theoph_all_sd$Q25  <- exp(prd_theoph_all_sd$Q25)
prd_theoph_all_sd$Q75 <- exp(prd_theoph_all_sd$Q75)
sd_mod_plot_lin <- ggplot(theophylline_data, aes(x = Time, y = conc)) + geom_point(alpha = 0.5) + facet_wrap(~Subject) + geom_line(data = prd_theoph_all_sd, aes(x = Time, y = Estimate), inherit.aes = FALSE) + geom_ribbon(data = prd_theoph_all_sd, aes(x = Time, y = Estimate, ymin = Q5, ymax = Q95), inherit.aes = FALSE, alpha = 0.2) + geom_ribbon(data = prd_theoph_all_sd, aes(x = Time, y = Estimate, ymin = Q25, ymax = Q75), inherit.aes = FALSE, alpha = 0.4) + theme_bw() + ggtitle('M2 Penalty, varying-error') + ylab('Theophylline concentration (mg/L)') + xlab('Time (h)') + coord_cartesian(ylim = c(0,15), xlim = c(0,24.5))
sd_mod_plot_log <- ggplot(theophylline_data, aes(x = Time, y = conc)) + geom_point(alpha = 0.5) + facet_wrap(~Subject) + geom_line(data = prd_theoph_all_sd, aes(x = Time, y = Estimate), inherit.aes = FALSE) + geom_ribbon(data = prd_theoph_all_sd, aes(x = Time, y = Estimate, ymin = Q5, ymax = Q95), inherit.aes = FALSE, alpha = 0.2) + geom_ribbon(data = prd_theoph_all_sd, aes(x = Time, y = Estimate, ymin = Q25, ymax = Q75), inherit.aes = FALSE, alpha = 0.4) + theme_bw() + ggtitle('M2 Penalty, varying-error') + ylab('Theophylline concentration (mg/L)') + xlab('Time (h)') + scale_x_log10()
sub_plots_theoph_sd <- ggarrange(sd_mod_plot_lin, sd_mod_plot_log, nrow = 2)
ggsave(sub_plots_theoph_sd, file = 'sub_plots_theoph_sd.svg', height = 220, width = 200, units = 'mm')
pp_check(theoph_all_mod_sd, ndraws = 30)
theoph_residuals <- rbind(theoph_residuals, data.frame(time = theophylline_data$Time, time_log1 = theophylline_data$Time_log1, residual = residuals(theoph_all_mod_sd)[,1], model = 'sd'))

# Now add between-subject variation in wiggliness (https://doi.org/10.7717/peerj.6876), which allows for greater flexibility in the form of the concentration-time function across subjects.
#     A subject-level hierarchical intercept is included to accommodate location shift (without it, at time 0 the predictions are identical).
#     This visibly improves the fit for some individuals compared with the shared-wiggliness model, and may be more realistic than the time-varying residual variation.
theoph_all_mod_by <- brm(Conc_log ~ 1 + s(Time_log1, k = 6) + s(Time_log1, by = Subject, m = 1, k = 6) + (1|Subject), data = theophylline_data, cores = 4, chains = 4, control = list(adapt_delta = 0.99, max_treedepth = 12), prior = set_prior('normal(0,0.5)', class = 'sigma') + set_prior('normal(0,0.5)', class = 'sd') + set_prior('normal(0,5)', class = 'Intercept'))
summary(theoph_all_mod_by)
bayes_R2(theoph_all_mod_by, probs = c(0.05,0.95))
conditional_effects(theoph_all_mod_by)
new_theoph_all <- expand.grid(Time_log1 = seq(0,round(log(25+1),1),0.01), Subject = unique(theophylline_data$Subject))
prd_theoph_all_by <- cbind(new_theoph_all, fitted(theoph_all_mod_by,newdata = new_theoph_all, probs = c(0.05,0.25,0.75,0.95)))
prd_theoph_all_by$Time <- exp(prd_theoph_all_by$Time_log1)-1
prd_theoph_all_by$Estimate <- exp(prd_theoph_all_by$Estimate)
prd_theoph_all_by$Q5 <- exp(prd_theoph_all_by$Q5)
prd_theoph_all_by$Q95 <- exp(prd_theoph_all_by$Q95)
prd_theoph_all_by$Q25 <- exp(prd_theoph_all_by$Q25)
prd_theoph_all_by$Q75 <- exp(prd_theoph_all_by$Q75)
by_mod_plot_lin <- ggplot(theophylline_data, aes(x = Time, y = conc)) + geom_point(alpha = 0.5) + facet_wrap(~Subject) + geom_line(data = prd_theoph_all_by, aes(x = Time, y = Estimate), inherit.aes = FALSE) + geom_ribbon(data = prd_theoph_all_by, aes(x = Time, y = Estimate, ymin = Q5, ymax = Q95), inherit.aes = FALSE, alpha = 0.2) + geom_ribbon(data = prd_theoph_all_by, aes(x = Time, y = Estimate, ymin = Q25, ymax = Q75), inherit.aes = FALSE, alpha = 0.4) + theme_bw() + ggtitle('Variable Wiggliness') + ylab('Theophylline concentration (mg/L)') + xlab('Time (h)') + coord_cartesian(ylim = c(0,15), xlim = c(0,24.5))
by_mod_plot_log <- ggplot(theophylline_data, aes(x = Time, y = conc)) + geom_point(alpha = 0.5) + facet_wrap(~Subject) + geom_line(data = prd_theoph_all_by, aes(x = Time, y = Estimate), inherit.aes = FALSE) + geom_ribbon(data = prd_theoph_all_by, aes(x = Time, y = Estimate, ymin = Q5, ymax = Q95), inherit.aes = FALSE, alpha = 0.2) + geom_ribbon(data = prd_theoph_all_by, aes(x = Time, y = Estimate, ymin = Q25, ymax = Q75), inherit.aes = FALSE, alpha = 0.4) + theme_bw() + ggtitle('Variable Wiggliness') + ylab('Theophylline concentration (mg/L)') + xlab('Time (h)') + scale_x_log10()
sub_plots_theoph_by <- ggarrange(by_mod_plot_lin, by_mod_plot_log, nrow = 2)
ggsave(sub_plots_theoph_by, file = 'sub_plots_theoph_by.svg', height = 220, width = 200, units = 'mm')
pp_check(theoph_all_mod_sd, ndraws = 30)
theoph_residuals <- rbind(theoph_residuals, data.frame(time = theophylline_data$Time, time_log1 = theophylline_data$Time_log1, residual = residuals(theoph_all_mod_by)[,1], model = 'by'))

# Generate residual analyses across models to demonstrate the relative performance of each.
theoph_sd_pred <- data.frame(time = seq(0,25,0.1), time_log1 = log(seq(0,25,0.1)+1), model = 'm2', Estimate = ((summary(theoph_all_mod_m2,probs = c(0.05,0.95)))$spec_pars[1]), `Est.Error` = as.numeric((summary(theoph_all_mod_m2,probs = c(0.05,0.95)))$spec_pars[2]), Q5 = as.numeric((summary(theoph_all_mod_m2,probs = c(0.05,0.95)))$spec_pars[3]), Q95 = as.numeric((summary(theoph_all_mod_m2,probs = c(0.05,0.95)))$spec_pars[4]))
theoph_sd_pred <- rbind(theoph_sd_pred, data.frame(time = seq(0,25,0.1), time_log1 = log(seq(0,25,0.1)+1), model = 'by', Estimate = ((summary(theoph_all_mod_by,probs = c(0.05,0.95)))$spec_pars[1]), `Est.Error` = as.numeric((summary(theoph_all_mod_by,probs = c(0.05,0.95)))$spec_pars[2]), Q5 = as.numeric((summary(theoph_all_mod_by,probs = c(0.05,0.95)))$spec_pars[3]), Q95 = as.numeric((summary(theoph_all_mod_by,probs = c(0.05,0.95)))$spec_pars[4])))
theoph_sd_pred <- rbind(theoph_sd_pred, cbind((data.frame(time = seq(0,25,0.1), time_log1 = log(seq(0,25,0.1)+1), model = 'sd')), fitted(theoph_all_mod_sd, newdata = data.frame(Time_log1 = log(seq(0,25,0.1)+1), Subject = NA), dpar = 'sigma', re_formula = NA, probs = c(0.05,0.95))))
theoph_sd_pred$model <- factor(theoph_sd_pred$model, levels = c('m2','sd','by'))
theoph_residuals$model <- factor(theoph_residuals$model, levels = c('m2','sd','by'))
theoph_sd_lin  <- ggplot(data = theoph_sd_pred, aes(x = time, y = Estimate)) + geom_line() + geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.3) + facet_wrap(~model) + theme_bw() + ylab('Conditional standard deviation') + xlab('Time (h)')
theoph_sd_log  <- ggplot(data = theoph_sd_pred, aes(x = time_log1, y = Estimate)) + geom_line() + geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.3) + facet_wrap(~model) + theme_bw() + ylab('Conditional standard deviation') + xlab('log(Time+1)')
theoph_res_lin <- ggplot(data = theoph_residuals, aes(x = time, y = residual)) + geom_point(alpha = 0.3) + facet_wrap(~model) + theme_bw() + ylim(c(-1.5,1.5)) + xlab('Time (h)') + ylab('Residual') + geom_line(data = theoph_sd_pred, aes(x = time, y = 0-Estimate), inherit.aes = FALSE, alpha = 0.6) + geom_line(data = theoph_sd_pred, aes(x = time, y = 0+Estimate), inherit.aes = FALSE, alpha = 0.6) + geom_line(data = theoph_sd_pred, aes(x = time, y = 0-2*(Estimate)), inherit.aes = FALSE, alpha = 0.6) + geom_line(data = theoph_sd_pred, aes(x = time, y = 0+2*(Estimate)), inherit.aes = FALSE, alpha = 0.6)
theoph_res_log <- ggplot(data = theoph_residuals, aes(x = time_log1, y = residual)) + geom_point(alpha = 0.3) + facet_wrap(~model) + theme_bw() + ylim(c(-1.5,1.5)) + xlab('log(Time+1)') + ylab('Residual')  + geom_line(data = theoph_sd_pred, aes(x = time_log1, y = 0-Estimate), inherit.aes = FALSE, alpha = 0.6) + geom_line(data = theoph_sd_pred, aes(x = time_log1, y = 0+Estimate), inherit.aes = FALSE, alpha = 0.6) + geom_line(data = theoph_sd_pred, aes(x = time_log1, y = 0-2*(Estimate)), inherit.aes = FALSE, alpha = 0.6) + geom_line(data = theoph_sd_pred, aes(x = time_log1, y = 0+2*(Estimate)), inherit.aes = FALSE, alpha = 0.6)
theoph_residuals_plots <- ggarrange(theoph_res_lin, theoph_sd_lin, theoph_res_log, theoph_sd_log, nrow = 4, ncol = 1)
ggsave(theoph_residuals_plots, file = 'theoph_residuals_plots.svg', width = 200, height = 220, units = 'mm')

# For each model, generate a graphic for the population smooth and the data.
theoph_conc_pop_sd <- cbind(data.frame(Time_log1 = seq(0,round(log(25+1),1),0.01), exp(fitted(theoph_all_mod_sd, newdata = data.frame(Time_log1 = seq(0,round(log(25+1),1),0.01), Subject = NA), re_formula = NA, probs = c(0.05,0.25,0.75,0.95)))))
theoph_conc_pop_sd$Time <- exp(theoph_conc_pop_sd$Time_log1)-1
theoph_conc_pop_m2 <- cbind(data.frame(Time_log1 = seq(0,round(log(25+1),1),0.01), exp(fitted(theoph_all_mod_m2, newdata = data.frame(Time_log1 = seq(0,round(log(25+1),1),0.01), Subject = NA), re_formula = NA, probs = c(0.05,0.25,0.75,0.95)))))
theoph_conc_pop_m2$Time <- exp(theoph_conc_pop_m2$Time_log1)-1
theoph_conc_pop_by <- cbind(data.frame(Time_log1 = seq(0,round(log(25+1),1),0.01), exp(fitted(theoph_all_mod_by, newdata = data.frame(Time_log1 = seq(0,round(log(25+1),1),0.01), Subject = NA), re_formula = NA, probs = c(0.05,0.25,0.75,0.95)))))
theoph_conc_pop_by$Time <- exp(theoph_conc_pop_m2$Time_log1)-1
theoph_pop_smooth_sd <- ggplot(data = theoph_conc_pop_sd, aes(x = Time, y = Estimate)) + geom_line() + geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.2) + geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = 0.3) + theme_bw() + geom_point(data = theophylline_data, aes(x = Time, y = conc), inherit.aes = FALSE, alpha = 0.2) + geom_line(data = theophylline_data, aes(x = Time, y = conc, group = Subject), inherit.aes = FALSE, alpha = 0.1) + scale_y_continuous(breaks = c(0,5,10)) + scale_x_continuous(breaks = c(0,10,20)) + theme(panel.grid.minor = element_blank()) + xlab('Time (h)') + ylab('Concentration (mg/L)') + coord_cartesian(ylim = c(0,12), xlim = c(0,24.5)) + ggtitle('M2 Penalty, varying-error')
theoph_pop_smooth_m2 <- ggplot(data = theoph_conc_pop_m2, aes(x = Time, y = Estimate)) + geom_line() + geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.2) + geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = 0.3) + theme_bw() + geom_point(data = theophylline_data, aes(x = Time, y = conc), inherit.aes = FALSE, alpha = 0.2) + geom_line(data = theophylline_data, aes(x = Time, y = conc, group = Subject), inherit.aes = FALSE, alpha = 0.1) + scale_y_continuous(breaks = c(0,5,10)) + scale_x_continuous(breaks = c(0,10,20)) + theme(panel.grid.minor = element_blank()) + xlab('Time (h)') + ylab('Concentration (mg/L)') + coord_cartesian(ylim = c(0,12), xlim = c(0,24.5)) + ggtitle('M2 Penalty, constant-error')
theoph_pop_smooth_by <- ggplot(data = theoph_conc_pop_by, aes(x = Time, y = Estimate)) + geom_line() + geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.2) + geom_ribbon(aes(ymin = Q25, ymax = Q75), alpha = 0.3) + theme_bw() + geom_point(data = theophylline_data, aes(x = Time, y = conc), inherit.aes = FALSE, alpha = 0.2) + geom_line(data = theophylline_data, aes(x = Time, y = conc, group = Subject), inherit.aes = FALSE, alpha = 0.1) + scale_y_continuous(breaks = c(0,5,10)) + scale_x_continuous(breaks = c(0,10,20)) + theme(panel.grid.minor = element_blank()) + xlab('Time (h)') + ylab('Concentration (mg/L)') + coord_cartesian(ylim = c(0,12), xlim = c(0,24.5)) + ggtitle('Variable Wiggliness')
theoph_pop_plots <- ggarrange(theoph_pop_smooth_m2, theoph_pop_smooth_sd, theoph_pop_smooth_by, ncol = 1, nrow = 3)
ggsave(theoph_pop_plots, file = 'theoph_pop_plots.svg', units = 'mm', height = 220, width = 200)

# For comparison, we will add the AUC_last as determined with the classical trapezoidal method.
AUC_est_total <- data.frame(AUClast_GAM_sd = numeric(length(unique(theophylline_data$Subject))))
AUC_est_total$AUClast_GAM_M2 <- 0
AUC_est_total$AUClast_trap <- 0
AUC_est_total$AUCinf_trap <- 0
AUC_est_total$Cmax_GAM_sd <- 0
AUC_est_total$Tmax_GAM_sd <- 0
AUC_est_total$Cmax_trap <- 0
AUC_est_total$Tmax_trap <- 0
for (i in 1:length(unique(theophylline_data$Subject))){
  subject_ind <- unique(theophylline_data$Subject)[i]
  AUC_est_total$AUClast_trap[i] <- (ncappc(obsFile = data.frame(ID = theophylline_data$Subject[theophylline_data$Subject == subject_ind], TIME = theophylline_data$Time[theophylline_data$Subject == subject_ind], DV = theophylline_data$conc[theophylline_data$Subject == subject_ind]), onlyNCA = TRUE, printOut = FALSE, extrapolate = TRUE, evid = FALSE))$ncaOutput$AUClast
  AUC_est_total$AUCinf_trap[i] <- (ncappc(obsFile = data.frame(ID = theophylline_data$Subject[theophylline_data$Subject == subject_ind], TIME = theophylline_data$Time[theophylline_data$Subject == subject_ind], DV = theophylline_data$conc[theophylline_data$Subject == subject_ind]), onlyNCA = TRUE, printOut = FALSE, extrapolate = TRUE, evid = FALSE))$ncaOutput$AUCINF_obs
  AUC_est_total$Cmax_trap[i] <- (ncappc(obsFile = data.frame(ID = theophylline_data$Subject[theophylline_data$Subject == subject_ind], TIME = theophylline_data$Time[theophylline_data$Subject == subject_ind], DV = theophylline_data$conc[theophylline_data$Subject == subject_ind]), onlyNCA = TRUE, printOut = FALSE, extrapolate = TRUE, evid = FALSE))$ncaOutput$Cmax
  AUC_est_total$Tmax_trap[i] <- (ncappc(obsFile = data.frame(ID = theophylline_data$Subject[theophylline_data$Subject == subject_ind], TIME = theophylline_data$Time[theophylline_data$Subject == subject_ind], DV = theophylline_data$conc[theophylline_data$Subject == subject_ind]), onlyNCA = TRUE, printOut = FALSE, extrapolate = TRUE, evid = FALSE))$ncaOutput$Tmax
}

# Generate point estimates of the AUC by numeric integration to 24h, using the standard M2-penalty, or the varying-error model.
single_AUC_int <- function(time_int, subject, model){
  time_int_data <- data.frame(Time_log1 = log(time_int+1), Subject = subject)
  concentration_estimate <- (exp(fitted(model, newdata = time_int_data)[,1]))
  return(concentration_estimate)
}
for (i in 1:length(unique(theophylline_data$Subject))){
  AUC_est_total$AUClast_GAM_sd[i] <- as.numeric(integrate(single_AUC_int,0,24, rel.tol = 0.001, abs.tol = 0.001, subject = unique(theophylline_data$Subject)[i], model = theoph_all_mod_sd)[1])
  AUC_est_total$AUClast_GAM_M2[i] <- as.numeric(integrate(single_AUC_int,0,24, rel.tol = 0.001, abs.tol = 0.001, subject = unique(theophylline_data$Subject)[i], model = theoph_all_mod_m2)[1])
  int_optim <- optimize(single_AUC_int,c(0,24), maximum = TRUE, subject = unique(theophylline_data$Subject)[i], model = theoph_all_mod_sd)
  AUC_est_total$Cmax_GAM_sd[i] <- int_optim$objective
  AUC_est_total$Tmax_GAM_sd[i] <- int_optim$maximum
}
AUC_est_total$subject <- unique(theophylline_data$Subject)

# Uncertainty statements about the AUC can be obtained simply from the samples. The integration is conducted on the parameters of each sample rather than the posterior mean parameters.
#     This computation is really expensive, probably partially because of the setup cost for calling to the Stan model to make predictions during integration.
#     Here a parallel operation is implemented via package 'future' (https://arxiv.org/abs/2008.00553) to improve time efficiency.
#     The analysis was tested on an AMD Ryzen 7 8-core 16-thread CPU, so 8 workers are requested in this example.
# Integration function to return the predicted concentration at any time.
single_AUC_theoph <- function(time_int, samples, subject, model){
  time_int_data <- data.frame(Time_log1 = log(time_int+1), Subject = subject)
  if(is.na(subject)){
    concentration_estimate <- (exp(fitted(model, newdata = time_int_data, draw_ids = samples, summary = FALSE, re.formula = NA, allow_new_levels = TRUE)))
  }
  else{
    concentration_estimate <- (exp(fitted(model, newdata = time_int_data, draw_ids = samples, summary = FALSE)))  
  }
  return(concentration_estimate)
}

# Call the parallel (12 workers) implementation of the per-subject, per-sample integration, and aggregate it into a sample-per-row table, incorporating the subject predictors.
#     Firstly for the varying-SD model.
plan(multisession, workers = 12)
AUC_est_list <- future_lapply((unique(theophylline_data$Subject)), function(subject_ind) {
  AUC_est_theoph_sub <- numeric(ndraws(pred_subj_brms))
  for (h in 1:ndraws(pred_subj_brms)){
    AUC_est_theoph_sub[h] <- as.numeric(integrate(single_AUC_theoph, 0, 24, samples = h, subject = subject_ind, model = theoph_all_mod_sd)[1])
  }
  AUC_est_theoph_sub
})
AUC_est_list_theoph <- AUC_est_list
for(i in 1: length(AUC_est_list_theoph)){
  AUC_est_list_theoph[[i]] <- data.frame(AUC = AUC_est_list_theoph[[i]], subject = (unique(theoph_data$Subject))[i])
}
AUC_est_samples_sd <- list.rbind(AUC_est_list)

# Visualize the posterior distributions, and the trapezoidal rule estimates of the AUC.
AUC_est_theoph_samples1 <- data.frame(AUC = melt(AUC_est_samples_sd, id.vars = NA)[,3])
AUC_est_theoph_samples1$subject <- 0
AUC_est_theoph_samples1$subject <- rep((unique(theophylline_data$Subject)),1)
AUC_est_theoph_samples1$subject <- factor(AUC_est_theoph_samples1$subject)
AUClast_theoph_plot_sd <- ggplot(data = AUC_est_theoph_samples1, aes(x = AUC)) + facet_wrap(~subject) + stat_histinterval(alpha = 0.8) + geom_vline(data = AUC_est_total, aes(xintercept = AUClast_trap), linetype = 'dashed', alpha = 0.7) + xlab(expression(AUC['24h']~(mg.hr/L))) + ylab('') + theme_bw() + coord_cartesian(xlim = c(0,200)) + scale_y_continuous(breaks = NULL) + ggtitle('M2 Penalty, varying-error')
ggsave(AUClast_theoph_plot_sd, file = 'AUClast_theoph_plot_sd.svg', units = 'mm', width = 200, height = 200)
rm(AUC_est_list,AUC_est_samples_sd)

# Call the parallel (12 workers) implementation of the per-subject, per-sample integration, and aggregate it into a sample-per-row table, incorporating the subject predictors.
#     Now for the M2 penalty model.
plan(multisession, workers = 12)
AUC_est_list <- future_lapply((unique(theophylline_data$Subject)), function(subject_ind) {
  AUC_est_theoph_sub <- numeric(ndraws(pred_subj_brms))
  for (h in 1:ndraws(pred_subj_brms)){
    AUC_est_theoph_sub[h] <- as.numeric(integrate(single_AUC_theoph, 0, 24, samples = h, subject = subject_ind, model = theoph_all_mod_m2)[1])
  }
  AUC_est_theoph_sub
})
AUC_est_list_theoph <- AUC_est_list
for(i in 1: length(AUC_est_list_theoph)){
  AUC_est_list_theoph[[i]] <- data.frame(AUC = AUC_est_list_theoph[[i]], subject = (unique(theoph_data$Subject))[i])
}
AUC_est_samples_m2 <- list.rbind(AUC_est_list)

# Visualize the posterior distributions, and the trapezoidal rule estimates of the AUC.
AUC_est_theoph_samples2 <- data.frame(AUC = melt(AUC_est_samples_m2, id.vars = NA)[,3])
AUC_est_theoph_samples2$subject <- 0
AUC_est_theoph_samples2$subject <- rep((unique(theophylline_data$Subject)),1)
AUC_est_theoph_samples2$subject <- factor(AUC_est_theoph_samples2$subject)
AUClast_theoph_plot_m2 <- ggplot(data = AUC_est_theoph_samples2, aes(x = AUC)) + facet_wrap(~subject) + stat_histinterval(alpha = 0.8) + geom_vline(data = AUC_est_total, aes(xintercept = AUClast_trap), linetype = 'dashed', alpha = 0.7) + xlab(expression(AUC['24h']~(mg.hr/L))) + ylab('') + theme_bw() + coord_cartesian(xlim = c(0,200)) + scale_y_continuous(breaks = NULL) + ggtitle('M2 Penalty, constant-error')
ggsave(AUClast_theoph_plot_m2, file = 'AUClast_theoph_plot_m2.svg', units = 'mm', width = 200, height = 200)
rm(AUC_est_list,AUC_est_samples_m2)

# Call the parallel (12 workers) implementation of the per-subject, per-sample integration, and aggregate it into a sample-per-row table, incorporating the subject predictors.
#     Now for the M2 penalty model.
plan(multisession, workers = 12)
AUC_est_list <- future_lapply((unique(theophylline_data$Subject)), function(subject_ind) {
  AUC_est_theoph_sub <- numeric(ndraws(pred_subj_brms))
  for (h in 1:ndraws(pred_subj_brms)){
    AUC_est_theoph_sub[h] <- as.numeric(integrate(single_AUC_theoph, 0, 24, samples = h, subject = subject_ind, model = theoph_all_mod_by)[1])
  }
  AUC_est_theoph_sub
})
AUC_est_list_theoph <- AUC_est_list
for(i in 1: length(AUC_est_list_theoph)){
  AUC_est_list_theoph[[i]] <- data.frame(AUC = AUC_est_list_theoph[[i]], subject = (unique(theoph_data$Subject))[i])
}
AUC_est_samples_by <- list.rbind(AUC_est_list)

# Visualize the posterior distributions, and the trapezoidal rule estimates of the AUC.
AUC_est_theoph_samples3 <- data.frame(AUC = melt(AUC_est_samples_by, id.vars = NA)[,3])
AUC_est_theoph_samples3$subject <- 0
AUC_est_theoph_samples3$subject <- rep((unique(theophylline_data$Subject)),1)
AUC_est_theoph_samples3$subject <- factor(AUC_est_theoph_samples2$subject)
AUClast_theoph_plot_by <- ggplot(data = AUC_est_theoph_samples3, aes(x = AUC)) + facet_wrap(~subject) + stat_histinterval(alpha = 0.8) + geom_vline(data = AUC_est_total, aes(xintercept = AUClast_trap), linetype = 'dashed', alpha = 0.7) + xlab(expression(AUC['24h']~(mg.hr/L))) + ylab('') + theme_bw() + coord_cartesian(xlim = c(0,200)) + scale_y_continuous(breaks = NULL) + ggtitle('Variable Wiggliness')
ggsave(AUClast_theoph_plot_by, file = 'AUClast_theoph_plot_by.svg', units = 'mm', width = 200, height = 200)
rm(AUC_est_list,AUC_est_samples_by)


# Aggregate the subject-level estimates and intervals of the AUC.
AUC_quantile_est <- data.frame(subject = numeric(length(unique(theophylline_data$Subject))), AUC_sd_q05 = numeric(length(unique(theophylline_data$Subject))), AUC_sd_q50 = numeric(length(unique(theophylline_data$Subject))), AUC_sd_q95 = numeric(length(unique(theophylline_data$Subject))), AUC_m2_q05 = numeric(length(unique(theophylline_data$Subject))), AUC_m2_q50 = numeric(length(unique(theophylline_data$Subject))), AUC_m2_q95 = numeric(length(unique(theophylline_data$Subject))))
for(i in 1:length(unique(theophylline_data$Subject))){
  AUC_quantile_est$subject[i] <- (unique(theophylline_data$Subject))[i]
  AUC_quantile_est[i,2:4] <- quantile(AUC_est_theoph_samples1$AUC[(AUC_est_theoph_samples1$subject == (unique(theophylline_data$Subject))[i])],c(0.05,0.50,0.95))
  AUC_quantile_est[i,5:7] <- quantile(AUC_est_theoph_samples2$AUC[(AUC_est_theoph_samples2$subject == (unique(theophylline_data$Subject))[i])],c(0.05,0.50,0.95))
}

# For comparison to the HGAM approach, here generate percentile bootstrap estimates of the 'population' (average) AUC.
boot_AUC_fun <- function(input_data, input_index){
  boot_data <- input_data[input_index]
  boot_out <- mean(log(boot_data))
}
AUC_trap_boot <- boot(AUC_est_total$AUClast_trap, boot_AUC_fun, 1000)
AUC_trap_CI <- boot.ci(AUC_trap_boot, type = 'perc', conf = 0.9)
c(exp(AUC_trap_CI$t0),exp(AUC_trap_CI$percent[4]),exp(AUC_trap_CI$percent[5]))

theoph_parameter_table <- data.frame(parameter = character(length = 9), model = character(length = 9), Q5 = numeric(length = 9), Q25 = numeric(length = 9), Q50 = numeric(length = 9), Q75 = numeric(length = 9), Q95 = numeric(length = 9))
# Credible intervals and the posterior median average AUC using the HGAM with the M2-penalty model.
AUC_est_theoph_pop_m2  <- numeric(ndraws(theoph_all_mod_m2))
Tmax_est_theoph_pop_m2 <- numeric(ndraws(theoph_all_mod_m2))
Cmax_est_theoph_pop_m2 <- numeric(ndraws(theoph_all_mod_m2))
for (h in 1:ndraws(theoph_all_mod_m2)){
  AUC_est_theoph_pop_m2[h] <- as.numeric(integrate(single_AUC_theoph, 0, 24, samples = h, subject = NA, model = theoph_all_mod_m2)[1])
  iter_max <- optimize(function(x) exp(fitted(theoph_all_mod_m2, newdata = data.frame(Time_log1 = log(x+1), Subject = NA), draw_ids = h)[1]), lower = 0, upper = 12, maximum = TRUE)
  Tmax_est_theoph_pop_m2[h] <- as.numeric(iter_max$maximum)
  Cmax_est_theoph_pop_m2[h] <- as.numeric(iter_max$objective)
}
theoph_parameter_table[1,3:7] <- quantile(AUC_est_theoph_pop_m2,c(0.05,0.25,0.50,0.75,0.95))
theoph_parameter_table[4,3:7] <- quantile(Tmax_est_theoph_pop_m2,c(0.05,0.25,0.50,0.75,0.95))
theoph_parameter_table[7,3:7] <- quantile(Cmax_est_theoph_pop_m2,c(0.05,0.25,0.50,0.75,0.95))

# Credible intervals and the posterior median average AUC using the HGAM with the varying error model.
AUC_est_theoph_pop_sd  <- numeric(ndraws(theoph_all_mod_sd))
Tmax_est_theoph_pop_sd <- numeric(ndraws(theoph_all_mod_sd))
Cmax_est_theoph_pop_sd <- numeric(ndraws(theoph_all_mod_sd))
for (h in 1:ndraws(theoph_all_mod_sd)){
  AUC_est_theoph_pop_sd[h] <- as.numeric(integrate(single_AUC_theoph, 0, 24, samples = h, subject = NA, model = theoph_all_mod_sd)[1])
  iter_max <- optimize(function(x) exp(fitted(theoph_all_mod_sd, newdata = data.frame(Time_log1 = log(x+1), Subject = NA), draw_ids = h)[1]), lower = 0, upper = 12, maximum = TRUE)
  Tmax_est_theoph_pop_sd[h] <- as.numeric(iter_max$maximum)
  Cmax_est_theoph_pop_sd[h] <- as.numeric(iter_max$objective)
}
theoph_parameter_table[2,3:7] <- quantile(AUC_est_theoph_pop_sd,c(0.05,0.25,0.50,0.75,0.95))
theoph_parameter_table[5,3:7] <- quantile(Tmax_est_theoph_pop_sd,c(0.05,0.25,0.50,0.75,0.95))
theoph_parameter_table[8,3:7] <- quantile(Cmax_est_theoph_pop_sd,c(0.05,0.25,0.50,0.75,0.95))

# Credible intervals and the posterior median average AUC using the HGAM with the varying-wiggliness model.
AUC_est_theoph_pop_by  <- numeric(ndraws(theoph_all_mod_sd))
Tmax_est_theoph_pop_by <- numeric(ndraws(theoph_all_mod_sd))
Cmax_est_theoph_pop_by <- numeric(ndraws(theoph_all_mod_sd))
for (h in 1:ndraws(theoph_all_mod_by)){
  AUC_est_theoph_pop_by[h] <- as.numeric(integrate(single_AUC_theoph, 0, 24, samples = h, subject = NA, model = theoph_all_mod_by)[1])
  iter_max <- optimize(function(x) exp(fitted(theoph_all_mod_by, newdata = data.frame(Time_log1 = log(x+1), Subject = NA), re_formula = NA, draw_ids = h)[1]), lower = 0, upper = 12, maximum = TRUE)
  Tmax_est_theoph_pop_by[h] <- as.numeric(iter_max$maximum)
  Cmax_est_theoph_pop_by[h] <- as.numeric(iter_max$objective)
}
theoph_parameter_table[3,3:7] <- quantile(AUC_est_theoph_pop_by,c(0.05,0.25,0.50,0.75,0.95))
theoph_parameter_table[6,3:7] <- quantile(Tmax_est_theoph_pop_by,c(0.05,0.25,0.50,0.75,0.95))
theoph_parameter_table[9,3:7] <- quantile(Cmax_est_theoph_pop_by,c(0.05,0.25,0.50,0.75,0.95))
theoph_parameter_table$parameter <- rep(c('AUC','Tmax','Cmax'), each = 3)
theoph_parameter_table$model <- rep(c('m2','sd','by'),3)
theoph_parameter_table$Q5  <- round(theoph_parameter_table$Q5,  2)
theoph_parameter_table$Q25 <- round(theoph_parameter_table$Q25,  2)
theoph_parameter_table$Q50 <- round(theoph_parameter_table$Q50, 2)
theoph_parameter_table$Q75 <- round(theoph_parameter_table$Q75,  2)
theoph_parameter_table$Q95 <- round(theoph_parameter_table$Q95, 2)
write.csv(theoph_parameter_table, file = 'theoph_parameter_table.csv', row.names = FALSE)

# Aggregate those population AUC estimates and visualize, including the bootstrap confidence interval.
AUC_est_pop_all <- rbind(data.frame(AUC = AUC_est_theoph_pop_m2, model = 'm2'), data.frame(AUC = AUC_est_theoph_pop_sd, model = 'sd'), data.frame(AUC = AUC_est_theoph_pop_by, model = 'by'), data.frame(AUC = NA, model = 'conf'))
AUC_est_pop_all$model <- factor(AUC_est_pop_all$model, levels = c('m2','sd','by','conf'))
AUC_est_pop_con <- data.frame(model = 'conf', AUC5 = exp(AUC_trap_CI$percent[4]), AUC50 = exp(AUC_trap_CI$t0), AUC95 = exp(AUC_trap_CI$percent[5]))
theoph_pop_est_plot <- ggplot(data = AUC_est_pop_all, aes(x = AUC, y = model)) + stat_dotsinterval(quantiles = 100, .width = c(0.5,0.9)) + theme_bw() + geom_errorbar(data = AUC_est_pop_con, aes(x = AUC50, xmin = AUC5, xmax = AUC95, y = model), width = 0.1) + geom_point(data = AUC_est_pop_con, aes(x = AUC50, y = model), size = 3) + scale_x_continuous(limits = c(0,150)) + scale_y_discrete(labels = c('Constant error', 'Time-varying error', 'Varying wiggliness', 'Two-stage bootstrap')) + ylab('') + xlab(expression(paste('AUC (mg', '\U00B7', 'hr', '\U00B7', L^-1, ')', sep = '')))
ggsave(theoph_pop_est_plot, file = 'theoph_pop_est_plot.svg', units = 'mm', height = 150, width = 130)