# for correct compilation of stan files (make.exe)
#Sys.setenv(PATH = paste("C:/rtools44/usr/bin", Sys.getenv("PATH"), sep = ";"))

library(readr)
library(cmdstanr)
library(posterior)
library(ggplot2)
library(loo)

####################################################################################
# FROM file arpav.R
find_data <- function(path) {
  row <- read_lines(path, locale = locale(encoding = "latin1", decimal_mark = "."))

  row_min <- grep("Temperatura", row)[1]
  index_1994_min <- row_min + grep("^1994;", row[row_min:length(row)])[1] - 1
  index_2024_min <- row_min + grep("^2024;", row[row_min:length(row)])[1] - 1

  row_mean <- grep("Temperatura", row)[2]
  index_1994_mean <- row_mean + grep("^1994;", row[row_mean:length(row)])[1] - 1
  index_2024_mean <- row_mean + grep("^2024;", row[row_mean:length(row)])[1] - 1

  row_max <- grep("Temperatura", row)[3]
  index_1994_max <- row_max + grep("^1994;", row[row_max:length(row)])[1] - 1
  index_2024_max <- row_max + grep("^2024;", row[row_max:length(row)])[1] - 1
  n_row <- index_2024_max - index_1994_max + 1

  data_min <- read_delim(
    path,
    delim = ";",
    col_names = c(
      "Anno", "GEN", "FEB", "MAR", "APR", "MAG", "GIU",
      "LUG", "AGO", "SET", "OTT", "NOV", "DIC", "Medio"
    ),
    skip = index_1994_min - 1,
    n_max = n_row,
    col_types = cols(.default = col_double()),
    locale = locale(encoding = "latin1", decimal_mark = "."),
    na = c("", "NA", ">>")
  )

  data_mean <- read_delim(
    path,
    delim = ";",
    col_names = c(
      "Anno", "GEN", "FEB", "MAR", "APR", "MAG", "GIU",
      "LUG", "AGO", "SET", "OTT", "NOV", "DIC", "Medio"
    ),
    skip = index_1994_mean - 1,
    n_max = n_row,
    col_types = cols(.default = col_double()),
    locale = locale(encoding = "latin1", decimal_mark = "."),
    na = c("", "NA", ">>")
  )

  data_max <- read_delim(
    path,
    delim = ";",
    col_names = c(
      "Anno", "GEN", "FEB", "MAR", "APR", "MAG", "GIU",
      "LUG", "AGO", "SET", "OTT", "NOV", "DIC", "Medio"
    ),
    skip = index_1994_max - 1,
    n_max = n_row,
    col_types = cols(.default = col_double()),
    locale = locale(encoding = "latin1", decimal_mark = "."),
    na = c("", "NA", ">>")
  )

  return(list(data_min, data_mean, data_max))
}


####################################################################
# TO DO: implement this as function (pick location, mins or means or maxes, and month(s) or yearly avg)

# NOTE: "." assumed project folder
data_montesummano <- find_data("./data/Monte_Summano.csv")
ms_means <- data_montesummano[[2]]
ms_means_yearavg <- ms_means$Medio
cat(ms_means_yearavg)
# 11 9.9 10 10.5 10.2 10.4 11.1 10.7 10.8 12.8 10.8 10.3 11.3 12 11.3 11.8 10.3 11.9 11.7 11.3 12 12.2 11.5 11.6 12.3 12.2 11.9 11.4 12.8 12.5 12.9
ms_m_ya_length <- length(ms_means_yearavg)
cat(ms_m_ya_length)
# 31
years <- ms_means$Anno
cat(years)
# 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024

# standardize x variable (years)
years_mean <- mean(years)
cat(years_mean)
# 2009
years_sd_ceiling <- ceiling(sd(years))
cat(years_sd_ceiling) # cleaner bc 10
# 10
years_stdd <- (years - years_mean) / years_sd_ceiling
N <- ms_m_ya_length
step <- abs(years_stdd[1] - years_stdd[2]) # for standardized data
N_fut_steps <- 10
time_future <- seq(1 + step, by = step, length.out = N_fut_steps) # for standardized data
stan_data_sdd <- list(N = N, time = years_stdd, T = ms_means_yearavg, N_future = N_fut_steps, time_future = time_future)

# function to compile and fit both stan models
build_fit_stan_model <- function(stan_file_path, my_data) {
  my_mod <- cmdstan_model(stan_file_path)
  cat("MODEL COMPILED SUCCESSFULLY\n")
  my_fit <- my_mod$sample(
    data = my_data,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 2000,
    seed = 222 # for reproducibility
  )
  return(list(mod = my_mod, fit = my_fit))
}

res <- build_fit_stan_model("./linear_model.stan", stan_data_sdd) # chain takes about 2 seconds
mod <- res$mod
fit <- res$fit
stan_T_param_names <- c("T0", "T_dot")
fit$summary(stan_T_param_names)
fit$summary("sigma")
fit$cmdstan_diagnose()

# function to plot chain samples on par1-par2 plane (main model only)
plot_2_correlations <- function(my_draws_df, params_list, axes_names = c(expression(T[0]), expression(dot(T)))) {
  cm <- cor(as.matrix(my_draws_df[, params_list])) # correlation matrix
  p <- ggplot(my_draws_df, aes(x = .data[[params_list[1]]], y = .data[[params_list[2]]])) +
    geom_point(alpha = 0.15, size = 0.8, color = "black") +
    coord_fixed(ratio = 1) + # to keep correct proportions (visual clue)
    labs(
      x = axes_names[1],
      y = axes_names[2]
    ) +
    theme_classic(base_size = 12, base_family = "serif") +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10, color = "black")
    )
  return(list(cor = cm, plot = p))
}

draws_df <- as_draws_df(fit$draws()) # draws from posterior: dataframe
lapply(draws_df[, stan_T_param_names], quantile, probs = c(0.025, 0.975)) # quantiles for 95% CI on predicted T

cor_plot <- plot_2_correlations(draws_df, stan_T_param_names)
print(cor_plot$plot)
cor_plot$cor # correlation matrix: abs value in 0,1 (corr = covariance / product of variances)
# TO DO: T0 vs sigma, T_dot vs sigma; or, in 3 dimensions?
cor_plot_sigma_0 <- plot_2_correlations(draws_df, c("T0", "sigma"))
print(cor_plot_sigma_0$plot)
cor_plot_sigma_0$cor
cor_plot_sigma_dot <- plot_2_correlations(draws_df, c("T_dot", "sigma"))
print(cor_plot_sigma_dot$plot)
cor_plot_sigma_dot$cor
# TO DO: all in one?

# function to  plot posterior: mean, 95% CI, comparison to data
plot_posterior <- function(my_draws_df, my_x_stdd, my_y_data,
                           my_x_grid = NULL, n_grid = N, n_pred_grid = N_fut_steps,
                           orig_mean = years_mean, orig_sd = years_sd_ceiling, my_step = step,
                           xb_l = 5, xl_l = 5, yg_b = 0.5) {
  if (is.null(my_x_grid)) { # custom grid optional
    x_grid <- seq(min(my_x_stdd), max(my_x_stdd) + my_step * n_pred_grid, length.out = n_grid + n_pred_grid)
  } else {
    x_grid <- my_x_grid
  }

  # matrix: rows correspond to (T0, T_dot) samples, columns to the x grid
  S <- nrow(my_draws_df)
  T_mat <- matrix(NA, nrow = S, ncol = length(x_grid))

  # to consider only T0 if model is the constant one
  if (!"T_dot" %in% colnames(my_draws_df)) {
    T_dot_vals <- rep(0, S) # zeros
  } else {
    T_dot_vals <- my_draws_df$T_dot
  }

  # predicted T value for each matrix element
  for (i in seq_along(x_grid)) {
    T_mat[, i] <- my_draws_df$T0 + T_dot_vals * x_grid[i]
  }

  # ... get its mean and quantiles; apply to columns ("2")
  T_mean <- apply(T_mat, 2, mean)
  T_lo <- apply(T_mat, 2, quantile, 0.025)
  T_hi <- apply(T_mat, 2, quantile, 0.975)
  T_pred_mat <- matrix(NA, S, length(x_grid))
  for (i in 1:S) {
    T_pred_mat[i, ] <- rnorm(length(x_grid), mean = T_mat[i, ], sd = my_draws_df$sigma[i])
  }
  T_pred_lo <- apply(T_pred_mat, 2, quantile, 0.025)
  T_pred_hi <- apply(T_pred_mat, 2, quantile, 0.975)

  plot_df <- data.frame(x_stdd = x_grid, mean = T_mean, lo = T_lo, hi = T_hi, p_lo = T_pred_lo, p_hi = T_pred_hi)

  # Define x axis limits based on observed data
  # x_min_obs <- min(my_x_stdd)
  # x_max_obs <- max(my_x_stdd)

  # Fine grid positions (every 0.2) for x and y (but restrict grid to observed range)
  # x_grid_pos <- seq(floor(x_min_obs), ceiling(x_max_obs), by = 0.2)
  # x_grid_pos <- x_grid_pos[x_grid_pos >= x_min_obs & x_grid_pos <= x_max_obs]
  # y_grid_pos <- seq(0, 1, by = 0.2)

  # Integer breaks for x-axis labels (only those inside the observed range)
  # integer_breaks <- seq(floor(x_min_obs), ceiling(x_max_obs), by = 1)
  # integer_breaks <- integer_breaks[integer_breaks >= x_min_obs & integer_breaks <= x_max_obs]

  p <- ggplot() +
    # Manual grid lines
    # geom_hline(yintercept = y_grid_pos, color = "grey85", linewidth = 0.3) +
    # geom_vline(xintercept = x_grid_pos, color = "grey85", linewidth = 0.3) +
    # data points
    geom_point(aes(x = my_x_stdd, y = my_y_data, color = "Observations")) +
    # Posterior mean line
    geom_line(data = plot_df, aes(x = x_stdd, y = mean, color = "Posterior mean"), linewidth = 1) +
    # Posterior credible ribbon
    geom_ribbon(data = plot_df, aes(x = x_stdd, ymin = lo, ymax = hi, fill = "95% CI"), alpha = 0.3) +
    # predictive
    # geom_ribbon(data=plot_df, aes(x= x_stdd[n_grid:n_pred_grid], y=p_lo[n_grid:n_pred_grid], ymax=p_hi[n_grid:n_pred_grid], fill="95% CI (pred.)"), alpha=0.3) +
    geom_ribbon(data = plot_df[(n_grid + 1):(n_grid + n_pred_grid), ], aes(x = x_stdd, ymin = p_lo, ymax = p_hi, fill = "95% CI (pred.)"), alpha = 0.3) +
    scale_color_manual(values = c("Observations" = "black", "Posterior mean" = "blue")) +
    scale_fill_manual(values = c("95% CI" = "grey70", "95% CI (pred.)" = "grey40")) +
    labs(x = "Year", y = "Average temperature (°C)") +
    guides(color = guide_legend(title = NULL), fill = guide_legend(title = NULL)) +
    scale_x_continuous(
      breaks = seq(min(x_grid), max(x_grid), length.out = xb_l),
      labels = round(orig_mean + orig_sd * seq(min(x_grid), max(x_grid), length.out = xl_l))
    ) +
    geom_hline(
      yintercept = seq(floor(min(my_y_data)), floor(max(T_pred_hi)), by = yg_b),
      color = "grey85", linewidth = 0.3
    ) + # y grid

    theme_classic(base_size = 12, base_family = "serif") +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10, color = "black"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_blank(),
      legend.text = element_text(size = 10, family = "serif"),
      legend.margin = margin(t = 0, b = 0, unit = "pt"),
      legend.spacing.x = unit(5, "pt")
    )

  return(p)
}
# TO DO: prettier?

plot <- plot_posterior(draws_df, my_x_stdd = stan_data_sdd$time, my_y_data = stan_data_sdd$T)
print(plot)

# function to plot trace plot(s), sampels histogram(s), ACF(s)
plot_various_graphs <- function(my_fit, params_list) {
  draws_array <- as.array(my_fit$draws())
  n_iters <- dim(draws_array)[1]
  n_chains <- dim(draws_array)[2]

  df_list <- list()
  idx <- 1 # index for chains
  for (p in params_list) { # T0, and maybe T_dot
    for (chain in 1:n_chains) {
      vals <- draws_array[, chain, p]
      df_list[[idx]] <- data.frame(
        Parameter = p,
        Chain = factor(chain),
        Iteration = seq_along(vals),
        Value = as.numeric(vals),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  df_long <- do.call(rbind, df_list) # all chains

  label_map <- c("T0" = "T[0]", "T_dot" = "dot(T)", "sigma" = "sigma")
  df_long$Parameter <- factor(df_long$Parameter,
    levels = params_list,
    labels = label_map[params_list]
  )
  # trace plots part
  trace_plot <- ggplot(df_long, aes(x = Iteration, y = Value, color = Chain)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~Parameter, scales = "free_y", ncol = 1, labeller = label_parsed) + # modify n_cols, free scales if necessary
    # label_parsed interprets the strings as plotmath expressions;
    # facet_wrap for 2-dimensional output graphs
    labs(x = "Iteration", y = "Parameter value") +
    scale_color_manual(values = c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73", "4" = "#CC79A7")) +
    theme_classic(base_size = 12, base_family = "serif") +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12, family = "serif"),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )

  # occurences histogram part
  hist_plot <- ggplot(df_long, aes(x = Value, fill = Chain, color = Chain)) +
    # relative frequency: each bar = count in bin / total count across all chains
    geom_histogram(aes(y = after_stat(count / sum(count))), alpha = 0.4, bins = 30, position = "identity") +
    # aesthetics only
    facet_wrap(~Parameter, scales = "free", ncol = 1, labeller = label_parsed) +
    labs(x = "Parameter value", y = "Relative frequency") +
    scale_color_manual(values = c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73", "4" = "#CC79A7")) + # colorblind-friendly palette
    scale_fill_manual(values = c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73", "4" = "#CC79A7")) +
    theme_classic(base_size = 12, base_family = "serif") +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12, family = "serif"),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )

  # ACFs part
  acf_list <- list()
  for (param in params_list) {
    param_draws <- as.vector(my_fit$draws(param))
    acf_obj <- acf(param_draws, plot = FALSE)
    acf_list[[param]] <- data.frame(
      Lag = as.numeric(acf_obj$lag),
      ACF = as.numeric(acf_obj$acf),
      Parameter = param
    )
  }
  acf_df <- do.call(rbind, acf_list)
  acf_df$Parameter <- factor(acf_df$Parameter,
    levels = params_list # ,
    # labels = paste0("beta[", sub("b", "", params_list), "]")
  )
  acf_df$Parameter <- factor(acf_df$Parameter,
    levels = params_list,
    labels = label_map[params_list]
  )
  acf_plot <- ggplot(acf_df, aes(x = Lag, y = ACF)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    geom_segment(aes(xend = Lag, yend = 0), color = "black", linewidth = 0.7) +
    facet_wrap(~Parameter, ncol = 1, labeller = label_parsed) +
    labs(x = "Lag", y = "ACF") +
    theme_classic(base_size = 12, base_family = "serif") +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12, family = "serif"),
      legend.position = "none",
      legend.text = element_text(size = 9, family = "serif")
    )
  return(list(acf = acf_plot, trace = trace_plot, density = hist_plot))
}
graphs <- plot_various_graphs(fit, c("T0", "T_dot", "sigma"))
print(graphs$trace)
print(graphs$density)
print(graphs$acf)

# CONSTANT MODEL
res <- build_fit_stan_model("./constant_model.stan", stan_data_sdd)
c_mod <- res$mod
c_fit <- res$fit
c_fit$summary(c("T0", "sigma"))
c_fit$cmdstan_diagnose()
c_draws_df <- as_draws_df(c_fit$draws()) # draws from posterior
plot_posterior(c_draws_df, my_x_stdd = stan_data_sdd$time, my_y_data = stan_data_sdd$T)
c_graphs <- plot_various_graphs(c_fit, c("T0", "sigma"))
print(c_graphs$trace)
print(c_graphs$density)
print(c_graphs$acf)

# LOO COMPARISON
# extract log-likelihood matrix (draws x N)
log_lik <- fit$draws("log_lik", format = "matrix")
# compute LOO object for both models
loo <- loo(log_lik)
c_log_lik <- c_fit$draws("log_lik", format = "matrix")
c_loo <- loo(c_log_lik)
comp <- loo_compare(loo, c_loo) # compare the two models based on ELPD statistics
print(loo)
print(c_loo)
print(comp) 

# minor issue (bc one observation weird):
# Pareto k diagnostic values:
#             Count Pct.   
# (-Inf, 0.7] 30    96.8%  
#    (0.7, 1]  1     3.2%      <<<<<<<<<<<---------------------
#    (1, Inf)  0     0.0%  