# PACKAGES 

library(readr)
library(cmdstanr)
library(posterior)
library(ggplot2)
library(loo)
library(bridgesampling)

# DATA READING FUNCTION

find_data <- function(path){
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
    col_names = c("Anno", "GEN", "FEB", "MAR", "APR", "MAG", "GIU", 
                  "LUG", "AGO", "SET", "OTT", "NOV", "DIC", "Medio"),
    skip = index_1994_min - 1,
    n_max = n_row,
    col_types = cols(.default = col_double()),
    locale = locale(encoding = "latin1", decimal_mark = "."),
    na = c("", "NA", ">>") 
  )
  
  data_mean <- read_delim(
    path, 
    delim = ";",
    col_names = c("Anno", "GEN", "FEB", "MAR", "APR", "MAG", "GIU", 
                  "LUG", "AGO", "SET", "OTT", "NOV", "DIC", "Medio"),
    skip = index_1994_mean - 1,
    n_max = n_row,
    col_types = cols(.default = col_double()),
    locale = locale(encoding = "latin1", decimal_mark = "."),
    na = c("", "NA", ">>") 
  )
  
  data_max <- read_delim(
    path, 
    delim = ";",
    col_names = c("Anno", "GEN", "FEB", "MAR", "APR", "MAG", "GIU", 
                  "LUG", "AGO", "SET", "OTT", "NOV", "DIC", "Medio"),
    skip = index_1994_max - 1,
    n_max = n_row,
    col_types = cols(.default = col_double()),
    locale = locale(encoding = "latin1", decimal_mark = "."),
    na = c("", "NA", ">>") 
  )
  
  return(list(data_min, data_mean, data_max))
}

# DATA PLOTTING FUNCTION

plot_hist <- function(temp_min, temp_mean, temp_max, month = "", city = "") {
  year <- if (city == "Oderzo") 1994:2023 else 1994:2024
  t_min <- min(c(temp_min, temp_mean, temp_max), na.rm = TRUE)
  t_max <- max(c(temp_min, temp_mean, temp_max), na.rm = TRUE)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par)) 
  par(mar = c(5, 4, 4, 7) + 0.1)
  
  if (month == "Mean") {
    title <- paste("Average Temperature -", city) 
  } else {
    title <- paste("Temperature in", month, "-", city)
  }
  
  plot(year, temp_max,
       type = "n",
       ylim = c(-5, 35),
       main = title,
       xlab = "Year",
       ylab = "Temperature (°C)")
  
  grid(nx = NULL,
       ny = NULL,
       col = "lightgray",
       lty = "dotted",
       lwd = 1)
  
  lines(year, temp_max, 
        type = "h", 
        lwd = 10,
        lend = 1,
        pch = 16, 
        col = "red", border = "white")
  
  lines(year, temp_mean, 
        type = "h", 
        lwd = 10,
        lend = 1,
        pch = 16, 
        col = "darkorange", border = "white")
  
  lines(year, temp_min, 
        type = "h", 
        lwd = 10,
        lend = 1,
        pch = 16, 
        col = "blue", border = "white")
  
  abline(h = 0, col = "gray50", lwd = 1, lty = 2)
  
  legend("topright", 
         legend = c("Max", "Mean", "Min"),
         fill = c("red", "darkorange", "blue"),
         border = NA,
         bty = "n",
         xpd = TRUE,
         inset = c(-0.22, 0))
}

# SUBPLOT FIT FUNCTION

plot_hist_all <- function(data, city = "") {
  
  months <- c("GEN", "FEB", "MAR", "APR", "MAG", "GIU",
              "LUG", "AGO", "SET", "OTT", "NOV", "DIC")
  
  month_name <- c("January", "February", "March", "April", "May", "June",
                  "July", "August", "September", "October", "November", "December")
  
  x <- if (city == "Oderzo") 1994:2023 else 1994:2024
  
  all_min  <- unlist(lapply(data[[1]][months], identity))
  all_mean <- unlist(lapply(data[[2]][months], identity))
  all_max  <- unlist(lapply(data[[3]][months], identity))
  
  y_lim <- range(c(all_min, all_mean, all_max), na.rm = TRUE)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(4, 3),
      mar = c(4, 4, 3, 1),
      oma = c(5, 5, 5, 1),
      
      cex.axis = 1.4,   
      font.axis = 2,    
      
      cex.lab = 1.5,    
      font.lab = 2,    
      
      mgp = c(3, 1.2, 0),
      tck = 0.03)
  
  for (i in seq_along(months)) {
    
    m <- months[i]
    
    temp_min  <- data[[1]][[m]]
    temp_mean <- data[[2]][[m]]
    temp_max  <- data[[3]][[m]]
    
    plot(x, temp_max,
         type = "n",
         ylim = y_lim,
         xlim = range(x),
         main = month_name[i],
         cex.main = 1.9,
         font.main = 2,
         xlab = "",
         ylab = "")
    
    grid(col = "lightgray", lty = "dotted")
    
    lines(x, temp_max,  type = "h", lwd = 9, col = "red", lend = 1)
    lines(x, temp_mean, type = "h", lwd = 9, col = "darkorange", lend = 1)
    lines(x, temp_min,  type = "h", lwd = 9, col = "blue", lend = 1)
    
    abline(h = 0, col = "gray50", lty = 2)
    
    if (i == 1) {
      legend("topright",
             legend = c("Max", "Mean", "Min"),
             fill = c("red", "darkorange", "blue"),
             border = NA,
             bty = "n",
             cex = 1.5,
             text.font = 2,
             pt.cex = 2)
    }
  }
  
  mtext(paste("Min-Mean-Max Temperature -", city),
        outer = TRUE,
        line = 2,
        cex = 1.8,
        font = 2)
  
  mtext("Year",
        side = 1,
        outer = TRUE,
        line = 2.5,
        cex = 1.6,
        font = 2)
  
  mtext("Temperature (°C)",
        side = 2,
        outer = TRUE,
        line = 2.5,
        cex = 1.6,
        font = 2)
}

# CORRELATIONS FUNCTION

plot_2_correlations <- function(my_draws_df, params_list, axes_names = c(expression(T[0]), expression(dot(T)))) {
  cm <- cor(as.matrix(my_draws_df[, params_list])) 
  p <- ggplot(my_draws_df, aes(x = .data[[params_list[1]]], y = .data[[params_list[2]]])) +
    geom_point(alpha = 0.15, size = 0.8, color = "black") +
    coord_fixed(ratio = 1) + 
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

# POSTERIOR PLOT

plot_posterior <- function(my_draws_df, my_x_stdd, my_y_data,
                           my_x_grid = NULL, n_grid = N, n_pred_grid = N_fut_steps,
                           orig_mean = years_mean, orig_sd = years_sd_ceiling, my_step = step,
                           xb_l = 5, xl_l = 5, yg_b = 0.5) {
  if (is.null(my_x_grid)) { 
    x_grid <- seq(min(my_x_stdd), max(my_x_stdd) + my_step * n_pred_grid, length.out = n_grid + n_pred_grid)
  } else {
    x_grid <- my_x_grid
  }
  
  S <- nrow(my_draws_df)
  T_mat <- matrix(NA, nrow = S, ncol = length(x_grid))
  
  if (!"T_dot" %in% colnames(my_draws_df)) {
    T_dot_vals <- rep(0, S) 
  } else {
    T_dot_vals <- my_draws_df$T_dot
  }
  
  for (i in seq_along(x_grid)) {
    T_mat[, i] <- my_draws_df$T0 + T_dot_vals * x_grid[i]
  }
  
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
  
  p <- ggplot() +
    geom_point(aes(x = my_x_stdd, y = my_y_data, color = "Observations")) +
    geom_line(data = plot_df, aes(x = x_stdd, y = mean, color = "Posterior mean"), linewidth = 1) +
    geom_ribbon(data = plot_df, aes(x = x_stdd, ymin = lo, ymax = hi, fill = "95% CI"), alpha = 0.3) +
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
    ) + 
    
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

# GRAPHS: Traces, ACFs, densities

plot_various_graphs <- function(my_fit, params_list) {
  draws_array <- as.array(my_fit$draws())
  n_iters <- dim(draws_array)[1]
  n_chains <- dim(draws_array)[2]
  
  df_list <- list()
  idx <- 1 
  for (p in params_list) { 
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
  df_long <- do.call(rbind, df_list) 
  
  label_map <- c("T0" = "T[0]", "T_dot" = "dot(T)", "sigma" = "sigma")
  df_long$Parameter <- factor(df_long$Parameter,
                              levels = params_list,
                              labels = label_map[params_list]
  )
  
  trace_plot <- ggplot(df_long, aes(x = Iteration, y = Value, color = Chain)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~Parameter, scales = "free_y", ncol = 1, labeller = label_parsed) + 
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
  
  hist_plot <- ggplot(df_long, aes(x = Value, fill = Chain, color = Chain)) +
    geom_histogram(aes(y = after_stat(count / sum(count))), alpha = 0.4, bins = 30, position = "identity") +
    facet_wrap(~Parameter, scales = "free", ncol = 1, labeller = label_parsed) +
    labs(x = "Parameter value", y = "Relative frequency") +
    scale_color_manual(values = c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73", "4" = "#CC79A7")) + 
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
                             levels = params_list 
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

# Linear Model Log-Posterior for BF comparison (via BridgeSampling)

log_posterior_linear <- function(samples.row, data) {
  T0 <- samples.row["T0"]
  T_dot <- samples.row["T_dot"]
  sigma <- samples.row["sigma"]
  
  # Log priors + Log likelihood
  dnorm(T0, 0, 20, log = TRUE) +
    dnorm(T_dot, 0, 20, log = TRUE) +
    dnorm(sigma, 0, 10, log = TRUE) - pnorm(0, 0, 10, lower.tail = FALSE, log.p = TRUE) +
    sum(dnorm(data$T_obs, T0 + T_dot * data$time_vec, sigma, log = TRUE))
}

# Constant Model Log-Posterior for BF comparison (via BridgeSampling)

log_posterior_constant <- function(samples.row, data) {
  T0 <- samples.row["T0"]
  sigma <- samples.row["sigma"]
  
  # Log priors + Log likelihood
  dnorm(T0, 0, 20, log = TRUE) +
    dnorm(sigma, 0, 10, log = TRUE) - pnorm(0, 0, 10, lower.tail = FALSE, log.p = TRUE) + 
    sum(dnorm(data$T_obs, T0, sigma, log = TRUE))
}

# Predictive pdf functions

linear_predictive_pdf <- function(chain, t_new) {
  
  t0_samples    <- chain$T0
  t_dot_samples <- chain$T_dot
  sigma_samples <- chain$sigma
  
  mu_components    <- t0_samples + t_dot_samples * t_new
  sigma_components <- sigma_samples
  
  # Return PDF function
  function(x) {
    sapply(x, function(val) {
      mean(dnorm(val, mean = mu_components, sd = sigma_components))
    })
  }
}

constant_predictive_pdf <- function(chain, t_new) {
  
  t0_samples    <- chain$T0
  sigma_samples <- chain$sigma
  
  mu_components    <- t0_samples
  sigma_components <- sigma_samples
  
  # Return PDF function
  function(x) {
    sapply(x, function(val) {
      mean(dnorm(val, mean = mu_components, sd = sigma_components))
    })
  }
}

# Pred PDF plot function

plot_predictive_comparison <- function(draws_lin, draws_const, time_future, years_sd, years_mean, temp_range = c(5, 20)) {
  
  layout_matrix <- matrix(c(
    1, 2,  3,  4,  5,  6,
    1, 7,  8,  9, 10, 11
  ), nrow = 2, byrow = TRUE)
  
  layout(layout_matrix, widths = c(0.5, 1, 1, 1, 1, 1))
  
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  legend("center",
         legend = c("Linear: Pred PDF", "Constant: Pred PDF"),
         col = c("darkblue", "red"),
         lwd = 2,
         lty = c(1, 2),
         bty = "n",
         cex = 0.8)
  
  temperature_grid <- seq(temp_range[1], temp_range[2], length.out = 300)
  years_future_real <- (time_future * years_sd) + years_mean
  
  n_steps <- length(time_future) 
  
  for (i in 1:n_steps) {
    par(mar = c(4, 4, 3, 1))
    
    pdf_lin_fun <- linear_predictive_pdf(draws_lin, time_future[i])
    pdf_const_fun <- constant_predictive_pdf(draws_const, time_future[i])
    
    density_lin <- pdf_lin_fun(temperature_grid)
    density_const <- pdf_const_fun(temperature_grid)
    
    y_max <- max(c(density_lin, density_const)) * 1.15
    
    plot(temperature_grid, density_lin, 
         type = "l", 
         lwd = 2, 
         col = "darkblue",
         xlab = "Temp (°C)", 
         ylab = "Density",
         main = paste0("Step ", i, "\n(", round(years_future_real[i]), ")"),
         cex.main = 0.9,
         cex.axis = 0.8,
         ylim = c(0, y_max)) 
    
    lines(temperature_grid, density_const, lwd = 2, col = "red", lty = 2)
    
    grid()
  }
  
  layout(1)
}

# Predictive summaries

predictive_summaries <- function(chain, time_vector, model_type = "linear") {
  
  results_list <- lapply(time_vector, function(t_val) {
    
    if (model_type == "linear") {
      pdf_fun <- linear_predictive_pdf(chain, t_val)
      mu_chain <- chain$T0 + chain$T_dot * t_val
    } else {
      pdf_fun <- constant_predictive_pdf(chain, t_val)
      mu_chain <- chain$T0
    }
    
    exact_mean <- mean(mu_chain)
    exact_sd <- sqrt(mean(chain$sigma^2) + mean((mu_chain - exact_mean)^2)) # variance decomposition for gaussian mixture
    
    temp_grid <- seq(exact_mean - 6 * exact_sd, exact_mean + 6 * exact_sd, length.out = 1200)
    dx <- temp_grid[2] - temp_grid[1]
    
    densities <- pdf_fun(temp_grid)
    densities <- densities / sum(densities * dx) # numerical normalization of the mixture
    
    map_estimate <- temp_grid[which.max(densities)]
    
    cdf <- cumsum(densities * dx) # numerical approx of the cdf
    ci_lower <- temp_grid[which.min(abs(cdf - 0.025))]
    ci_upper <- temp_grid[which.min(abs(cdf - 0.975))]
    
    data.frame(
      Time_Step = t_val,
      MAP = round(map_estimate, 3),
      Mean = round(exact_mean, 3),
      SD = round(exact_sd, 3),
      CI_2.5 = round(ci_lower, 3),
      CI_97.5 = round(ci_upper, 3)
    )
  })
  
  return(do.call(rbind, results_list))
}