# Load required libraries
library(tidyverse)
library(data.table)
library(ggplot2)
library(magrittr)
library(dplyr)
library(xtable)
library(gridExtra)
library(reshape2)

# Source functions from external file
source('functions.R')

# Function to extract parameters from the filename
extract_settings <- function(filename) {
  basename <- basename(filename)
  params <- regmatches(basename, gregexpr("[-0-9]+", basename))[[1]]
  params <- as.numeric(params)
  names(params) <- c("n", "Q", "B", "C", "trial")
  return(params)
}

# Read all files
files <- list.files("simulation_results", pattern = "*.rds", full.names = TRUE)
nsim <- length(files)

param_names <- c('dem', 'dcm', 'deh', 'dmv')
stat_names <- c('est', 'se', 'lower', 'upper')

df_names <- outer(param_names, stat_names, FUN = function(x, y) paste(x, y, sep = "_")) %>% as.vector()
mat <- matrix(NA, nrow = nsim, ncol = length(df_names))
df <- data.frame(mat)
df <- df %>% setNames(df_names)

settings_list <- lapply(files, extract_settings)
results_list <- lapply(files, readRDS)

# Calculate 'true' values for each parameter combination
true_values <- lapply(settings_list, function(setting) {
  list(
    dcm = calc_dcm(setting[['Q']]/10, setting[['B']], setting[['C']])$dcm,
    deh = calc_deh(setting[['Q']]/10, setting[['B']], setting[['C']])$deh,
    dem = calc_dem(setting[['Q']]/10, setting[['B']], setting[['C']])$dem,
    dmv = calc_dmv(setting[['Q']]/10, setting[['B']], setting[['C']])$dmv
  )
})

# Add columns for settings
df$n <- NA
df$Q <- NA
df$B <- NA
df$C <- NA

# Populate df with estimates and other data
for(p in param_names){
    df[[paste0(p, '_SIMCHECK')]] <- rep(NA, nsim)
    
    for(i in 1:nsim){
        df$n[i] <- settings_list[[i]]['n']
        df$Q[i] <- settings_list[[i]]['Q']
        df$B[i] <- settings_list[[i]]['B']
        df$C[i] <- settings_list[[i]]['C']
        df[[paste0(p, '_SIMCHECK')]][i] <- results_list[[i]][[paste0('sim_', p)]]
        df[[paste0(p, '_est')]][i] <- results_list[[i]][['bindecomp']][[p]]$est
        df[[paste0(p, '_lower')]][i] <- results_list[[i]][['bindecomp']][[p]]$lower
        df[[paste0(p, '_upper')]][i] <- results_list[[i]][['bindecomp']][[p]]$upper
        df[[paste0(p, '_se')]][i] <- results_list[[i]][['bindecomp']][[p]]$se
    }
}

# Initialize summary dataframe
summary <- data.frame(matrix(NA, nrow = 0, ncol = 16))
colnames(summary) <- c("n", "Q", "B", "C", "dcm", "deh", "dem", "dmv", 
                       "MSE_dcm", "MSE_deh", "MSE_dem", "MSE_dmv", 
                       "Cov_dcm", "Cov_deh", "Cov_dem", "Cov_dmv")

# Get unique settings
settings <- unique(df[, c("n", "Q", "B", "C")])

# Compute MSE and coverage probability
for(j in 1:nrow(settings)){
    setting <- settings[j,]
    n_val <- setting[["n"]]; Q_val <- setting[["Q"]]; B_val <- setting[["B"]]; C_val <- setting[["C"]]
    
    idx <- which(df$n == n_val & df$Q == Q_val & df$B == B_val & df$C == C_val)
    
    df_filt <- df[idx, ]
    
    true_val_list <- true_values[[idx[1]]]
    
    summary_row <- list(n = n_val, Q = Q_val/10, B = B_val, C = C_val)
    
    for(p in param_names){
        true_val <- true_val_list[[p]]
        summary_row[[p]] <- true_val

        mse <- mean((df_filt[[paste0(p, '_est')]] - true_val)^2)
        summary_row[[paste0('MSE_', p)]] <- mse

        coverage <- mean(true_val >= df_filt[[paste0(p, '_lower')]] & true_val <= df_filt[[paste0(p, '_upper')]])
        summary_row[[paste0('Cov_', p)]] <- coverage
    }
    
    summary <- rbind(summary, summary_row)
}

# Order summary by n
summary <- summary[order(summary$n),]

# Create Scenario variable
summary$Scenario <- as.numeric(factor(paste0("Q", summary$Q, "B", summary$B, "C", summary$C)))
df$Scenario <- as.numeric(factor(paste0("Q", df$Q, "B", df$B, "C", df$C)))

# Save summary
saveRDS(summary, 'summary.rds')

# Compute n*MSE
summary$nMSE_dcm <- summary$n * summary$MSE_dcm
summary$nMSE_deh <- summary$n * summary$MSE_deh
summary$nMSE_dem <- summary$n * summary$MSE_dem
summary$nMSE_dmv <- summary$n * summary$MSE_dmv

# Reshape data for plotting n*MSE
nMSE_data <- summary %>% select(n, Q, B, C, Scenario, nMSE_dcm, nMSE_deh, nMSE_dem, nMSE_dmv)
nMSE_data_long <- melt(nMSE_data, id.vars = c('n', 'Q', 'B', 'C', 'Scenario'),
                       variable.name = 'Parameter', value.name = 'nMSE')
nMSE_data_long$Parameter <- sub('nMSE_', '', nMSE_data_long$Parameter)

# Reshape data for plotting coverage probability
Cov_data <- summary %>% select(n, Q, B, C, Scenario, Cov_dcm, Cov_deh, Cov_dem, Cov_dmv)
Cov_data_long <- melt(Cov_data, id.vars = c('n', 'Q', 'B', 'C', 'Scenario'),
                      variable.name = 'Parameter', value.name = 'Coverage')
Cov_data_long$Parameter <- sub('Cov_', '', Cov_data_long$Parameter)

# Merge data for plotting
plot_data <- merge(nMSE_data_long, Cov_data_long, by = c('n', 'Q', 'B', 'C', 'Scenario', 'Parameter'))

# Modify the Scenario labels 
scenario_labels <- paste0("Scenario ", sort(unique(plot_data$Scenario)))
names(scenario_labels) <- sort(unique(plot_data$Scenario))

# Define mapping from parameter names to expressions
parameter_labels <- c(
  dem = expression(delta[EM]),
  dcm = expression(delta[CM]),
  deh = expression(delta[EH]),
  dmv = expression(delta[MV])
)

# Define plotting theme
my_theme <- theme_bw(base_size = 14) +
    theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
    )

# Plot n*MSE
nMSE_plot <- ggplot(plot_data, aes(x = n, y = nMSE, color = Parameter)) +
    geom_point(size = 4) +
    geom_line(size = 1.2) +
    facet_wrap(~ Scenario, scales = "free", ncol = 2, labeller = labeller(Scenario = scenario_labels)) +
    labs(title = "n*MSE vs Sample Size",
         x = "Sample Size (n)",
         y = expression(n * MSE)) +
    scale_x_log10() +
    scale_color_manual(
      name = "Parameter",
      values = c("dem" = "red", "dcm" = "blue", "deh" = "green", "dmv" = "purple"),
      labels = parameter_labels
    ) +
    my_theme

# Plot Coverage Probability 
Coverage_plot <- ggplot(plot_data, aes(x = n, y = Coverage, color = Parameter)) +
    geom_point(size = 4) +
    geom_line(size = 1.2) +
    facet_wrap(~ Scenario, scales = "free", ncol = 2, labeller = labeller(Scenario = scenario_labels)) +
    labs(title = "Coverage Probability vs Sample Size",
         x = "Sample Size (n)",
         y = "Coverage Probability") +
    scale_x_log10() +
    scale_color_manual(
      name = "Parameter",
      values = c("dem" = "red", "dcm" = "blue", "deh" = "green", "dmv" = "purple"),
      labels = parameter_labels
    ) +
    my_theme

# Save the plots
ggsave("nMSE_plot.png", plot = nMSE_plot, width = 12, height = 8)
ggsave("Coverage_plot.png", plot = Coverage_plot, width = 12, height = 8)

# Create LaTeX table for scenarios and parameter values
scenario_table <- summary %>%
  select(Scenario, Q, B, C, dcm, deh, dem, dmv) %>%
  distinct()

latex_scenario_table <- xtable(scenario_table, 
                               caption = "Scenarios and Corresponding Parameter Values",
                               label = "tab:scenarios")

# Save LaTeX table
sink("scenario_table.txt")
print(latex_scenario_table, include.rownames = FALSE, floating = TRUE,
      table.placement = "htbp",
      caption.placement = "top",
      hline.after = c(-1, 0, nrow(scenario_table)),
      sanitize.text.function = function(x){x})
sink()

# Generate caption for the panel of figures
figure_caption <- "Simulation results showing \\( n \\times \\text{Mean Squared Error (MSE)} \\) and Coverage Probability for the parameters \\( \\delta_{CM} \\), \\( \\delta_{EM} \\), \\( \\delta_{MV} \\), and \\( \\delta_{EH} \\) across different sample sizes and scenarios. Each scenario represents a unique combination of Q, B, and C values, numbered from 1 to 8."

# Save the caption to a file
writeLines(figure_caption, "figure_caption.txt")

