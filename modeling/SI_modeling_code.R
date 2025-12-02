rm(list = ls())

setwd("C:/Users/rpolivie/OneDrive - Syracuse University/Tao DATA") 

getwd()
list.files()
install.packages(c('lme4', 'lmerTest', 'broom.mixed', 'performance', 'car', 'ggplot2'))
install.packages('reformulas')

library(reformulas)
library(tidyverse)
library(lme4)         
library(lmerTest)     
library(broom.mixed)  
library(performance)  
library(car)         
library(ggplot2)
library(gridExtra)
library(broom)
library(FSA)        
library(ggsignif)  


PANetworkAnalysis <- read.table("BMI_Network_CooccurrenceAnalysis.txt", 
                                header = TRUE, 
                                sep = "\t", 
                                fill = TRUE, 
                                quote = "",
                                stringsAsFactors = FALSE)


df <- read.csv("BugModelingData.csv")

#idiot check to make sure all the numbers are still the same
table(df$Small)  #  0-1839, 1-4987 
table(df$has_aml)#  0-5125, 1-1701
table(df$season) #  fall-860, spring-5966 
table(df$US_L3NAME)

df_spring <- df %>% filter(season == "spring")  %>%
  mutate(
    DEV_z = scale(DEV)[,1],
    uncon_density_z = scale(unconventional_density)[,1],
    conv_density_z = scale(conventional_density)[,1]
  ) %>%  mutate(across(c(CG, FC, PR, SC, SH), ~ . * 100))

#-------------------------------------------------------------------------------

#OGD intensity model
numeric_vars <- df_spring %>% select(where(is.numeric)) %>% select(-all_of(c('DEV_z','uncon_density_z','conv_density_z')))

# Create histograms for each numeric variable
hist_list <- lapply(names(numeric_vars), function(var) {
  ggplot(df_spring, aes_string(x = var)) +
    geom_histogram(fill = "#2C3E50", color = "white", bins = 30) +
    theme_minimal() +
    labs(title = paste("Histogram of", var), x = var, y = "Frequency") +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
})

# Arrange histograms in a grid (adjust ncol as needed)
do.call(grid.arrange, c(hist_list, ncol = 3))

#ogd presence model
numeric_vars <- PANetworkAnalysis %>% select(where(is.numeric))
# Create histograms for each numeric variable
hist_list <- lapply(names(numeric_vars), function(var) {
  ggplot(PANetworkAnalysis, aes_string(x = var)) +
    geom_histogram(fill = "#2C3E50", color = "white", bins = 30) +
    theme_minimal() +
    labs(title = paste("Histogram of", var), x = var, y = "Frequency") +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
})

# Arrange histograms in a grid (adjust ncol as needed)
do.call(grid.arrange, c(hist_list, ncol = 3))


#-------------------------------------------------------------------------------
N_model <- lm(C_pos ~ N_pos, data = PANetworkAnalysis)
C_model <- lm(S_pos ~ C_pos, data = PANetworkAnalysis)
summary(N_model)
summary(C_model)

NetworkModel <- lm(MeanIBI ~  M_pos + N_pos + C_pos + S_pos, data = PANetworkAnalysis) 
summary(NetworkModel) #explains 42% of variance in ibi score

PANetworkAnalysis$PredictedIBI <- predict(NetworkModel)

# Add predicted values to the dataframe
PANetworkAnalysis$PredictedIBI <- predict(NetworkModel)

# Extract model statistics
model_summary <- summary(NetworkModel)
r_squared <- round(model_summary$r.squared, 2)
p_value <- round(glance(NetworkModel)$p.value, 4)

# Create scatter plot with best fit line and annotations
ggplot(PANetworkAnalysis, aes(x = MeanIBI, y = PredictedIBI)) +
  geom_point(color = "#1F77B4", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "#FF7F0E", linetype = "dashed") +
  annotate("text", x = min(PANetworkAnalysis$MeanIBI), y = max(PANetworkAnalysis$PredictedIBI),
           label = paste("R² =", r_squared, "\nP =", p_value),
           hjust = 0, vjust = 1, size = 5, fontface = "bold") +
  theme_minimal() +
  labs(
    title = "Predicted vs Actual IBI",
    x = "Actual IBI",
    y = "Predicted IBI"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  )
#42 % of variance


##------------------------------------------------------------------------------

df_small <- df %>% filter(Small == 1)  %>%
  mutate(
    DEV_z = scale(DEV)[,1],
    uncon_density_z = scale(unconventional_density)[,1],
    conv_density_z = scale(conventional_density)[,1]
  )

df_small <- df_small %>% mutate(
  Richness_z = scale(Richness)[,1],
  SHANdivers_z = scale(SHANdivers)[,1],
  richEPT_z = scale(richEPT)[,1],
  IBI_z = scale(IBI)[,1]
)

#predicting taxonomic metrics
response_vars <- c("Richness_z", "SHANdivers_z", "richEPT_z", "IBI_z")

model_results <- list()

for (resp in response_vars) {
  
  #some formatting so i can tell which results go to what
  cat("\n=========================\n")
  cat("Response:", resp, "\n")
  cat("=========================\n")
  
  #formula- dev land cover, uogd density and cogd density as fixed effects, ecoregion (lv3) and amd_presence as random effects
  formula <- as.formula(paste0(
    resp, " ~ DEV_z + uncon_density_z + conv_density_z + (1 |US_L3NAME) + (1 |has_aml)"
  ))
  
  #REML false for model comparison
  model <- lmer(formula, data = df_small, REML = FALSE)
  
  #check model residuals
  x <- check_model(model)
  
  #check vif (no multicolinearity)
  vif_model <- lm(as.formula(paste0(resp, " ~ DEV_z + uncon_density_z + conv_density_z")),
                  data = df_small)
  vif_vals <- vif(vif_model)
  
  #show fixed effects with p-values
  fixed <- summary(model)$coefficients
  fixed <- cbind(fixed, `Pr(>|t|)` = summary(model)$coefficients[, "Pr(>|t|)"])
  
  #show random effects
  rand_var <- as.data.frame(VarCorr(model))[, c("grp", "vcov")]
  colnames(rand_var) <- c("Group", "Variance")
  
  #get r2 values
  r2_vals <- r2_nakagawa(model)
  
  #summary
  print(x)
  print(fixed)
  cat("\n-- Random Effect Variance --\n")
  print(rand_var)
  cat("\n-- R² Values --\n")
  print(r2_vals)
  cat("\n-- VIF --\n")
  print(vif_vals)
  
  model_results[[resp]] <- list(
    model = model,
    fixed = fixed,
    rand_var = rand_var,
    r2 = r2_vals,
    vif = vif_vals
  )
}
##################### plotting results

library(broom.mixed)
library(dplyr)
library(ggplot2)

#use saved model results
coef_df <- bind_rows(
  lapply(names(model_results), function(name) {
    broom.mixed::tidy(model_results[[name]]$model, effects = "fixed") %>%
      mutate(model = name)
  })
)

#just use fixed effects
coef_df <- coef_df %>% filter(term %in% c("DEV_z","uncon_density_z","conv_density_z"))

#add significance stars
coef_df <- coef_df %>%
  mutate(
    model = factor(model, levels = c("Richness_z", "richEPT_z", "SHANdivers_z", "IBI_z")),
    term = factor(term, levels = c("DEV_z", "conv_density_z", "uncon_density_z")),
    # Add significance stars
    p.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    # y-position for stars
    y_star = estimate + 1.96 * std.error + 0.05 * max(abs(estimate + 1.96 * std.error))
  )
#plot each model result and sd 
ggplot(coef_df, aes(
  x = term,
  y = estimate,
  ymin = estimate - 1.96 * std.error,
  ymax = estimate + 1.96 * std.error
)) +
  geom_point(color = "steelblue", size = 5,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                    ymax = estimate + 1.96 * std.error),
                width = 0.2, color = "steelblue",
                position = position_dodge(width = 0.4)) +
  geom_text(aes(y = y_star, label = p.signif),
            position = position_dodge(width = 0.4),
            vjust = 0, size = 6, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~model, nrow = 1, scales = "free_x") +
  labs(
    y = "Estimate ± 95% CI",
    x = "Predictor",
    title = "Effect sizes across models for wadeable streams only"
  ) +
  theme_bw(base_size = 20) +   
  theme(
    panel.grid = element_blank(),  
    strip.text = element_text(size = 11, face = "bold"),
    panel.spacing = unit(0.5, "lines"),   # restore spacing between facets
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", linewidth = 0.75, fill = NA) 
  )

##############predicting functional proprtions

#predicting funtional metrics
response_vars <- c("CG", "FC", "PR", "SC", "SH")

model_results <- list()

for (resp in response_vars) {
  
  #some formatting so i can tell which results go to what
  cat("\n=========================\n")
  cat("Response:", resp, "\n")
  cat("=========================\n")
  
  #formula- dev land cover, uogd density and cogd density as fixed effects, ecoregion (lv3) and amd_presence as random effects
  formula <- as.formula(paste0(
    resp, " ~ DEV_z + uncon_density_z + conv_density_z + (1 |US_L3NAME) + (1 |has_aml)"
  ))
  
  #REML false for model comparison
  model <- lmer(formula, data = df_small, REML = FALSE)
  
  #check model residuals
  x <- check_model(model)
  
  #check vif (no multicolinearity)
  vif_model <- lm(as.formula(paste0(resp, " ~ DEV_z + uncon_density_z + conv_density_z")),
                  data = df_small)
  vif_vals <- vif(vif_model)
  
  #show fixed effects with p-values
  fixed <- summary(model)$coefficients
  fixed <- cbind(fixed, `Pr(>|t|)` = summary(model)$coefficients[, "Pr(>|t|)"])
  
  #show random effects
  rand_var <- as.data.frame(VarCorr(model))[, c("grp", "vcov")]
  colnames(rand_var) <- c("Group", "Variance")
  
  #get r2 values
  r2_vals <- r2_nakagawa(model)
  
  #summary
  print(x)
  print(fixed)
  cat("\n-- Random Effect Variance --\n")
  print(rand_var)
  cat("\n-- R² Values --\n")
  print(r2_vals)
  cat("\n-- VIF --\n")
  print(vif_vals)
  
  model_results[[resp]] <- list(
    model = model,
    fixed = fixed,
    rand_var = rand_var,
    r2 = r2_vals,
    vif = vif_vals
  )
}
##################### plotting functional results

library(broom.mixed)
library(dplyr)
library(ggplot2)

#use saved model results
coef_df <- bind_rows(
  lapply(names(model_results), function(name) {
    broom.mixed::tidy(model_results[[name]]$model, effects = "fixed") %>%
      mutate(model = name)
  })
)

#just use fixed effects
coef_df <- coef_df %>% filter(term %in% c("DEV_z","uncon_density_z","conv_density_z"))

#add significance stars
coef_df <- coef_df %>%
  mutate(
    model = factor(model, levels = c("CG", "FC", "PR", "SC", "SH")),
    term = factor(term, levels = c("DEV_z", "conv_density_z", "uncon_density_z")),
    # Add significance stars
    p.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    # y-position for stars
    y_star = estimate + 1.96 * std.error + 0.05 * max(abs(estimate + 1.96 * std.error))
  )
#plot each model result and sd 
ggplot(coef_df, aes(
  x = term,
  y = estimate,
  ymin = estimate - 1.96 * std.error,
  ymax = estimate + 1.96 * std.error
)) +
  geom_point(color = "steelblue", size = 5,
             position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                    ymax = estimate + 1.96 * std.error),
                width = 0.2, color = "steelblue",
                position = position_dodge(width = 0.4)) +
  geom_text(aes(y = y_star, label = p.signif),
            position = position_dodge(width = 0.4),
            vjust = 0, size = 6, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~model, nrow = 1, scales = "free_x") +
  labs(
    y = "Estimate ± 95% CI",
    x = "Predictor",
    title = "Effect sizes across models for wadeable streams only"
  ) +
  theme_bw(base_size = 20) +   
  theme(
    panel.grid = element_blank(),  
    strip.text = element_text(size = 11, face = "bold"),
    panel.spacing = unit(0.5, "lines"),   # restore spacing between facets
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", linewidth = 0.75, fill = NA) 
  )

#-------------------------------------------------------------------------------


PANetworkAnalysis <- PANetworkAnalysis %>%
  mutate(U_C_group = factor(U_C_group, levels = c("U0_C0", "U1_C0", "U0_C1", "U1_C1")))

# Kruskal-Wallis test
kruskal_res <- kruskal.test(MeanDEV ~ U_C_group, data = PANetworkAnalysis)
print(kruskal_res)

# Dunn's post-hoc test
dunn_res <- dunnTest(MeanDEV ~ U_C_group, data = PANetworkAnalysis, method="bonferroni")
print(dunn_res)

# Extract significant comparisons
sig_pairs <- dunn_res$res %>%
  filter(P.adj < 0.05) %>%
  select(Comparison) %>%
  mutate(Comparison = strsplit(as.character(Comparison), " - ")) %>%
  tidyr::unnest_wider(Comparison, names_sep = "_")

colnames(sig_pairs) <- c("group1", "group2")
# Convert to list of pairs for geom_signif
sig_list <- split(sig_pairs, seq(nrow(sig_pairs)))
sig_list <- lapply(sig_list, function(x) c(x$group1, x$group2))

ggplot(PANetworkAnalysis, aes(x = U_C_group, y = MeanDEV, fill = U_C_group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  labs(x = "Group (U_bool_C_bool)", y = "MeanDEV") +
  geom_signif(
    comparisons = sig_list,
    map_signif_level = TRUE,
    step_increase = 0.1
  )

#-------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(FSA)        # for dunnTest()
library(ggsignif)
library(tidyr)
library(purrr)
library(patchwork)
library(RColorBrewer)

# Ensure HUC8 is a factor
PANetworkAnalysis <- PANetworkAnalysis %>%
  mutate(HUC8 = as.factor(HUC8))

# --- Variables to analyze ---
vars <- c("N_pos", "C_pos", "S_pos", "M_pos")

# --- Loop through variables and create plots ---
plots <- map(vars, function(var_name) {
  
  # Kruskal–Wallis test
  formula_kw <- as.formula(paste(var_name, "~ HUC8"))
  kw <- kruskal.test(formula_kw, data = PANetworkAnalysis)
  cat("\n---------------------------------------------\n")
  cat("Kruskal-Wallis test for", var_name, "\n")
  print(kw)
  
  # Dunn’s post-hoc test
  dunn_res <- dunnTest(formula_kw, data = PANetworkAnalysis, method = "bonferroni")
  
  # Extract significant comparisons
  sig_pairs <- dunn_res$res %>%
    filter(P.adj < 0.05) %>%
    select(Comparison) %>%
    mutate(Comparison = strsplit(as.character(Comparison), " - ")) %>%
    unnest_wider(Comparison, names_sep = "_")
  
  # Prepare list of comparisons
  sig_list <- if (nrow(sig_pairs) > 0) {
    split(sig_pairs, seq(nrow(sig_pairs))) %>%
      lapply(function(x) c(x$Comparison_1, x$Comparison_2))
  } else {
    list()
  }
  
  # Compute y-limits (auto-scale for significance lines)
  y_max <- max(PANetworkAnalysis[[var_name]], na.rm = TRUE)
  y_lim <- y_max * (1 + 0.25 + 0.05 * length(sig_list))
  
  # Plot
  p <- ggplot(PANetworkAnalysis, aes(x = HUC8, y = .data[[var_name]], fill = HUC8)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal(base_size = 13) +
    labs(
      x = "HUC8 Group",
      y = var_name,
      title = paste("Distribution of", var_name, "by HUC8 Group")
    ) +
    coord_cartesian(ylim = c(NA, y_lim)) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Add significance annotations if there are significant differences
  if (length(sig_list) > 0) {
    p <- p + geom_signif(
      comparisons = sig_list,
      map_signif_level = TRUE,
      step_increase = 0.1
    )
  }
  
  return(p)
})

# Combine all plots in a 2x2 layout
final_plot <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
final_plot