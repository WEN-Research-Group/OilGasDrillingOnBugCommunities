
#this code uses mixed effect modeling to test the association between potential stressors (cogd density, uogd density, and DLC)
#as fixed effects and ecoregions and amd as random effects on 1) taxonomic metrics and 2) ffg proportions

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
  )

#table(df_spring$SMALL)
#table(df_spring$has_aml)
#table(df_spring$season)
#table(df_spring$US_L3NAME)

#turning proportions into percentages
df_spring <- df_spring %>%
  mutate(across(c(CG, FC, PR, SC, SH), ~ . * 100))


df_spring <- df_spring %>% mutate(
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
  model <- lmer(formula, data = df_spring, REML = FALSE)
  
  #check model residuals
  x <- check_model(model)
  
  #check vif (no multicolinearity)
  vif_model <- lm(as.formula(paste0(resp, " ~ DEV_z + uncon_density_z + conv_density_z")),
                  data = df_spring)
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
    title = "Effect sizes across models"
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
  model <- lmer(formula, data = df_spring, REML = FALSE)
  
  #check model residuals
  x <- check_model(model)
  
  #check vif (no multicolinearity)
  vif_model <- lm(as.formula(paste0(resp, " ~ DEV_z + uncon_density_z + conv_density_z")),
                  data = df_spring)
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
    title = "Effect sizes across models"
  ) +
  theme_bw(base_size = 20) +   
  theme(
    panel.grid = element_blank(),  
    strip.text = element_text(size = 11, face = "bold"),
    panel.spacing = unit(0.5, "lines"),   # restore spacing between facets
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", linewidth = 0.75, fill = NA) 
  )