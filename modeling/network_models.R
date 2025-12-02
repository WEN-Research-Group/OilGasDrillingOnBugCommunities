#network modeling
rm(list = ls())

setwd("C:/Users/rpolivie/OneDrive - Syracuse University/Tao DATA") 

getwd()
list.files()

install.packages(c('lme4','influence.ME',"performance","DHARMa","car","sjPlot"))
library(lme4)
library(dunn.test)
library(performance) 
library(DHARMa)   
library(car)         
library(influence.ME)
library(sjPlot)
library(lmerTest)
library(vegan)
library(ggplot2)
library(tidyverse)
library(stringr)
library(dunn.test)   
library(ggpubr)     

PANetworkAnalysis <- read.table("BMI_Network_CooccurrenceAnalysis.txt", 
                                header = TRUE, 
                                sep = "\t", 
                                fill = TRUE, 
                                quote = "",
                                stringsAsFactors = FALSE)


#########################################################################

#mixed modeling of netowrk metrics and composition using huc8 as random effect

#scale variables so results from different models can be compared to each other
PANetworkAnalysis <- PANetworkAnalysis %>% mutate(
  M_pos_z = scale(M_pos)[,1],
  N_pos_z = scale(N_pos)[,1],
  C_pos_z = scale(C_pos)[,1],
  S_pos_z = scale(S_pos)[,1],
  MeanPTV_z = scale(MeanPTV)[,1]
)


response_vars <- c("M_pos_z", "N_pos_z", "C_pos_z", "S_pos_z", "MeanPTV_z")
model_results <- list()

for (resp in response_vars) {
  
  #some formatting so i can tell which results go to what
  cat("\n=========================\n")
  cat("Response:", resp, "\n")
  cat("=========================\n")
  
  #formula- u_bool and C_bool as fixed effects with HUC8 as random effect
  formula <- as.formula(paste0(
    resp, " ~ U_bool + C_bool + (1 | HUC8)"
  ))
  
  #REML false for model comparison
  model <- lmer(formula, data = PANetworkAnalysis, REML = FALSE)
  
  #check model residuals
  x <- check_model(model)
  
  #check vif (no multicolinearity)
  vif_model <- lm(as.formula(paste0(resp, " ~ U_bool + C_bool")),
                  data = PANetworkAnalysis)
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
coef_df <- coef_df %>% filter(term %in% c("U_bool", "C_bool"))

coef_df <- coef_df %>%
  mutate(
    model = factor(model, levels = c("N_pos_z","C_pos_z","M_pos_z","S_pos_z","MeanPTV_z")),
    term = factor(term, levels = c("C_bool", "U_bool")),
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

######################################################################  #funtional composition of networks


ffg_data <- read.csv('PADEP Aquatic Macroinvertebrate Taxa.csv')
unique_taxa_lists <- PANetworkAnalysis %>%
  select(Taxa_List, U_bool, C_bool, HUC8) %>%
  mutate(List_ID = row_number())   # keep a unique ID for each row

# Step 2: Split taxa into long format
taxa_long <- unique_taxa_lists %>%
  mutate(Taxa = str_split(Taxa_List, ",\\s*")) %>%
  unnest(Taxa) %>%
  rename(Taxon = Taxa)

# Step 3: Join with FFG lookup
# Keep only the most frequent FFG per taxon
ffg_data_unique <- ffg_data %>%
  group_by(ASSESSMENT_ID, DEP_LAB_ID_FFG) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ASSESSMENT_ID) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%  # keep the most frequent FFG
  ungroup() %>%
  select(ASSESSMENT_ID, DEP_LAB_ID_FFG)

taxa_ffg <- taxa_long %>%
  left_join(ffg_data_unique, by = c("Taxon" = "ASSESSMENT_ID")) %>%
  filter(!is.na(DEP_LAB_ID_FFG))

# Step 4: Count FFGs and compute proportions per original row
ffg_proportions <- taxa_ffg %>%
  group_by(List_ID, DEP_LAB_ID_FFG) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ffg_proportions <- ffg_proportions %>%
  filter(DEP_LAB_ID_FFG != "")

# Now pivot wider
ffg_wide <- ffg_proportions %>%
  select(-n) %>%
  pivot_wider(names_from = DEP_LAB_ID_FFG,
              values_from = prop,
              values_fill = 0)

# Step 6: Add back metadata (HUC8, Taxa_List, U/C flags)
ffg_wide <- ffg_wide %>%
  left_join(unique_taxa_lists, by = "List_ID") %>%
  mutate(U_C_group = paste0("U", U_bool, "_C", C_bool))

ffg_wide <- ffg_wide %>%
  mutate(across(c(CG, FC, PR, SC, SH), ~ . * 100))

#####
response_vars <- c("CG", "FC", "PR", "SC", "SH")
model_results <- list()

for (resp in response_vars) {
  
  #some formatting so i can tell which results go to what
  cat("\n=========================\n")
  cat("Response:", resp, "\n")
  cat("=========================\n")
  
  #formula- u_bool and C_bool as fixed effects with HUC8 as random effect
  formula <- as.formula(paste0(
    resp, " ~ U_bool + C_bool + (1 | HUC8)"
  ))
  
  #REML false for model comparison
  model <- lmer(formula, data = ffg_wide, REML = FALSE)
  
  #check model residuals
  x <- check_model(model)
  
  #check vif (no multicolinearity)
  vif_model <- lm(as.formula(paste0(resp, " ~ U_bool + C_bool")),
                  data = ffg_wide)
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


#use saved model results
coef_df <- bind_rows(
  lapply(names(model_results), function(name) {
    broom.mixed::tidy(model_results[[name]]$model, effects = "fixed") %>%
      mutate(model = name)
  })
)

#just use fixed effects
coef_df <- coef_df %>% filter(term %in% c("U_bool", "C_bool"))

#add significance stars
coef_df <- coef_df %>%
  mutate(
    model = factor(model, levels = c("CG", "FC", "PR", "SC", "SH")),
    term = factor(term, levels = c("C_bool", "U_bool")),
    p.signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    
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



