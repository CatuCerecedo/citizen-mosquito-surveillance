############################## BG Model (GLMM version) #########################

# NOT RUN the Script
# This scripts represent the followed process to select variables
# It should not include in the whole modelling process once we decided the explanatory variables

library(tidyverse)
library(lme4)
library(MuMIn)
library(ggplot2)
library(glmmTMB)
library(caret)
library(parallel)

rm(list = ls())

loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")

tiger <- readRDS(paste0(loc.output, "bg_tiger_spain.rds"))
  
# Functions to check colinearity -----------------------------------------------
# Checking collinearity --------------------------------------------------------

#Function to obtain the variable with the highest VIF

get_variable_with_highest_vif <- function(tmb_full) {
  
  # Obtaining VIF Values
  vif <- performance::check_collinearity(tmb_full)
  
  variable <- vif$Term[which.max(vif$VIF)]
  return(variable)
}

# Function to delete variables with large VIF
delate_variables <- function(tmb_full, m.vif = 4){
  
  # Loop to delete variables with the highest vif
  while (TRUE) {
    vif_model <- performance::check_collinearity(tmb_full)
    max_vif <- max(vif_model$VIF)
    
    if (max_vif <= m.vif) {
      print("stop")
      break
    }
    
    variable_to_remove <- get_variable_with_highest_vif(tmb_full)
    
    print(paste0(variable_to_remove, " VIF:", max_vif)) # VIF > 4, change it
    
    formula_updated <- as.formula(paste(". ~ . -", variable_to_remove))
    
    # Ajustar el modelo con la nueva fórmula
    tmb_full <- update(tmb_full, formula = formula_updated)
  }
  return(tmb_full)
}

# Correlations -----------------------------------------------------------------
correlations <- data.frame(
  variable = names(tiger[9:30]),
  cor = as.vector(cor(tiger$females, tiger[9:30], method = "spearman"))
) 
correlations <- correlations %>% mutate(
  positive = ifelse(cor > 0, "TRUE", "FALSE")
)

# Full model -------------------------------------------------------------------
func_f <- females ~ poly(mean_temperature, 2) +  poly(precipitation, 2) + 
  mean_relative_humidity + wind_speed +
  offset(log(trapping_effort)) + offset(log(n_traps)) +
  (1 | prov_name)

tmb_full <- glmmTMB(func_f,
                    data = tiger, family = nbinom2)

tmb_full <- delate_variables(tmb_full)

options(na.action = "na.fail")
tmb_full_dredge <- dredge(tmb_full, cluster = 4, rank = "AICc", trace = 2, # fixed = "cond(log(trapping_effort))", # getAllTerms(tmb_full)
                          m.lim = c(2, 7))
options(na.action = "na.omit")

func_f <- females ~ poly(mean_temperature, 2) +  poly(precipitation, 2) + 
  mean_relative_humidity + wind_speed +
  offset(log(trapping_effort)) + offset(log(n_traps)) +
  (1 | prov_name)

tmb_full <- glmmTMB(func_f,
                    data = tiger, family = nbinom2)

tmb_full <- delate_variables(tmb_full)

options(na.action = "na.fail")
tmb_full_dredge <- dredge(tmb_full, cluster = 4, rank = "AICc", trace = 2, # fixed = "cond(log(trapping_effort))", # getAllTerms(tmb_full)
                          m.lim = c(2, 7))
options(na.action = "na.omit")

func_f <- females ~ poly(mean_temperature, 2) +  
 artificial_green_urban + forests + inland_water + inland_wet +
  offset(log(trapping_effort)) + offset(log(n_traps)) +
  (1 | prov_name)

tmb_full <- glmmTMB(func_f,
                    data = tiger, family = nbinom2)

options(na.action = "na.fail")
tmb_full_dredge <- dredge(tmb_full, cluster = 4, rank = "AICc", trace = 2, # fixed = "cond(log(trapping_effort))", # getAllTerms(tmb_full)
                          m.lim = c(2, 7))
options(na.action = "na.omit")
                   