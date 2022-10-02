library(causalTree)
library(grf)
library(mgcv)
library(caret)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(mmtable2)
library(remotes)
library(simputation)


# Import Other R Functions
source("R/DTR.R")
source("R/evaluation.R")
source("Data/mimic_preparation.R")

# ---------------------------------------------------------------------------------------
# Import Data
load("Data/mimic_imputed.Rdata")

# colnames(vitals_labs_mean)   # dynamic data
# colnames(baseline)           # static data
# colnames(interventions)      # treatments
# colnames(outcome)            # outcomes


# Preparation
K <- 10
lag <- 3
preparation_vent <- prepare(K, lag, "vent")
preparation_vaso <- prepare(K, lag, "vaso")


# Fit DTRs under different approaches for 2 treatments
fitmimic_qlearn1 <- DTR(Y, 
                        outcome.mod = preparation_vent$formulae, 
                        treatment.mod = preparation_vent$formulaeA,
                        data = preparation_vent$MIMIC, 
                        approach = "qlearn", 
                        simulation = F)
fitmimic_dwols1 <- DTR(Y, 
                       outcome.mod = preparation_vent$formulae, 
                       treatment.mod = preparation_vent$formulaeA,
                       data = preparation_vent$MIMIC, 
                       approach = "dwols", 
                       simulation = F)
fitmimic_gest1 <- DTR(Y, 
                      outcome.mod = preparation_vent$formulae, 
                      treatment.mod = preparation_vent$formulaeA,
                      data = preparation_vent$MIMIC, 
                      approach = "gestimation", 
                      simulation = F)
fitmimic_knn1 <- DTR(Y, 
                     outcome.mod = preparation_vent$formulae, 
                     treatment.mod = preparation_vent$formulaeA,
                     data = preparation_vent$MIMIC, 
                     approach = "knn", 
                     simulation = F)
fitmimic_cart1 <- DTR(Y, 
                      outcome.mod = preparation_vent$formulae, 
                      treatment.mod = preparation_vent$formulaeA,
                      data = preparation_vent$MIMIC, 
                      approach = "cart", 
                      simulation = F)
fitmimic_tree1 <- DTR(Y, 
                      outcome.mod = preparation_vent$formulae, 
                      treatment.mod = preparation_vent$formulaeA,
                      data = preparation_vent$MIMIC, 
                      approach = "causaltree", 
                      simulation = F)
fitmimic_forest1 <- DTR(Y, 
                        outcome.mod = preparation_vent$formulae, 
                        treatment.mod = preparation_vent$formulaeA,
                        data = preparation_vent$MIMIC, 
                        approach = "causalforest", 
                        simulation = F)
fitmimic_qlearn2 <- DTR(Y, 
                        outcome.mod = preparation_vaso$formulae, 
                        treatment.mod = preparation_vaso$formulaeA,
                        data = preparation_vaso$MIMIC, 
                        approach = "qlearn", 
                        simulation = F)
fitmimic_dwols2 <- DTR(Y, 
                       outcome.mod = preparation_vaso$formulae, 
                       treatment.mod = preparation_vaso$formulaeA,
                       data = preparation_vaso$MIMIC, 
                       approach = "dwols", 
                       simulation = F)
fitmimic_gest2 <- DTR(Y, 
                      outcome.mod = preparation_vaso$formulae, 
                      treatment.mod = preparation_vaso$formulaeA,
                      data = preparation_vaso$MIMIC, 
                      approach = "gestimation", 
                      simulation = F)
fitmimic_cart2 <- DTR(Y, 
                      outcome.mod = preparation_vaso$formulae, 
                      treatment.mod = preparation_vaso$formulaeA,
                      data = preparation_vaso$MIMIC, 
                      approach = "cart", 
                      simulation = F)
fitmimic_knn2 <- DTR(Y, 
                     outcome.mod = preparation_vaso$formulae, 
                     treatment.mod = preparation_vaso$formulaeA,
                     data = preparation_vaso$MIMIC, 
                     approach = "knn", 
                     simulation = F)
fitmimic_tree2 <- DTR(Y, 
                      outcome.mod = preparation_vaso$formulae, 
                      treatment.mod = preparation_vaso$formulaeA,
                      data = preparation_vaso$MIMIC, 
                      approach = "causaltree", 
                      simulation = F)
fitmimic_forest2 <- DTR(Y, 
                        outcome.mod = preparation_vaso$formulae, 
                        treatment.mod = preparation_vaso$formulaeA,
                        data = preparation_vaso$MIMIC, 
                        approach = "causalforest", 
                        simulation = F)

Y.opt.hat1 <- c(mean(preparation_vent$MIMIC$Y),
                mean(fitmimic_qlearn1$Y.opt),
                mean(fitmimic_dwols1$Y.opt),
                mean(fitmimic_gest1$Y.opt),
                mean(fitmimic_cart1$Y.opt),
                mean(fitmimic_knn1$Y.opt),
                mean(fitmimic_tree1$Y.opt),
                mean(fitmimic_forest1$Y.opt))
Y.opt.hat2 <- c(mean(preparation_vaso$MIMIC$Y),
                mean(fitmimic_qlearn2$Y.opt),
                mean(fitmimic_dwols2$Y.opt),
                mean(fitmimic_gest2$Y.opt),
                mean(fitmimic_cart2$Y.opt),
                mean(fitmimic_knn2$Y.opt),
                mean(fitmimic_tree2$Y.opt),
                mean(fitmimic_forest2$Y.opt))

approaches_names <- c("Q-Learning", 
                      "DWOLS", 
                      "G-Estimation",
                      "CART",
                      "k-NN Regression",
                      "DTR-CT", 
                      "DTR-CF")

cumulative_rewards <- as.data.frame(cbind(c("observed", approaches_names, "observed", approaches_names),
                                          c(round(Y.opt.hat1-mean(preparation_vent$MIMIC$Y),3), round(Y.opt.hat2-mean(preparation_vaso$MIMIC$Y),3)),
                                          c(rep("Ventilation",8), rep("Vasopressin",8))))

# Evaluation
table_reward(cumulative_rewards)

# Explainability Plots
## Causal Tree Illustration
venttree <- get_tree_plot(fitmimic_tree1, 
                          t = 10)
vasotree <- get_tree_plot(fitmimic_tree2, 
                          t = 10)

## TOP N important variables by Stage
ventvip <- get_variable_importance(fitmimic_forest1, 
                                   t = 10, 
                                   top = 10,
                                   treatment = "Ventilation")$plot

vasovip <- get_variable_importance(fitmimic_forest2, 
                                   t = 10, 
                                   top = 10,
                                   treatment = "Vasopressin")$plot 
