# Libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(causalTree)
library(grf)
library(mmtable2)
library(gt)
library(forcats)
library(tidyr)
library(purrr)
library(FNN)

source("R/data_generating_process.R")
source("R/DTR.R")
source("R/evaluation.R")

set.seed(2021)

## Simulation Function
simulate <- function(iter = 100, n = 5000,
                     scenarios, scenarios_names, 
                     approaches, approaches_names)
{
  
  
  # initialize data frames for training and testing regret
  evaluation <- data.frame(matrix(NA, nrow = iter*length(approaches)*length(scenarios)*2, ncol = 5))
  colnames(evaluation) <- c("Accuracy", "Regret", "Approach", "Scenario", "Data")
  
  indices <- seq(1,length(approaches)*length(scenarios)*iter)
  index1 <- indices[1]
  index2 <- length(approaches)*length(scenarios)*iter + index1
  
  for(i in 1:length(scenarios)){
    for(j in 1:length(approaches)){
      for(it in 1:iter){
        # declarations
        evaluation[index1, "Scenario"] <- (scenarios_names[i])
        evaluation[index2, "Scenario"] <- (scenarios_names[i])
        
        evaluation[index1, "Approach"] <- (approaches_names[j])
        evaluation[index2, "Approach"] <- (approaches_names[j])
        
        evaluation[index1, "Data"] <- "train"
        evaluation[index2, "Data"] <- "test"
        
        # evaluate approaches
        evaluated <- evaluate.DTR(scenario = scenarios[i], approach = approaches[j])
        
        # fill in evaluations of training and test data
        evaluation[index1, "Accuracy"] <- evaluated$train$accuracy
        evaluation[index2, "Accuracy"] <- evaluated$test$accuracy
        
        evaluation[index1, "Regret"] <- evaluated$train$regret
        evaluation[index2, "Regret"] <- evaluated$test$regret
        
        # update
        index1 <- index1 + 1
        index2 <- index2 + 1
      }
    }
  }
  
  
  # factor declarations 
  evaluation$Approach <- as.factor(evaluation$Approach)
  evaluation$Scenario <- as.factor(evaluation$Scenario)
  evaluation$Data <- as.factor(evaluation$Data)
  
  
  out <- list(evaluation = evaluation)
  return(out)
}

# --------------------------------------------------------------------------------------------------------
# Initialization
## Population Size
n <- 5000

## Data Scenarios
scenarios <- c("linear_interaction", 
               "nonlinear_interaction")
scenarios_names <- c("linear interactions", 
                     "non-linear interactions")

## approaches
approaches <- c("qlearn", 
                "dwols", 
                "gestimation",
                "cart",
                "knn",
                "causaltree", 
                "causalforest")
approaches_names <- c("Q-Learning", 
                      "DWOLS", 
                      "G-Estimation",
                      "CART",
                      "k-NN Regression",
                      "DTR-CT", 
                      "DTR-CF")

## Number of iterations
iter <- 100

# Run simulation to get evaluations
simulation <- simulate(iter = iter, 
                       n = n, 
                       scenarios = scenarios, 
                       scenarios_names = scenarios_names, 
                       approaches = approaches, 
                       approaches_names = approaches_names)

evaluation <- simulation$evaluation


# Results and Outcomes
write.csv(evaluation, 
          "Results/evaluation_simulation.csv", 
          row.names = F)

#evaluation <- read.csv("Results/evaluation_simulation.csv")

table_eval(evaluation, 
           metric = "Regret", 
           title_metric = "Cumulative Regret")
table_eval(evaluation, 
           metric = "Accuracy", 
           title_metric = "Accuracies")