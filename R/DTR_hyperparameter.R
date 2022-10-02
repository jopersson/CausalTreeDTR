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

source("R/data_generating_process.R")
source("R/DTR.R")
source("R/evaluation.R")

## Simulation Function
hyperparameter_simulation <- function(iter = 10, n = 5000,
                                      scenarios, scenarios_names, 
                                      approaches, approaches_names, 
                                      parameter = c(20, 100))
{
  
  
  # initialize data frames for training and testing regret
  evaluation <- data.frame(matrix(NA, nrow = iter*length(approaches)*length(scenarios)*length(parameter)*2, ncol = 6))
  colnames(evaluation) <- c("Accuracy", "Regret", "Approach", "Scenario", "Data", "Parameter")
  
  indices <- seq(1,length(approaches)*length(scenarios)*length(parameter)*iter)
  index1 <- indices[1]
  index2 <- length(approaches)*length(scenarios)*length(parameter)*iter + index1
  
  for(i in 1:length(scenarios)){
    for(j in 1:length(approaches)){
      for(k in 1:length(parameter)){
        for(it in 1:iter){
          # declarations
          evaluation[index1, "Scenario"] <- (scenarios_names[i])
          evaluation[index2, "Scenario"] <- (scenarios_names[i])
          
          evaluation[index1, "Approach"] <- (approaches_names[j])
          evaluation[index2, "Approach"] <- (approaches_names[j])
          
          evaluation[index1, "Parameter"] <- (parameter[k])
          evaluation[index2, "Parameter"] <- (parameter[k])
          
          evaluation[index1, "Data"] <- "train"
          evaluation[index2, "Data"] <- "test"
          
          # evaluate approaches
          evaluated <- evaluate.DTR(scenario = scenarios[i], 
                                    approach = approaches[j],
                                    minsize = parameter[k])
          
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
  }
  
  
  # factor declarations 
  evaluation$Approach <- as.factor(evaluation$Approach)
  evaluation$Scenario <- as.factor(evaluation$Scenario)
  evaluation$Data <- as.factor(evaluation$Data)
  evaluation$Parameter <- as.factor(evaluation$Parameter)
  
  
  out <- list(evaluation = evaluation)
  return(out)
}


# --------------------------------------------------------------------------------------------------------

# Initialization
## Data Scenarios
scenarios <- c("linear", 
               "linear_interaction", 
               "nonlinear", 
               "nonlinear_interaction")
scenarios_names <- c("linear", 
                     "linear interactions", 
                     "non-linear", 
                     "non-linear interactions")

## Approaches
approaches <- c(#"qlearn", 
                #"dwols", 
                #"gestimation", 
                "knn",
                "cart",
                "causaltree", 
                "causalforest")
approaches_names <- c(#"Q-Learning", 
                      #"DWOLS", 
                      #"G-Estimation",
                      "k-NN Regression",
                      "CART",
                      "DTR-CT", 
                      "DTR-CF")


## Number of iterations and pupulation size and hyperparameter
iter <- 100
n <- 5000
minsize <- c(10,20,30,40,50,60,70,80)

# Run simulation to get evaluations
simulation <- hyperparameter_simulation(iter = iter, n = n, 
                                        scenarios = scenarios, scenarios_names = scenarios_names, 
                                        approaches = approaches,  approaches_names = approaches_names,
                                        parameter = minsize)

evaluation <- simulation$evaluation

# Results and Outcomes
write.csv(evaluation, 
          "Results/evaluation_hyperparameter.csv", 
          row.names = F)

#evaluation <- read.csv("Results/evaluation_hyperparameter.csv")

plot_tune(evaluation, metric = "Regret")
plot_tune(evaluation, metric = "Accuracy")

plot_tune_big(evaluation, metric = "Regret", data = "train")
plot_tune_big(evaluation, metric = "Accuracy", data = "train")
plot_tune_big(evaluation, metric = "Regret", data = "test")
plot_tune_big(evaluation, metric = "Accuracy", data = "test")





