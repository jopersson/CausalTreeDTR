# Libraries
library(dplyr)
library(causalTree)
library(grf)
library(tree)
library(FNN)
library(rpart)
library(tidyverse)

source("R/data_generating_process.R")
source("R/DTR.R")

# Initialization
n <- 5000

# Inspect Data
data1 <- data_generation(n, model = "linear_interaction")
mean(Y.opt1 <- data1$data$Y.opt)
mean(Y.obs1 <- data1$data$Y)

data2 <- data_generation(n, model = "nonlinear_interaction")
mean(Y.opt2 <- data2$data$Y.opt)
mean(Y.obs2 <- data2$data$Y)

# Run DTR
mean(data1$data$Y.opt)
mean(data1$data$Y)
evaluate.DTR(scenario = "linear_interaction", approach = "qlearn", n)
evaluate.DTR(scenario = "linear_interaction", approach = "dwols", n)
evaluate.DTR(scenario = "linear_interaction", approach = "gestimation", n)
evaluate.DTR(scenario = "linear_interaction", approach = "knn", n)
evaluate.DTR(scenario = "linear_interaction", approach = "cart", n)
evaluate.DTR(scenario = "linear_interaction", approach = "causaltree", n)
evaluate.DTR(scenario = "linear_interaction", approach = "causalforest", n)


mean(data2$data$Y.opt)
mean(data2$data$Y)
evaluate.DTR(scenario = "nonlinear_interaction", approach = "qlearn", n)
evaluate.DTR(scenario = "nonlinear_interaction", approach = "dwols", n)
evaluate.DTR(scenario = "nonlinear_interaction", approach = "gestimation", n)
evaluate.DTR(scenario = "nonlinear_interaction", approach = "knn", n)
evaluate.DTR(scenario = "nonlinear_interaction", approach = "cart", n)
evaluate.DTR(scenario = "nonlinear_interaction", approach = "causaltree", n)
evaluate.DTR(scenario = "nonlinear_interaction", approach = "causalforest", n)