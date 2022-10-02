library(mgcv)
library(caret)
library(ggplot2)
library(ggpubr)
library(causalTree)
library(grf)
library(dplyr)
library(tidyverse)
library(remotes)
library(lubridate)
library(simputation)
library(mmtable2)
library(gt)
library(forcats)
library(tidyr)
library(purrr)
library(FNN)

source("R/data_generating_process.R")
source("R/DTR.R")


# ------------------------------ SIMULATION -------------------------------
table_eval <- function(evaluation, metric = "Regret", title_metric = "Regret"){
  filename <- paste0("table_", metric, ".pdf", collapse = "")
  
  # data frames
  data_wrangled <- evaluation %>%
    group_by(Approach, Data, Scenario) %>%
    summarise(across(.cols = c(metric), .fns = c(mean, sd))) %>%
    ungroup() 
  
  
  data_wrangled <- data_wrangled  %>% 
    rename_with(.cols=c(4,5), ~c("mean", "sd")) %>% 
    filter(Scenario == "linear interactions" | Scenario == "non-linear interactions") %>% 
    mutate(Value = paste(sprintf("%.3f",round(mean,3))," (", sprintf("%.3f",round(sd,3)), ") ", sep = "")) 
  
  # order levels appropriately
  data_wrangled <- data_wrangled %>% mutate(Data = fct_relevel(Data, "train", "test"),
                                            Approach = fct_relevel(Approach, "Q-Learning", "DWOLS", "G-Estimation", "CART", "k-NN Regression", "DTR-CT", "DTR-CF"))
  
  # tables
  table_mm <- data_wrangled %>% 
    mmtable(table_data = Value, table_name = "Values") +
    header_top(Data)  +   
    header_left(Approach) + 
    header_left_top(Scenario)
  
  obj <- table_mm %>% 
    gt::tab_header(title = paste(title_metric, "of DTR across Approaches and Scenarios for Training and Test Data"),
                   subtitle = "Mean (SD)") %>% 
    tab_options(table.width = 700, container.overflow.x = F, container.overflow.y = F) %>% 
    tab_style(style = list(cell_text(weight = "bold")), locations = cells_title("title")) 
  
  gtsave(obj, filename, path = "Results", zoom = .5)
  return(obj)
}



# ----------------------- HYPERPARAMETER SIMULATION -----------------------

plot_tune <- function(evaluation, 
                      metric = "Regret")
{
  filename <- paste0("tune_", metric, ".pdf", collapse = "", sep = "")
  
  df <- evaluation %>% select_at(vars(metric, Scenario, Approach, Data, Parameter)) %>% 
    filter(Data == "test", 
           Approach == "DTR-CT" | Approach == "DTR-CF",      # | Approach == "k-NN Regression"
           Scenario != "non-linear", Scenario != "linear") %>%      
    mutate(Data = fct_relevel(Data, "train", "test"),
           Approach = fct_relevel(Approach, approaches_names)) %>% 
    group_by_at(vars(Approach, Scenario, Data, Parameter)) %>% 
    summarise(mean=mean(.data[[metric]]), sd=sd(.data[[metric]])) 
  
  df2 <- df %>% select(Approach, Scenario, Data, Parameter, mean, sd) %>%
    mutate(lwr = mean-sd, upr = mean + sd) %>% select(-sd)
  
  df2$Parameter <- as.numeric(df2$Parameter)
  
  position <- "none"
  stript <-  element_text()
  height <- 2.5
  if(metric == "Accuracy"){
    position <- "bottom"
    stript <-  element_blank()
    height <- 2.8
  }
  
  plot <- ggplot(df2, aes(x = Parameter, y = mean, color = Approach)) +
    geom_line(aes(x = Parameter, y = mean, color = Approach, group = Approach), size = 1) + 
    geom_ribbon(aes(y = mean, ymin = lwr, ymax = upr, group = Approach, colour = Approach, fill = Approach), 
                alpha = 0.2, colour = NA) +
    facet_grid(. ~ Scenario) +
    theme_minimal() +
    theme(legend.position= position,
          legend.title=element_blank(),
          axis.title.y = element_blank(),
          strip.text = stript,
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "pt")) + 
    xlab ("") + 
    scale_color_manual(values=c(#"#A73561", 
      "#56B4E9", "#08407E")) + 
    scale_fill_manual(values=c(#"#A73561", 
      "#56B4E9","#08407E")) +
    scale_x_continuous(breaks = seq(10,100,10))
  
  print(plot)
  ggsave(filename, 
         plot = plot, 
         device = "pdf",
         width = 5, height = height,
         path = "Results")
}

plot_tune_big <- function(evaluation, 
                          metric = "Regret",
                          data = "test")
{
  filename <- paste0("tune_big_", metric, "_", data,".pdf", collapse = "", sep = "")
  
  df <- evaluation %>% select_at(vars(metric, Scenario, Approach, Data, Parameter)) %>% 
    filter(Scenario != "non-linear", Scenario != "linear",
           Data == data) %>%      
    mutate(Data = fct_relevel(Data, "train", "test"),
           Approach = fct_relevel(Approach, approaches_names)) %>% 
    group_by_at(vars(Approach, Scenario, Data, Parameter)) %>% 
    summarise(mean=mean(.data[[metric]]), sd=sd(.data[[metric]])) 
  
  df2 <- df %>% select(Approach, Scenario, Data, Parameter, mean, sd) %>%
    mutate(lwr = mean-sd, upr = mean + sd) %>% select(-sd)
  
  df2$Parameter <- as.numeric(df2$Parameter)*10
  
  position <- "none"
  stript <-  element_text()
  height <- 2.5
  if(data == "test"){
    position <- "bottom"
    stript <-  element_blank()
    height <- 2.8
  }
  
  plot <- ggplot(df2, aes(x = Parameter, y = mean, color = Approach)) +
    geom_line(aes(x = Parameter, y = mean, color = Approach, group = Approach), size = 1) + 
    geom_ribbon(aes(x = Parameter, y = mean, ymin = lwr, ymax = upr, group = Approach, colour = Approach, fill = Approach), 
                alpha = 0.2, colour = NA) +
    facet_grid(. ~ Scenario) +
    theme_minimal() +
    theme(legend.position= position,
          legend.title=element_blank(),
          axis.title.y = element_blank(),
          strip.text = stript,
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "pt")) + 
    xlab ("") + 
    scale_color_manual(values=c("#A73561", "#8E6713", "#56B4E9", "#08407E")) + 
    scale_fill_manual(values=c("#A73561", "#8E6713", "#56B4E9","#08407E")) +
    scale_x_continuous(breaks = seq(10,80,10))
  
  print(plot)
  ggsave(filename, 
         plot = plot, 
         device = "pdf",
         width = 5, height = height,
         path = "Results")
}




# -------------------------------- MIMIC ----------------------------------
table_reward <- function(df){
  filename <- paste0("table_reward.pdf", collapse = "")
  
  colnames(df) <- c("Approach", "Value", "Treatment")
  df$Treatment <- as.factor(df$Treatment)
  
  # order levels appropriately
  df <- df %>% mutate(Treatment = fct_relevel(Treatment, "Ventilation", "Vasopressin"),
                      Approach = fct_relevel(Approach, "observed", "Q-Learning", "DWOLS", "G-Estimation", "CART", "k-NN Regression", "DTR-CT", "DTR-CF"))
  
  # # tables
  table_mm <- df %>% mmtable(table_data = Value, table_name = "Cumulative Rewards") +
    header_top(Treatment)  +  header_left(Approach) 
  
  obj <- table_mm %>% 
    gt::tab_header(title = paste("Cumulative Rewards of DTRs across Approaches and Treatments")) %>% 
    tab_options(table.width = 500, container.overflow.x = F, container.overflow.y = F) %>% 
    tab_style(style = list(cell_text(weight = "bold")), locations = cells_title("title")) 
  
  gtsave(obj, filename, path = "Results", zoom = .3)
  return(obj)
}


get_tree_plot <- function(fitmimic, t = 10){
  t1 <- rpart.plot(fitmimic_tree1$trees[[t]], 
                   cex = 1, type = 0, 
                   tweak =.8, fallen.leaves =F, 
                   faclen = 4, varlen = 9, 
                   clip.facs = T, 
                   branch=.9,
                   box.palette = "RdBu")
  return(t1)
}

get_variable_importance <- function(fitmimic = fitmimic_forest2, t = 10, top = 10, plot = T, treatment = "Vasopressin"){
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  #fitmimic <- fitmimic_forest2
  Variables <- colnames(fitmimic$forests[[t]]$X.orig)
  Importance <- variable_importance(fitmimic$forests[[t]])
  
  df <- as.data.frame(Variables)
  df$Importance <- (Importance)
  df <- df %>% arrange(desc(Importance)) %>% top_n(top)
  
  if(plot == T){
    vip <- ggplot(df) + 
      geom_col(aes(y = reorder(Variables, Importance), x = Importance), fill = "#08407E", color = NA) + 
      theme(axis.title.y = element_blank(),
            axis.text=element_text(size=18),
            axis.title.x = element_text(size=18),
            title = element_text(size = 18),
            plot.margin = unit(c(0, 15, 0, 0), "pt")) + 
      ggtitle(paste0("Treatment: ", treatment)) +  
      scale_y_discrete(expand = c(0,0)) +
      scale_x_continuous(labels = scaleFUN)
    
    print(vip)
  }
  out <- list(importance = df, 
              plot = vip)
  return(out)
}
