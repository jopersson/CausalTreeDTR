library(dplyr)
library(simputation)


prepare <- function(K = 10, lag = 3, treatment_name = "vent"){
  
  # select relevant patient information
  treatment_name <- treatment_name
  baseline_names <- c("age", "ethnicity", "gender", "insurance", "admission_type", "first_careunit")
  moderator_names <- c("diastolic blood pressure", "heart rate", "mean blood pressure", "oxygen saturation", "respiratory rate", "systolic blood pressure", "temperature")
  
  Y <- outcome %>% select(subject_id, los_icu) %>% mutate(Y = -los_icu) %>% select(-los_icu)
  X0 <- baseline %>% select(subject_id, all_of(baseline_names))
  df_patient <- merge(Y,X0, by = "subject_id")
  
  # get vita signs and treatment from first K time steps at the ICU
  Xt_selection <- vitals_labs_mean %>% select(subject_id, all_of(moderator_names))
  colnames(Xt_selection) <- c("subject_id", "blood_pressure_diastolic", "heart_rate", "blood_pressure_mean", "oxygen_saturation", "respiratory_rate", "blood_pressure_systolic", "temperature")
  Xt_selection <- Xt_selection %>% filter(subject_id %in% df_patient$subject_id)
  Xt_selection_reduced <- Xt_selection %>% group_by(subject_id) %>% filter(row_number() == c(1:(K)))
  Xt_selection_reduced <- Xt_selection_reduced %>% group_by(subject_id) %>% slice(1:n())
  
  treatments <- interventions %>% select(subject_id, all_of(treatment_name))
  treatments_reduced <- treatments %>% filter(subject_id %in% df_patient$subject_id) %>% group_by(subject_id) %>% filter(row_number() == 1:(K))
  treatments_reduced <- treatments_reduced %>% group_by(subject_id) %>% slice(1:n())
  
  datagen <- df_patient
  moderators <- list()
  for(t in 1:K){
    
    Atemp <- treatments_reduced %>% group_by(subject_id) %>% filter(row_number() == t)
    colnames(Atemp)[2] <- paste0("A",t)
    
    Xtemp <- Xt_selection_reduced %>% group_by(subject_id) %>% filter(row_number() == t)
    colnames(Xtemp)[2:length(colnames(Xtemp))] <- paste0(colnames(Xtemp)[2:length(colnames(Xtemp))],t)
    moderators[[t]] <- colnames(Xtemp)[2:length(colnames(Xtemp))]
    
    temp <- merge(Atemp, Xtemp, by = "subject_id")
    datagen <- merge(datagen, temp, by = "subject_id")
  }
  
  
  # data set formatted for DTRensemble function
  MIMIC <- datagen
  missing <- as.data.frame(apply(MIMIC, 2,function(x) sum(is.na(x))))
  
  # impute missing values in vital signs
  colnames <- colnames(MIMIC)
  which <- colnames(MIMIC)[8:dim(MIMIC)[2]]
  remove <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10")
  temperatures <- paste0("temperature", c(1:10))
  
  MIMICimputed <- MIMIC
  subset <- colnames[which(!colnames %in% remove)]
  subset <- subset[9:length(subset)]
  toimpute <- subset[which(!subset %in% temperatures)]
  
  MIMIC <- sapply(MIMIC[, subset], as.numeric)
  MIMIC <- as.data.frame(MIMIC)
  
  impute_median <- function(x){
    ifelse(is.na(x), median(x, na.rm = TRUE), x)
  }
  
  MIMICimputed <- MIMICimputed %>% mutate_at(vars(toimpute), impute_median)
  MIMICimputed <- MIMICimputed %>% impute_lm(temperature1 ~ blood_pressure_diastolic1 + heart_rate1 + blood_pressure_mean1 + oxygen_saturation1 + respiratory_rate1 + blood_pressure_systolic1)
  MIMICimputed <- MIMICimputed %>% impute_lm(temperature2 ~ blood_pressure_diastolic2 + heart_rate2 + blood_pressure_mean2 + oxygen_saturation2 + respiratory_rate2 + blood_pressure_systolic2)
  MIMICimputed <- MIMICimputed %>% impute_lm(temperature3 ~ blood_pressure_diastolic3 + heart_rate3 + blood_pressure_mean3 + oxygen_saturation3 + respiratory_rate3 + blood_pressure_systolic3)
  MIMICimputed <- MIMICimputed %>% impute_lm(temperature4 ~ blood_pressure_diastolic4 + heart_rate4 + blood_pressure_mean4 + oxygen_saturation4 + respiratory_rate4 + blood_pressure_systolic4)
  MIMICimputed <- MIMICimputed %>% impute_lm(temperature5 ~ blood_pressure_diastolic5 + heart_rate5 + blood_pressure_mean5 + oxygen_saturation5 + respiratory_rate5 + blood_pressure_systolic5)
  MIMICimputed <- MIMICimputed %>% impute_lm(temperature6 ~ blood_pressure_diastolic6 + heart_rate6 + blood_pressure_mean6 + oxygen_saturation6 + respiratory_rate6 + blood_pressure_systolic6)
  MIMICimputed <- MIMICimputed %>% impute_lm(temperature7 ~ blood_pressure_diastolic7 + heart_rate7 + blood_pressure_mean7 + oxygen_saturation7 + respiratory_rate7 + blood_pressure_systolic7)
  MIMICimputed <- MIMICimputed %>% impute_lm(temperature8 ~ blood_pressure_diastolic8 + heart_rate8 + blood_pressure_mean8 + oxygen_saturation8 + respiratory_rate8 + blood_pressure_systolic8)
  MIMICimputed <- MIMICimputed %>% impute_lm(temperature9 ~ blood_pressure_diastolic9 + heart_rate9 + blood_pressure_mean9 + oxygen_saturation9 + respiratory_rate9 + blood_pressure_systolic9)
  MIMICimputed <- MIMICimputed %>% impute_lm(temperature10 ~ blood_pressure_diastolic10 + heart_rate10 + blood_pressure_mean10 + oxygen_saturation10 + respiratory_rate10 + blood_pressure_systolic10)
  
  formulas <- matrix(NA, nrow=K, ncol=K)
  for(t in 1:K){
    start <- t 
    end <- min(K,t + lag - 1)
    
    print(c(start, end))
    formulas[start:end,t] <- paste(moderators[[t]], collapse = "+", sep = " ")
    
    formulae <- apply(formulas, 1, function(x) paste(x[!is.na(x)], collapse = " + ", sep = " "))
  }
  
  Aresponse <- paste0("A",1:K)
  formulae0 <- paste(baseline_names, collapse = "+")
  formulae0 <- paste("~", formulae0)
  formulae <- apply(expand.grid(formulae0, formulae), 1, paste, collapse="+")
  formulaeY <- paste("Y", formulae)
  formulaeA <- apply(cbind(Aresponse, formulae), 1, paste, collapse="")
  
  formulaeA <- sapply(formulaeA, as.formula)
  formulae <- sapply(formulae,  as.formula)
  formulaeY <- sapply(formulaeY,  as.formula)
  
  out <- list(Y = Y,
              MIMIC = MIMICimputed,
              formulae = formulae, 
              formulaeA = formulaeA)
  
  return(out)
}
