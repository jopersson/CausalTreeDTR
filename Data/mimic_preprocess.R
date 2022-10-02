# Libraries
library(mgcv)
library(caret)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(causalTree)
library(grf)
library(tidyverse)
library(remotes)
library(lubridate)
library(mmtable2)
library(gt)
library(boot)
library(htmlTable)
library(table1)
library(knitr)
library(gtsummary)
library(huxtable)
library(formattable)
library(htmltools)
library(webshot)
library(simputation)

# Initialization
codes <- read.csv("codes.csv", header = T)
interventions <- read.csv("interventions.csv", header = T)
intervention_name <- read.csv("intervention_name.csv", header = F)

outcome <- read.csv("outcome.csv", header = T)
baseline <- read.csv("baseline.csv", header = T)

vitals_labs_mean <- read.csv("vitals_labs_mean.csv", header = T)
vitals_labs_mean_name <- read.csv("dynamic_names.csv", header = T)


#------------------------------------ TREATMENTS -----------------------------------------
colnames(interventions)[1] <- "subject_id"
colnames(interventions)[2:15] <- gsub('^.|.$', '', intervention_name$V1)
colnames(interventions)[2:15] <- gsub('^.', '', colnames(interventions)[2:15])

interventions$subject_id <- as.factor(interventions$subject_id)
patient_interventions <- interventions %>% group_by(subject_id) %>% summarise(K = n(), pi_vent = sum(vent)/n(), vaso = sum(vaso), vent = sum(vent), adenosine = sum(adenosine), dobutamine = sum(dobutamine), epinephrine = sum(epinephrine), isuprel = sum(isuprel), milrinone = sum(milrinone), norepinephrine = sum(norepinephrine), phenylephrine = sum(phenylephrine), vasopressin = sum(vasopressin), colloid_bolus = sum(colloid_bolus), crystalloid_bolus = sum(crystalloid_bolus), nivdurations = sum(nivdurations))
intervention_mean <- patient_interventions %>% summarise(pi_vaso = mean(vaso/K), pi_vent = mean(vent/K), vaso = mean(vaso), vent = mean(vent), adenosine = mean(adenosine), dobutamine = mean(dobutamine), epinephrine = mean(epinephrine), isuprel = mean(isuprel), milrinone = mean(milrinone), norepinephrine = mean(norepinephrine), phenylephrine = mean(phenylephrine), vasopressin = mean(vasopressin), colloid_bolus = mean(colloid_bolus), crystalloid_bolus = mean(crystalloid_bolus), nivdurations = mean(nivdurations))

# Filtering
patient_interventions_reduced <- patient_interventions %>% select(subject_id, vent, pi_vent, K) %>% filter(K >= 20) %>% filter(pi_vent < 0.8) %>% filter(pi_vent > 0.2) #%>% filter(K < 50)
ind <- as.numeric(patient_interventions_reduced$subject_id)-1
interventions_reduced <- interventions[interventions$subject_id %in% ind,]
interventions_reduced %>% group_by(subject_id) %>% summarize(n = n())


# focus on ventilation first, then on vasopression as treatment
# get time series of treatments per row
get_trt_by_t <- function(data, K = 20){
  t <- seq(1,K)
  patients <- unique(data$subject_id)
  treatments <- data.frame(subject_id = patients)
  
  for(i in t){
    colname <- paste("A", i, sep = "")
    treatment <- data %>% group_by(subject_id) %>% filter(row_number()==i) %>% select(subject_id, vent)
    treatments[colname] <- treatment$vent
  }
  
  return(treatments)
}
treatments <- get_trt_by_t(data = interventions_reduced)



# ------------------------------- OUTCOME ----------------------------------------------
summary(outcome)
outcome$mort_icu <- as.factor(outcome$mort_icu)
outcome$mort_hosp <- as.factor(outcome$mort_hosp)

qplot(outcome$los_icu) 
qplot(outcome$mort_icu)
qplot(outcome$mort_hosp)

# ------------------------------- BASELINE COVARIATES ----------------------------------
summary(baseline)

baseline$age[baseline$age > 89] <- 90 #instead of 300
qplot(baseline$age)
qplot(baseline$gender)
levels(baseline$ethnicity) <- c("AMERICAN INDIAN", "AMERICAN INDIAN", "ASIAN", "ASIAN", "ASIAN", "ASIAN", "ASIAN", "ASIAN", "ASIAN", "ASIAN", "ASIAN", "ASIAN", 
                                "BLACK", "BLACK","BLACK","BLACK", "CARIBBEAN ISLAND", "HISPANIC", "HISPANIC", "HISPANIC", "HISPANIC", "HISPANIC", "HISPANIC", "HISPANIC", "HISPANIC", "HISPANIC", "HISPANIC", 
                                "MIDDLE EASTERN", "MULTI RACE ETHNICITY", "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", "UNKNOWN", "UNKNOWN", 
                                "PORTUGUESE", "SOUTH AMERICAN", "UNKNOWN", "UNKNOWN", "WHITE", "WHITE", "WHITE", "WHITE", "WHITE")
qplot(baseline$ethnicity)
qplot(baseline$insurance)
qplot(baseline$admission_type)
qplot(baseline$first_careunit)

baseline$admittime <- ymd_hms(baseline$admittime)
baseline$dischtime <- ymd_hms(baseline$dischtime)
baseline$intime <- ymd_hms(baseline$intime)
baseline$outtime <- ymd_hms(baseline$outtime)

time_in_hospital <- difftime(baseline$dischtime, baseline$admittime, units = "hours")
time_in_icu <- difftime(baseline$outtime, baseline$intime, units = "hours") 
time_to_icu <- difftime(baseline$intime, baseline$admittime, units = "mins")

qplot(time_in_hospital)
qplot(time_in_icu)
qplot(time_to_icu)

# Releveling
levels(baseline$ethnicity) <- c("Other", "Asian", "Black", "Other", "Hispanic", "Other", "Other", "Other", "Other", "Other", "Other", "White")
levels(baseline$admission_type) <- c("Elective", "Emergency", "Urgent")
levels(baseline$gender) <- c("Female", "Male")
# ---------------------------------- MODERATORS ------------------------------------------
colnames(vitals_labs_mean)[1] <- "subject_id"
colnames(vitals_labs_mean)[2:105] <- as.character(vitals_labs_mean_name$LEVEL2)

vitals_labs_mean$subject_id <- as.factor(vitals_labs_mean$subject_id)

# TODO DESCRIPTIVE STATISTICS OF MODERATORS
var <- colnames(vitals_labs_mean)[2:105]
vitals_labs_mean_reduced <- vitals_labs_mean %>% group_by(subject_id) %>% summarise_at(vars(var), mean, na.rm = TRUE)

nasum <- sapply(vitals_labs_mean, function(x) sum(is.na(x)))
naport <- round(nasum / nrow(vitals_labs_mean)*100,2)
naport <- as.data.frame(t(rbind(colnames(vitals_labs_mean), as.numeric(naport))))

selection <- c("diastolic blood pressure", "heart rate", "mean blood pressure", "oxygen saturation", "respiratory rate", "systolic blood pressure", "temperature")

vitals_labs_mean_selection <- vitals_labs_mean_reduced %>% select(subject_id, selection)


save(interventions,
     patient_interventions,
     interventions_reduced,
     patient_interventions_reduced,
     outcome, baseline,
     vitals_labs_mean, vitals_labs_mean_reduced,
     vitals_labs_mean_selection, file = "mimic_extended.RData")


# ----------------------------------- equal number of observations per time step ---------------
treatment_name <- "vent"
baseline_names <- c("age", "ethnicity", "gender", "insurance", "admission_type", "first_careunit")
moderator_names <- c("diastolic blood pressure", "heart rate", "mean blood pressure", "oxygen saturation", "respiratory rate", "systolic blood pressure", "temperature")

# patient level information
Y <- outcome %>% select(subject_id, los_icu) %>% mutate(Y = -los_icu) %>% select(-los_icu)
X0 <- baseline %>% select(subject_id, all_of(baseline_names))
df_patient <- merge(Y,X0, by = "subject_id")

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


# data set formatted for DTR function
MIMIC <- datagen

missing <- as.data.frame(apply(MIMIC, 2,function(x) sum(is.na(x))))

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

apply(MIMICimputed, 2,function(x) sum(is.na(x)))

save(MIMIC,
     MIMICimputed, file = "Data/mimic_imputed.Rdata")
