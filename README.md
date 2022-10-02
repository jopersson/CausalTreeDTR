# Learning Optimal Dynamic Treatment Regimes Using Causal Tree Methods with Application to Medicine
This is an anonymous git repository for the paper "Learning Optimal Dynamic Treatment Regimes Using Causal Tree Methods with Application to Medicine".

## Scripts
The following scripts are used to reproduce the results:
- R/data_generating_process.R (data generating processes for two scenarios: linear_interactions, nonlinear_interactions)
- R/evaluation.R (functions to create plots and tables for evaluation)
- R/DTR.R (main function to fit, predict and evaluate DTRs)
- R/DTR_fit.R (test script to fit a DTR once)
- R/DTR_simulation.R (simulation and evaluation of DTR fits with different parameter settings)
- R/DTR_hyperparameter.R (simulation and evaluation of DTR fits with different parameter settings)
- Data/mimic_preprocess.py (create and save single data frames from the "all_hourly_data.h5"-file)
- Data/mimic_preprocess.R (select treatments, outcomes and patient characteristics)
- Data/mimic_prepare.R (prepare a single data frame that our function can deal with)
- Data/DTR_mimic.R (reproduce the results with preprocessed MIMIC-data)

The first six scripts are self contained. We adapt the structure of our DTR function from the [DTRreg function](https://github.com/cran/DTRreg/blob/master/R/DTRreg.R), introduced by Wallace et al. (2017). The last four use the pre-processed MIMIC-III data set (Johnson et al., 2016). The we use the pre-processing pipeline provided by Wang et. al. (2020). See the references below.

## Requirements
### python 3.6 
h5py, csv, numpy
### R 3.6.3
causalTree, grf, rpart, FNN, dplyr, tidyr, ggplot2, mmtable2, gt, forcats, purrr

## Reproducing Simulation Results
Run the script DTR_simulation.R. You can specify approaches and scenarios of interest. The script contains synthetic data simulation and code that produces and saves the results reported. It outputs tables to the directory "~/Results/".

## Reproducing Hyperparameter Tuning Simulation Results
Run the script DTR_hyperparameter.R to learn about the performance of DTR-CT and DTR-CF over varying minimum number of observations required per leaf. You can specify approaches and scenarios of interest. The script contains synthetic data simulation and code that produces and saves the results reported. It outputs plots to the directory "~/Results/".

## Reproducing MIMIC-III Results
Run the script DTR_mimic.R. The script requires the pre-processed MIMIC-III data as input (see references below). It outputs an illustration of a causal tree and a variable importance plot of a causal forest, indicating the time stage of interest. 

MIMIC-III is a freely accessible database. However, access must be requested at https://physionet.org/content/mimiciii/1.4/. When MIMIC-III access is granted, the pre-processed data by Wang et. al. (2020) is accessible with instructions in the respective paper. We use this data set to study the effect of mechanical ventilation and vasopressin on the outcome variable of length of stay in the intensive care unit (ICU).

We extract 8059 patients with 10 time steps each. As baseline covariates, we consider age, gender, race, insurance, admission type, and first care unit. As time-vayring covariates, we use heart rate, diastolic blood pressure, mean blood pressure, oxygen saturation, respiratory rate, systolic blood pressure, and temperature. 

## References
Johnson, A. E., Pollard, T. J., Shen, L., Li-Wei, H. L., Feng, M., Ghassemi, M., ... & Mark, R. G. (2016). MIMIC-III, a freely accessible critical care database. Scientific data, 3(1), 1-9.

Wang, S., McDermott, M. B., Chauhan, G., Ghassemi, M., Hughes, M. C., & Naumann, T. (2020, April). MIMIC-extract: A data extraction, preprocessing, and representation pipeline for MIMIC-III. In Proceedings of the ACM Conference on Health, Inference, and Learning (pp. 222-235).

Wallace, M. P., Moodie, E. E., & Stephens, D. A. (2017). Dynamic treatment regimen estimation via regression-based techniques: Introducing r package dtrreg. Journal of Statistical Software, 80(1), 1-20.
