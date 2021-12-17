# Wisdom of Crowds: Hypothesis Test Algorithm

This repository contains data and R code for replicating simulations, empirical analysis, power analysis, 
tables and figures described in the working paper "A Hypothesis Test Algorithm for Determining 
When Weighting Individual Judgments Reliably Improves Collective Accuracy or Just Adds Noise" 
by Shu Huang, Stephen B. Broomell, and Russell Golman. 

## How to use
1. Clone this repository to your local computer
2. For each script to run, be sure to change the working directory in each script (searching "setwd()")

## 0algorithm
1. The file "weighting functions.R" contains all related weighting models mentioned in the paper. 
2. The file "test algorithm.R" contains the hypothesis test algorithm proposed in the working paper. 
3. The file "cross validation.R" contains the commonly used cross validation method. 

## 1simulation
1. The file "Sim_EqualCase.R" contains the sample code for the simulation with equal-weighting true states. 
2a. The file "Sim_SpfBoots.R" contains the sample code for the simulation with bootstrapped true states based on two SPF data sets.
2b. The file "Sim_SpfBoots_OW.R" also contains the sample code for the simulation with bootstrapped true states based on two SPF data sets, but only for the optimal weighting methods without or with the non-negative constraint.
3. The file "Sim_CovidEst.R" contains the sample code for the simulation with estimated true states from covid-19 new cases forecasts.
4. The files "EqualStates.csv", "InfStates.csv", "UnempStates.csv", "InfStates_OW.csv", "UnempStates_OW.csv", and "CovidStates.csv" record the simulation results for the corresponding true states cases. Those results include (for each given sample size and for each weighting method) the p-values, cross validation errors, out-of-sample MSEs and the t-test results for comparing the out-of-sample MSEs between weighted averages and the simple average. 

## 2empirical_analysis
1. The file "EmpiricalAnalysis_inf.R" contains the sample code for analyzing the inflation rate forecasts from ECB SPF by using the data file "inf_filtered.csv". 
2. The file "EmpiricalAnalysis_unemp.R" contains the sample code for analyzing the unemployment rate forecasts from ECB SPF by using the data file "unemp_filtered.csv".
3. The file "EmpiricalAnalysis_cases.R" contains the sample code for analyzing the covid-19 new cases forecasts from CDC US by using the data file "all_data_case_county.csv". The file "all_data_case_county.csv" is too large to upload on Github. Please contact the author for requesting the dataset if necessary. 
4. The file "EmpiricalAnalysis_deaths.R" contains the sample code for analyzing the covid-19 cumulative deaths forecasts from CDC US by using the data file "all_data_death.csv".
5. The file "EmpiricalAnalysis_keck.R" contains the sample code for analyzing the experimental data from Keck and Tang (2020) by using the data file "KeckDataStudy1".
6. The series of files titled with "res_" record all empirical analysis results. Those results include the p-values and out-of-sample MSEs for the OW and CWM weighting methods for each given sample size. 

## 3power_analysis
1. The file "PowerAnalysis.R" contains the sample code for doing the power analysis in this paper. 
2. The files "power_inf.csv" and "power_unemp.csv" are the recorded p-values for all weighting methods with different given sample sizes. 

## 4results_analysis
1. The file "ResultsAnalysis.R" contains the sample code for loading all results from simulation, empirical analysis and power analysis and then making the tables and figures in the paper. 
