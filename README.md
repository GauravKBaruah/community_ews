# Repository for the paper: 

**"Community structure determines the predictability of population collapse"** 2022.Journal of Animal Ecology.
**Author: Gaurav Baruah, Arpat Ozgul, Christopher F Clements**


Contact: gbaruahecoevo@gmail.com 

## 01 mutualism data

This folder contains data from the simulations of the mutualism community module and R modelling scripts.
`01_functions.R` contains all the functions for simulation and EWS analyses. 
`05_Mutualism_constant_compvalues.R` is the R script that contains the function that simulates eco-evo dynamics and collapse of species in the community.
`plotting_mutualism.R` is the R script for plotting EWS and stability-resilience metrics and fitness benefits for the mutualism module.

`DATA_MUTUALISM_CONSTANT_COMPETITION.RData` is the data file that contains all the parameter values that leads to feasible communities including competition matrix, initial mean trait values etc.
`Mutualism_data_potential_results.RData` is the data file that constains the potential curves that could be plotted. Although this is not required to reproduce the figures of the main-text.

`Mutualism_data_potential_curve_ests_v3.RData` is the Rdata file that contains all the metrics estimated from 1-dimensional effective potential curve for the mutualism module.

`Mutualism_data_results_1.RData` contains the data of results for all the EWS metrics which could be used to reproduce EWS results of the mutualism module.

`Wentzell_potential.RData` is the data that contains stability-landscape metrics calculated for the quasi-potential function.

## 02 predator prey data

This folder contains R modelling scripts and simulated data for the predator-prey module.

`01_fweb_functions_analyses.R` is the script for all the functions used in modelling species collapse and EWS analyses. 

`food_web_constant_Competition.R` contains the function that simulates eco-evolutionary dynamics for species collapse within the predator-prey module. 

`sims_foodweb_constant_comp.R` script is the one that does all the replicate simulation and estimates all the potential metrics and EWS.

`plotting_foodweb_r.R` script is used to plot and reproduce figures for the food-web module.

 `FOOD_WEB_CONSTANT_COMPETITION.RData` s the data file that contains all the parameter values that leads to feasible communities including competition matrix, initial mean trait values etc.
 
 `foodweb_data_potential_curve_ests_v2.RData` is the Rdata file that contains all the metrics estimated from 1-dimensional effective potential curve for the mutualism module.
 
 `foodweb_data_potential_results_v2.RData` is the data file that constains the potential curves that could be plotted. Although this is not required to reproduce the figures of the main-text.
 
 `foodweb_data_results.RData` ontains the data of results for all the EWS metrics which could be used to reproduce EWS results of the predator-prey module.

 `Wentzell_potential_estimates.RData`is the data that contains stability-landscape metrics calculated for the quasi-potential function.



