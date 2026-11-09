# population-trends-dml-iml

Caroline D. Cappello   
Nov 2026  
  

code provided for peer review


## scripts

__run_dml_hpc.R:__ Runs the double machine learning (DML) trend model. Saves out model results and intermediate products to assess model performance.  

__model_scenarios.R:__ Script is read into run_dml_hpc.R and provides lists of specific environmental variables to be included in the propensity, conditional mean (i.e., marginal), and treatment models

__run_gam.R:__ Runs the IML GAM analyses using output and env variables from the DML trend model. 