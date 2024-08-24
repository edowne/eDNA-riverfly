# eDNA-riverfly
Data and analysis for MSc dissertation project and eDNA and Riverfly monitoring (2023-2024)

FILES: 
- "PROJECT Data Collection Form.xlsx": This is my main Excel for collecting field and lab data, including adding in Cq values from the qPCR data and adding in hoglouse count data from Riverfly dataset, so that it all sits together in one place.
- "COMP_Proj_Data.csv": This is a CSV version of the above dataset, which was used in the analyses in R.
- "QUESTIONNAIRE RESPONSES_Anonymised_15.07.24.xlsx": This is the data from the questionnaire to Riverfly participants. Any mentions of specific locations and names have been removed.
- "ANALYSIS_PART_1_LOD_LOQ.R": R script for calculating the Limit of Detection and Limit of Quantification, generating the standard/calibration curve and fitting a linear model to it to use for the unknown samples to estimate target DNA concentration based on their Cq.
- "ANALYSIS_PART_2_qPCR_vs_kick.R": R script for exploring the relationships between the qPCR data and the kick-sampling data, and comparing the three sampling methods.
- "ANALYSIS_PART_3_Env_vars_Model.R": R script for exploring the environmental variables and fitting a model to predict abundance (kick-sample counts) from qPCR/eDNA data.

