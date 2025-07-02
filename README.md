# Can k-fold cross validation improve survival model selection within health technology assessments? An exploratory study  

This repository includes code for the article on the use of k-fold cross validation for the selection of survival curves for health technology assessments.

## Authors

Inigo Bermejo, Data Science Institute (DSI), Hasselt University, Diepenbeek, Belgium

Sabine Grimm,  Klinische Epidemiologie en Medical Technology Assessment (KEMTA), Maastricht University Medical Centre (MUMC+), Maastricht, Limburg, Netherlands

## Background

The selection of survival models for informing economic evaluations of innovative therapies with limited long-term data traditionally relies on model selection algorithms applied in the full trial data. Models selected based on full trial data might underperform in the target population due to overfitting. K-fold cross validation (CV) is commonly used in machine learning for better estimation of fit in unseen data. We explore whether k-fold CV improves model selection.

# Instructions

This repository does not contain the data used, but most datasets are available through the survival package. The other datasets are accessible through SEER and other sources linked in the script. Make sure the datasets are placed in the locations the scripts expects them to be.

Run the script in SurvivalFitCrossValidation.R from beginning to end. This will create a number of results RDS files as well as PNG files with the Kaplan-Meier plots included in the paper. 

Running  the script in Result_processing.R will process these RDS files and generate the result tables shown in the article. 

