# STAT_638_Bayes_Project  

### Group Members: Aysegul Buber, Christopher Bowry, and Aditya Hoskere Rajiv

### Goal: Characterizing temporal patterns in daily household electricity consumption and developing a Bayesian analysis that can be used to predict future daily electricity demand.  
The dataset is the UCI Individual Household Electric Power Consumption dataset. 

Project 17: Bayesian Forecasting of Household Electricity Demand
Background. Residential electricity demand varies over time because of household activity, recurring
weekly patterns, seasonal changes, and other factors. Accurate forecasting of electricity demand is
important for energy planning and management. In this project, you will analyze household electricity-
consumption data and quantify uncertainty in future electricity demand.
Data. Use the UCI Individual Household Electric Power Consumption dataset. The data contain one-
minute measurements of household electricity use from December 2006 through November 2010.
Aggregate the measurements to obtain a daily electricity-consumption outcome for the analysis.
Data source. UCI Machine Learning Repository: Individual Household Electric Power Consumption.
Scientific objective. Characterize temporal patterns in daily household electricity consumption and
develop a Bayesian analysis that can be used to predict future daily electricity demand.
Seasonal component. Your model must account for annual seasonality using sine and cosine terms
with a period of 365 days. Begin with the first harmonic:
sin(2πt/365) and cos(2πt/365).
Investigate whether including the second harmonic improves the representation of the seasonal
pattern:
sin(4πt/365) and cos(4πt/365).
Scientific questions.
How does daily household electricity consumption vary over the study period?
Is there evidence of a long-term increase or decrease in electricity consumption?
What annual seasonal pattern is present, and does the second harmonic provide an important
improvement over the first-harmonic representation?
Are there systematic weekly patterns in electricity demand?
How accurately can future daily electricity consumption be predicted, and how does predictive
uncertainty change with the forecasting horizon?
Are there periods of unusually high or low electricity consumption that are not adequately explained
by the fitted model?
Analysis expectation. Develop and justify an appropriate Bayesian model that addresses the scientific
questions above. Clearly describe the conclusions in terms of household electricity demand and its
uncertainty.
