# STAT 638 Bayesian Project

## Group Members

- Aysegul Buber
- Christopher Bowry
- Aditya Hoskere Rajiv

## Goal

Characterize temporal patterns in daily household electricity consumption and
develop a Bayesian analysis that can be used to predict future daily
electricity demand.

## Project 17: Bayesian Forecasting of Household Electricity Demand

### Background

Residential electricity demand varies over time because of household activity,
recurring weekly patterns, seasonal changes, and other factors. Accurate
forecasting of electricity demand is important for energy planning and
management. In this project, we analyze household electricity-consumption data
and quantify uncertainty in future electricity demand.

### Data

This project uses the **UCI Individual Household Electric Power Consumption**
dataset. The data contain one-minute measurements of household electricity use
from December 2006 through November 2010. These measurements are aggregated to
obtain a daily electricity-consumption outcome for the analysis.

**Data source:** [UCI Machine Learning Repository — Individual Household Electric Power Consumption](https://archive.ics.uci.edu/dataset/235/individual+household+electric+power+consumption)

### Scientific Objective

Characterize temporal patterns in daily household electricity consumption and
develop a Bayesian analysis that can be used to predict future daily
electricity demand.

### Seasonal Component

The model must account for annual seasonality using sine and cosine terms with
a period of 365 days. The first harmonic is

$$
\sin\left(\frac{2\pi t}{365}\right)
\quad\text{and}\quad
\cos\left(\frac{2\pi t}{365}\right).
$$

We will investigate whether including the second harmonic improves the
representation of the seasonal pattern:

$$
\sin\left(\frac{4\pi t}{365}\right)
\quad\text{and}\quad
\cos\left(\frac{4\pi t}{365}\right).
$$

### Scientific Questions

1. How does daily household electricity consumption vary over the study
   period?
2. Is there evidence of a long-term increase or decrease in electricity
   consumption?
3. What annual seasonal pattern is present, and does the second harmonic
   provide an important improvement over the first-harmonic representation?
4. Are there systematic weekly patterns in electricity demand?
5. How accurately can future daily electricity consumption be predicted, and
   how does predictive uncertainty change with the forecasting horizon?
6. Are there periods of unusually high or low electricity consumption that are
   not adequately explained by the fitted model?

### Analysis Expectations

Develop and justify an appropriate Bayesian model that addresses the scientific
questions above. Clearly describe the conclusions in terms of household
electricity demand and its uncertainty.
