
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Modeling Age Pattern of Under-5 Mortality

<!-- badges: start -->
<!-- badges: end -->

**Working Paper**: [Guillot M., J. Romero Prieto, A. Verhulst, P.
Gerland. *Modeling Age Patterns of Under-5 Mortality: Results From a
Log-Quadratic Model Applied to High-Quality Vital Registration
Data*](https://repository.upenn.edu/psc_publications/54/).

## Overview

This repository contains the R package **logquad5q0** associated with
the above-mentioned article (accepted in the journal *Demography*).

The package uses the method of Lagrange to implement a log-quadratic
model able to estimate the age pattern of under-5 mortality by detailed
age. A variety of mortality inputs between 0 and 5 years can be used in
order to predict a series of 22 cumulative probabilities of dying (q(x))
and mortality rates (nMx) for the first 5 years of life.

Data and R code needed to replicate the coefficients of the
log-quadratic model are available on the website of the Under-5
Mortality Database (U5MD):
<https://web.sas.upenn.edu/global-age-patterns-under-five-mortality/>

To report issues or suggest improvements, please contact
verhulst@sas.upenn.edu.

#### Authors of the package

-   **Andrea Verhulst** - *Development of the R code and package.* -
    [web page](https://www.pop.upenn.edu/bio/andrea-verhulst)

-   **Julio Romero** - *Analytical development of the method of Langrage
    for solving the log-quadratic model.* - [web
    page](https://www.lshtm.ac.uk/aboutus/people/romero-prieto.julio)

## Installation

You can install the released version of logquad5q0 with:

``` r
install.packages("devtools")
library(devtools)

install_github(repo = "verhulsta/logquad5q0")
library(logquad5q0)
```

## Main function: lagrange5q0

The function **lagrange5q0** predicts a series of 22 cumulative
probabilities of dying (q(x)) and mortality rates (nMx) for the first 5
years of life, based on one or more mortality inputs between 0 and 5
years.

This function allows four types of mortality inputs:

1.  A single mortality input (either nMx and nqx) for any age interval.
    The value of k will be assumed equal to 0 (average outcome of the
    model) and the value of the mortality input will be matched exactly.

2.  A single mortality input (either nMx and nqx) for any age interval
    and a value of k. Both the mortality input and the value of k and
    will be matched exactly.

3.  Two mortality inputs (either nMx and nqx) for any age interval. The
    value of k will be estimated and the value of the two mortality
    inputs will be matched exactly.

4.  More than two mortality inputs (either nMx and nqx) for any age
    interval. One of them must be selected for matching. The value of k
    will be estimated, the selected input will be matched exactly, and
    the root mean square error will be minimized over the remaining
    mortality inputs.

The computation of the root mean square error (RMSE) can be weighted.
When minimizing a series of q(x), i.e. cumulative probabilities starting
at age 0, we recommend using weights proportional to the length of the
last age interval (see fourth example below). This approach was used to
produce the results of the paper.

\_

Use the function **format\_data** to prepare and verify the mortality
inputs.

Age intervals must be defined with two exact ages *in days* (1 year =
365.25 days, 1 month = 30.4375 days). The log-quadradic model deals with
the following age breakdowns:

| Age       | in Days   |
|:----------|-----------|
| 0         | 0         |
| 7 days    | 7         |
| 14 days   | 14        |
| 21 days   | 21        |
| 28 days   | 28        |
| 2 months  | 60.8750   |
| 3 months  | 91.3125   |
| 4 months  | 121.7500  |
| 5 months  | 152.1875  |
| 6 months  | 182.6250  |
| 7 months  | 213.0625  |
| 8 months  | 243.5000  |
| 9 months  | 273.9375  |
| 10 months | 304.3750  |
| 11 months | 334.8125  |
| 12 months | 365.2500  |
| 15 months | 456.5625  |
| 18 months | 547.8750  |
| 21 months | 639.1875  |
| 2 years   | 730.5000  |
| 3 years   | 1095.7500 |
| 4 years   | 1461.0000 |
| 5 years   | 1826.2500 |

\_

The function lagrange5q0 returns a list including:

-   The predicted value of q(5y).
-   The predicted value of k.
-   A data frame with predicted values of q(x) and nMx.

Different age intervals are provived for the values of q(x), i.e. the
cumulative probabilities of dying between 0 and age x, and for the
values nMx, i.e. the central death rates between x and x+n.

Two types of 95% confidence intervals (CIs) are automatically provided:

-   CIs based on the patterns of prediction errors observed in the
    underlying data when using a single mortality input (k = 0).
-   CIs based on the uncertainty of k when minimizing the RMSE over a
    full series of 22 q(x).

#### Caveat

The log-quadratic model is based on a set of vital records from Western
countries. Predictions based on values of q(5y) above 0.150 and on value
of k below -1.1 or above +1.5 are extrapolations. A warning message will
be produced when such cases occur. Beyond these limits, it also might be
impossible to converge to a solution. An error messages will be produced
in such case.

## Examples

``` r
#1. One input (k = 0): q(28d,5y) (Jordan 2015)
df <- format_data(
   rate      = 0.00804138,
   lower_age = 28,
   upper_age = 365.25*5,
   type      = "qx",
   fit       = "match",
   sex       = "total")

lagrange5q0(df)


#2. One input and k = 0.5: q(28d,5y) (Jordan 2015)
df <- format_data(
   rate      = 0.00804138,
   lower_age = 28,
   upper_age = 365.25*5,
   type      = "qx",
   fit       = "match",
   sex       = "total")

lagrange5q0(df, k = 0.5)


#3. Two inputs: q(0,1y) and q(1y,4y) (Belgium 1984)
df <- format_data(
  lower_age = c(0,365.25),
  upper_age = c(365.25,365.25*5),
  rate      = c(0.00996356, 0.00207788),
  type      = c("qx", "qx"),
  sex       = c("male", "male"))

lagrange5q0(df)


#4. 22 inputs + 1 match: q(x) (Finland 1933)
data(fin1933)
fin1933$weight <- c(fin1933$n[1:22]/(365.25*5), NA) 
df <- format_data(
  lower_age = fin_1933$lower_age,
  upper_age = fin_1933$upper_age,
  rate      = fin_1933$rate,
  type      = fin_1933$type,
  sex       = fin_1933$sex,
  fit       = fin_1933$fit,
  weight    = fin_1933$weight)

lagrange5q0(df)
```
