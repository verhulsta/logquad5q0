
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Modeling Age Patterns of Under-5 Mortality

<!-- badges: start -->
<!-- badges: end -->

**Papers**:

[Guillot M., J. Romero Prieto, A. Verhulst, P. Gerland. Modeling Age
Patterns of Under-5 Mortality: Results From a Log-Quadratic Model
Applied to High-Quality Vital Registration Data.
*Demography*](https://doi.org/10.1215/00703370-9709538).

[J. Romero Prieto J. , A. Verhulst, and M. Guillot. Estimating 1a0 and
4a1 in a Life Table: A Model Approach Based on Newly Collected Data.
*Demography*](https://doi.org/10.1215/00703370-11330227).

## Overview

This repository contains the R package **logquad5q0** associated with
the above-mentioned article published in the journal *Demography*.

The package uses the method of Lagrange to implement a log-quadratic
model able to estimate the age pattern of under-5 mortality by detailed
age. A variety of mortality inputs between 0 and 5 years can be used in
order to predict a series of 22 cumulative probabilities of dying (q(x))
and mortality rates (nMx) for the first 5 years of life.

Based on these outputs, the package also computes the average number of
years lived in the age intervals 0-1, 1-4 and 0-5 (i.e., 1a0, 4a1, and
5a0).

Data and R code needed to replicate the coefficients of the
log-quadratic model are available on the website of the Under-5
Mortality Database (U5MD):
<https://web.sas.upenn.edu/global-age-patterns-under-five-mortality/>

<span style="color:red"> BETA VERSION: The package includes new
coefficients derived from Demographic and Health Surveys (DHS) collected
between 1985 and 2022 in sub-Saharan Africa and south Asia. These
coefficients will be updated with the most recent DHS available before
publication of the paper. </span>

To report issues or suggest improvements, please contact
<andrea.verhulst@ined.fr>

#### Authors of the package

- **Andrea Verhulst** - *Development of the R code and package.* - [web
  page](https://www.ined.fr/en/research/researchers/Verhulst+Andrea)

- **Julio Romero** - *Analytical development of the method of Langrage
  for solving the log-quadratic model (under MATLAB).* - [web
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

1.  A single mortality input (either nqx or nMx) for any age interval.
    The value of k will be assumed equal to 0 (average outcome of the
    model) and the value of the mortality input will be matched exactly.

2.  A single mortality input (either nqx or nMx) for any age interval
    and a value of k. Both the mortality input and the value of k and
    will be matched exactly.

3.  Two mortality inputs (either nqx or nMx) for any age interval. For
    the purpose of estimating nax–the average number of years lived in
    the age interval x,x+n–with high precision, the proportion of infant
    deaths before 28 days or 3 months (z(28d) or z(3m)) can be used as
    second input along with another mortality input (either nqx or nMx)
    starting at age 0. The two mortality inputs will be matched exactly
    and the value of k will be estimated.

4.  More than two mortality inputs (either nqx or nMx) for any age
    interval. One of them must be selected for matching. The selected
    input will be matched exactly, the root mean square error will be
    minimized over the remaining mortality inputs, and the value of k
    will be estimated.

The computation of the root mean square error (RMSE) can be weighted.
When minimizing a series of q(x), i.e. cumulative probabilities starting
at age 0, we recommend using weights proportional to the length of the
last age interval (see example 5 below). This approach was used to
produce the results of the paper.

\_

Use the function **format_data** to prepare and verify the mortality
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

- The predicted value of q(5y).
- The predicted value of k.
- A data frame with predicted values of q(x) and nMx.
- The predicted value of 1a0, 4a1, and 5a0.

Different age intervals are provived for the values of q(x), i.e. the
cumulative probabilities of dying between 0 and age x, and for the
values nMx, i.e. the central death rates between x and x+n.

Two types of 95% confidence intervals (CIs) are automatically provided
when:

- CIs are based on the patterns of prediction errors observed in the
  underlying data when using a single mortality input (k = 0).
- CIs are based on the uncertainty of k when minimizing the RMSE over a
  full series of 22 q(x).

#### Caveat

The log-quadratic model is based on a set of vital records from Western
countries. Predictions based on values of q(5y) above 0.150 and on
values of k below -1.1 or above +1.5 are extrapolations. A warning
message will be produced when such cases occur. Beyond these limits, it
also might be impossible to converge to a solution. An error message
will be produced in such case.

<span style="color:red">

#### BETA VERSION

To use the new coefficients for sub-Saharan Africa and south Asia,
specify the ‘B’ model in the function **lagrange5q0** (see example 6
below).

In the case scenario of having a single mortality input and no
information on the value of k, the user can specify a region (‘Eastern
Africa’, ‘Middle Africa’, ‘Western Africa’, or ‘South Asia’) in order to
benefit from a regional prior of k instead of assuming k = 0 (see
example 6 below).

</span>

## Examples

``` r

#1. One input (k = 0): q(28d,5y) (VR Jordan 2015)
input <- format_data(
   rate      = 0.00804138,
   lower_age = 28,
   upper_age = 365.25*5,
   type      = "qx",
   sex       = "total")

lagrange5q0(data = input)


#2. One input and k = 0.5: q(28d,5y) (VR Jordan 2015)
input <- format_data(
   rate      = 0.00804138,
   lower_age = 28,
   upper_age = 365.25*5,
   type      = "qx",
   sex       = "total")

lagrange5q0(data = input, k = 0.5)


#3. Two inputs: q(0,1y) and q(1y,4y) (VR Belgium 1984)
input <- format_data(
  lower_age = c(0,365.25),
  upper_age = c(365.25,365.25*5),
  rate      = c(0.00996356, 0.00207788),
  type      = c("qx", "qx"),
  sex       = c("male", "male"))

lagrange5q0(data = input)


#4. Two inputs: M(0,1y) and z(28d) (VR Australia 1935)
input <-  format_data(
  lower_age = c(0,0),
  upper_age = c(365.25, 28),
  rate      = c(0.04196866, 0.68614291),
  type      = c("mx", "zx"),
  sex       = c("total", "total"))

lagrange5q0(data = input)


#5. 22 inputs + 1 match: q(x) (VR Finland 1933)
data(fin1933)
fin1933$weight <- c(fin1933$n[1:22]/(365.25*5), NA) 
input <- format_data(
  lower_age = fin1933$lower_age,
  upper_age = fin1933$upper_age,
  rate      = fin1933$rate,
  type      = fin1933$type,
  sex       = fin1933$sex,
  fit       = fin1933$fit,
  weight    = fin1933$weight)

lagrange5q0(data = input)


#6. One input (Model B & regional k) : q(5y) (DHS DR Congo 2013-14)
input <- format_data(
   rate      = 0.10924,
   lower_age = 0,
   upper_age = 365.25*5,
   type      = "qx",
   sex       = "total")

lagrange5q0(data = input, model = 'B', region = "Middle Africa")
```
