#' format_data
#'
#' @description \code{format_data} prepares and verifies mortality inputs for the implementation of \code{lagrange5q0}.

#' @param lower_age Lower boundary of the age interval in days for each input (see details below for values).
#' @param upper_age Upper boundary of the age interval in days for each input (see details below for values).
#' @param rate  Mortality rate for each input.
#' @param type  Character: 'qx' or 'mx'. Type of rate for each input.
#' @param sex   Character: 'female', 'male', or 'total'. Same for all inputs.
#' @param fit   Character: 'match' or 'min'. Needed only when there are more than 2 inputs. Select the input to be matched using 'match'. Use 'min' for all the other inputs.
#' @param weigth Optional when there are more than 2 inputs: weight for inputs whose root mean square error is minimized. Use 0 or NA for the matched input. The sum of weights must be equal to 1.
#'
#' @return \code{format_data} returns a data frame with the formated mortality inputs ready for the implementation of \code{lagrange5q0}.
#'
#' @details
#'
#' Possible boudaries for the age intervals (1 year = 365.25 days, 1 month = 30.4375 days):
#'
#'  | Age           | in Days       |
#'  | :------------ |---------------|
#'  |0              | 0             |
#'  |7 days         | 7             |
#'  |14 days        | 14            |
#'  |21 days        | 21            |
#'  |28 days        | 28            |
#'  |2 months       | 60.8750       |
#'  |3 months       | 91.3125       |
#'  |4 months       | 121.7500      |
#'  |5 months       | 152.1875      |
#'  |6 months       | 182.6250      |
#'  |7 months       | 213.0625      |
#'  |8 months       | 243.5000      |
#'  |9 months       | 273.9375      |
#'  |10 months      | 304.3750      |
#'  |11 months      | 334.8125      |
#'  |12 months      | 365.2500      |
#'  |15 months      | 456.5625      |
#'  |18 months      | 547.8750      |
#'  |21 months      | 639.1875      |
#'  |2 years        | 730.5000      |
#'  |3 years        | 1095.7500     |
#'  |4 years        | 1461.0000     |
#'  |5 years        | 1826.2500     |
#'
#'
#'
#'
#' @export
#'
#' @examples
#'#One input: q(28d,5y) (Jordan 2015)
#'input <- format_data(
#'  rate      = 0.00804138,
#'  lower_age = 28,
#'  upper_age = 365.25*5,
#'  type      = "qx",
#'  fit       = "match",
#'  sex       = "total")
#'
#'#Two inputs: q(0,1y) and q(1y,4y) (Belgium 1984)
#'input <- format_data(
#'  lower_age = c(0,365.25),
#'  upper_age = c(365.25,365.25*5),
#'  rate      = c(0.00996356, 0.00207788),
#'  type      = c("qx", "qx"),
#'  sex       = c("male", "male"))
#'
#'
#'#22 inputs + 1 match: q(x) (Finland 1933)
#'data(fin1933)
#'fin1933$weight <- c(fin1933$n[1:22]/(365.25*5), NA)
#'input <- format_data(
#'  lower_age = fin1933$lower_age,
#'  upper_age = fin1933$upper_age,
#'  rate      = fin1933$rate,
#'  type      = fin1933$type,
#'  sex       = fin1933$sex,
#'  fit       = fin1933$fit,
#'  weight    = fin1933$weight)


format_data <- function(lower_age, upper_age, rate, type, sex, fit, weight){

  if(missing(lower_age))         stop('"lower_age" is missing.')
  if(missing(upper_age))         stop('"upper_age" is missing.')
  if(missing(rate))              stop('"rate"      is missing.')
  if(missing(type))              stop('"type"      is missing.')
  if(missing(sex))               stop('"sex"       is missing.')

  if(length(unique(c(length(lower_age),
                     length(upper_age),
                     length(rate),
                     length(type),
                     length(sex)))) !=  1)      stop('Variables must have the same length.')

  if(missing(fit) == F){
    if(length(lower_age) != length(fit))        stop('Variables must have the same length.')
  }

  if(missing(weight) == F){
    if(length(lower_age) != length(weight))     stop('Variables must have the same length.')
  }

  if(!is.numeric(lower_age))    stop('"lower_age" must be numeric.')
  if(!is.numeric(upper_age))    stop('"upper_age" must be numeric.')
  if(!is.numeric(rate))         stop('"rates" must be numeric.')
  if(!is.character(type))       stop('"type" must be a character string.')
  if(!is.character(sex))        stop('"sex" must be a character string.')

  if(length(unique(sex)) != 1)  stop('"sex" must be identical for all inputs.')

  if(missing(fit) == F){
    if(!is.character(fit))      stop('"fit" must be a character string.')
  }

  if((F %in% (lower_age  %in%       c(0.0000,    7.0000,    14.0000,   21.0000,   28.0000,
                                      60.8750,   91.3125,   121.7500,  152.1875,  182.6250,
                                      213.0625,  243.5000,  273.9375,  304.3750,  334.8125,
                                      365.2500,  456.5625,  547.8750,  639.1875,  730.5000,
                                      1095.7500, 1461.0000))))                               stop('Wrong value for "lower_age".')
  if((F %in% (upper_age  %in%       c(7.0000,    14.0000,   21.0000,   28.0000,
                                      60.8750,   91.3125,   121.7500,  152.1875,  182.6250,
                                      213.0625,  243.5000,  273.9375,  304.3750,  334.8125,
                                      365.2500,  456.5625,  547.8750,  639.1875,  730.5000,
                                      1095.7500, 1461.0000, 1826.2500))))                   stop('Wrong value for "upper_age".')
  if(F %in% (lower_age < upper_age))                                                        stop('"lower_age"   >  "upper_age".')
  if(T %in% (rate <= 0))                                                                    stop('"rates" must be > 0.')
  if(F %in% (type   %in%   c("qx",  "mx")))                                                 stop('Wrong value for "type".')
  if(F %in% (sex    %in%   c("female",  "male", "total")))                                  stop('Wrong value for "sex".')

  if(length(lower_age) == 1 | length(lower_age)  == 2){

    if(missing(fit) == F){
      if((F %in%   (fit   %in%   c("min", "match"))))  warning('Wrong value for "fit". "match" used instead')
      if("min" %in% fit)                               warning('Wrong value for "fit". "match" used instead')
      fit <- "match"
    }else{
      fit <- "match"
    }

    if(missing(weight) == F)                           warning('"weight" not used with matching.')
  }

  df <-data.frame(
    "lower_age" = lower_age,
    "upper_age" = upper_age,
    "rate"      = rate,
    "type"      = type,
    "sex"       = sex,
    "fit"       = fit)

  if(length(lower_age) > 2){

    if(missing(fit))                                            stop('"fit" is missing.')

    if(F %in% (df$fit  %in%  c("min", "match")))                stop('Wrong value for fit.')
    if(F %in% ("match" %in% df$fit))                            stop('"match" is missing in "fit".')
    if(sum(    "match" == df$fit) > 1)                          stop('Only one input can be matched.')
    if(F %in% ("min"   %in% df$fit))                            stop('"min" is missing in "fit".')


    if(missing(weight) == F){
      df$weight <- weight
      if(is.character(df$weight[df$fit == "min"]))              stop('"weight" must be numeric.')
      if(!is.na(df$weight[df$fit == "match"]))                  warning('Weight for matched input not used. NA used instead')
      df$weight[df$fit == "match"] <- NA
      if(sum(df$weight[df$fit == "min"], na.rm = T) != 1)       stop('Sum of weights different from 1.')
    }else{weight <- rep(1/length(lower_age),length(lower_age))}

  }

  if(sum(duplicated(df[,c("lower_age", "upper_age", "type")])) !=
     sum(duplicated(df[,c("lower_age", "upper_age", "type", "rate")])))  stop('Different rates for a same age group.')

  return(list("input" = df,
              "format" = "checked"))


}



