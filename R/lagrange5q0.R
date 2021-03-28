#' lagrange5q0
#'
#' @description \code{lagrange5q0} predicts a series of 22 cumulative probabilities of dying (q(x)) and mortality rates (nMx) for the first 5 years of life, based on one or more mortality inputs between 0 and 5 years.
#'
#' @references \code{lagrange5q0} uses the  method of Langrage in order to solve the log-quadratic model presented in \url{https://repository.upenn.edu/psc_publications/54/}.
#' @author Andrea Verhulst, \email{verhulst@sas.upenn.edu}, Julio Romero \email{Julio.Romero-Prieto@lshtm.ac.uk}
#' @param data  mortality inputs formated with \code{format_data}.
#' @param k     constrained value of k (only when using a single mortality input).
#' @return \code{lagrange5q0} returns a list including three elements: the predicted values of q(5y) and k, and a data frame with predicted values of q(x) and nMx.
#'
#' The boundaries of the age intervals are given in days (1 year = 365.25 days, 1 month = 30.4375 days).
#'
#' Different age intervals are provived for the values of q(x), i.e. the cumulative probabilities of dying between 0 and age x, and for the values nMx,  i.e. the central death rates between x and x+n.
#'
#' @return Two types of 95% confidence intervals (CIs) are automatically estimated: (1)	CIs based on the patterns of prediction errors observed in the underlying data when using a single mortality input (k = 0), (2) CIs based on the uncertainty of k when minimizing the RMSE over a full series of 22 q(x).
#'
#'
#' @details \code{lagrange5q0} allows four types of mortality inputs:
#'
#' 1.	A single mortality input (either nqx or nMx) for any age interval. The value of k will be assumed equal to 0 (average outcome of the model) and the value of the mortality input will be matched exactly.
#'
#' 2.	A single mortality input (either nqx or nMx) for any age interval and a value of k. Both the mortality input and the value of k and will be matched exactly.
#'
#' 3.	Two mortality inputs (either nqx or nMx) for any age interval. The value of k will be estimated and the value of the two mortality inputs will be matched exactly.
#'
#' 4.	More than two mortality inputs (either nqx or nMx) for any age interval. One of them must be selected for matching. The value of k will be estimated, the selected input will be matched exactly, and the root mean square error will be minimized over the remaining mortality inputs.
#'
#'The computation of the root mean square error (RMSE) can be weighted. When minimizing a series of q(x), i.e. cumulative probabilities of dying starting at age 0, we recommend using weights proportional to the length of the last age interval (see fourth example below).
#'
#'Use \code{format_data} to prepare and verify the mortality inputs.
#'
#' @export
#'
#' @examples
#'#1. One input (k = 0): q(28d,5y) (Jordan 2015)
#'input <- format_data(
#'  rate      = 0.00804138,
#'  lower_age = 28,
#'  upper_age = 365.25*5,
#'  type      = "qx",
#'  sex       = "total")
#'
#'lagrange5q0(data = input)
#'
#'
#'#2. One input and k = 0.5: q(28d,5y) (Jordan 2015)
#'input <- format_data(
#'  rate      = 0.00804138,
#'  lower_age = 28,
#'  upper_age = 365.25*5,
#'  type      = "qx",
#'  sex       = "total")
#'
#'lagrange5q0(data = input, k = 0.5)
#'
#'
#'#3. Two inputs: q(0,1y) and q(1y,4y) (Belgium 1984)
#'input <- format_data(
#'  lower_age = c(0,365.25),
#'  upper_age = c(365.25,365.25*5),
#'  rate      = c(0.00996356, 0.00207788),
#'  type      = c("qx", "qx"),
#'  sex       = c("male", "male"))
#'
#'lagrange5q0(data = input)
#'
#'
#'#4. 22 inputs + 1 match: q(x) (Finland 1933)
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
#'
#'lagrange5q0(data = input)



lagrange5q0 <- function(data,k){

  if(F == ("format" %in% names(data)))             stop('Inputs were no formated with "format_data".')
  data <- data$input

  data(coef)
  pred <- subset(coef, sex == unique(data$sex))


  if(nrow(data) > 1){
  if(sum(c(0,7,14,21) %in% data$lower_age) == 0)    warning('No mortality input starting at age 0. Select a single mortality input instead for stable results.')
  }


  if(missing(k)) {
    par <- newton(data,pred,0)
  } else {
    if(nrow(data) > 1  )                stop('k cannot be constrained with more than one input')
    if(k < -1.1 | k > 1.5 )                warning('Input value of k extrapolated. k < -1.1 or k > 1.5.')
    par <- newton(data,pred,k)
  }

  if(par[1] == "error")                 stop("Model cannot find sensical solution.")

  if(pk < -1.1| k > 1.5)        warning('Predicted value of k extrapolated. k < -1.1 or k > 1.5.')
  if(exp(par$h) > 0.150)                warning('Predicted value of q(5y) extrapolated. q(5y) > 0.150.')

  pred$p_qx <- logquad(pred,par$h, par$k)
  if(is.unsorted(pred$p_qx))            stop('Model cannot converge to a solution.')

  pred$p_mx <- qx_to_mx(pred$p_qx)
  if(is.unsorted(rev(pred$p_mx)))       warning('Increase with age in the force of mortality (nMx). Prediction extrapolated.')




  check <- 0

  if(nrow(data[data$fit == "min" & data$type== "qx",]) ==  22){
    check          <- 1
    pred$qx        <- data$rate[data$fit == "min"]
    pred$wx        <- data$weight[data$fit == "min"]
    pred$p_qx_k0   <- logquad(pred,par$h,0)
    pred$ex_k0     <- log(pred$qx) - log(pred$p_qx_k0)
    var_k          <- sum(pred$wx*pred$ex_k0^2, na.rm = T)/sum(pred$wx*pred$vx^2, na.rm = T)-par$k^2

    pred$lower_p_qx  <- logquad(pred, par$h, (par$k+1.96*sqrt(var_k*22/21)))
    pred$upper_p_qx  <- logquad(pred, par$h, (par$k-1.96*sqrt(var_k*22/21)))
    pred$lower_p_mx  <- qx_to_mx(pred$lower_p_qx)
    pred$upper_p_mx  <- qx_to_mx(pred$upper_p_qx)

    if(is.unsorted(pred$lower_p_qx) |
       is.unsorted(pred$upper_p_qx))        warning('Model cannot converge to a solution for confidence interval.')

    if(is.unsorted(rev(pred$lower_p_mx)) |
       is.unsorted(rev(pred$upper_p_mx)))   warning('Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.')

    if(is.unsorted(pred$lower_p_qx) |
       is.unsorted(pred$upper_p_qx)){
      pred <- pred[,c("lower_age", "lower_age_m", "upper_age", "p_qx", "p_mx")]
      names(pred)[1] <- "low_age_q"
    }else{
      pred <- pred[,c("lower_age", "lower_age_m", "upper_age",
                      "p_qx", "lower_p_qx", "upper_p_qx",
                      "p_mx", "lower_p_mx", "upper_p_mx")]
      names(pred)[1] <- "low_age_q"
    }
  }


  if(par$k == 0 & nrow(data) == 1){
    check <- 2

    par1<- newton(data,pred,  .9327)
    par2<- newton(data,pred, -.6309)

    if(par1[1] == "error" | par2[1] == "error")      warning('Model cannot converge to a solution for confidence interval.')
    if(par1[1] == "error" | par2[1] == "error"){
      pred <- pred[,c("lower_age", "lower_age_m", "upper_age", "p_qx", "p_mx")]
      names(pred)[1] <- "low_age_q"

    }else{

      pred$lower_p_qx  <- logquad(pred,par1$h,  .9327)
      pred$upper_p_qx  <- logquad(pred,par2$h, -.6309)
      pred$lower_p_mx  <- qx_to_mx(pred$lower_p_qx)
      pred$upper_p_mx  <- qx_to_mx(pred$upper_p_qx)

      if(is.unsorted(pred$p.qx_lowerCI) |
         is.unsorted(pred$p.qx_upperCI))             warning('Model cannot converge to a solution for confidence interval.')

      if(is.unsorted(rev(pred$p.mx_lowerCI)) |
         is.unsorted(rev(pred$p.mx_upperCI)))        warning('Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.')

      if(is.unsorted(pred$p.qx_lowerCI) |
         is.unsorted(pred$p.qx_upperCI)){
        pred <- pred[,c("lower_age", "lower_age_m", "upper_age", "p_qx", "p_mx")]
        names(pred)[1] <- "low_age_q"
      }else{
        pred <- pred[,c("lower_age", "lower_age_m", "upper_age",
                        "p_qx", "lower_p_qx", "upper_p_qx",
                        "p_mx", "lower_p_mx", "upper_p_mx")]
        names(pred)[1] <- "low_age_q"
      }
    }
  }


  if(check == 0){
    pred <- pred[,c("lower_age", "lower_age_m", "upper_age", "p_qx", "p_mx")]
    names(pred)[1] <- "low_age_q"
  }

  return(list( "q(5y)"       = exp(par$h),
               "k"           = par$k,
               "predictions" = pred))

}
