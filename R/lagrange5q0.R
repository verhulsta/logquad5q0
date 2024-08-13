#' lagrange5q0
#'
#' @description \code{lagrange5q0} predicts a series of 22 cumulative probabilities of dying (q(x)) and mortality rates (nMx) for the first 5 years of life, based on one or more mortality inputs between 0 and 5 years.
#'
#' @references \code{lagrange5q0} uses the  method of Langrage in order to solve the log-quadratic model presented in \url{https://doi.org/10.1215/00703370-9709538}.
#' @author Andrea Verhulst, \email{andrea.verhulst@ined.fr}, Julio Romero \email{Julio.Romero-Prieto@lshtm.ac.uk}
#' @param data    mortality inputs formated with \code{format_data}.
#' @param k       constrained value of k (only when using a single mortality input).
#' @param model   Character: 'A' (default) or 'B'.
#' @param region  Character: 'Eastern Africa', 'Middle Africa', 'Western Africa', or 'South Asia' (only for Model B).
#' @return \code{lagrange5q0} returns a list including six elements: the predicted values of q(5y) and k, a data frame with predicted values of q(x) and nMx, and the predicted value of 1a0, 4a1, and 5a0.
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
#' 3.	Two mortality inputs (either nqx or nMx) for any age interval. For the purpose of estimating nax--the average number of years lived in the age interval x,x+n--with high precision, the proportion of infant deaths before 28 days or 3 months (z(28d) or z(3m)) can be used as second input along with another mortality input (either nqx or nMx) starting at age 0. The two mortality inputs will be matched exactly and the value of k will be estimated.
#'
#' 4.	More than two mortality inputs (either nqx or nMx) for any age interval. One of them must be selected for matching. The selected input will be matched exactly, the root mean square error will be minimized over the remaining mortality inputs, and the value of k will be estimated.
#'
#'The computation of the root mean square error (RMSE) can be weighted. When minimizing a series of q(x), i.e. cumulative probabilities of dying starting at age 0, we recommend using weights proportional to the length of the last age interval (see example 5 below).
#'
#'Use \code{format_data} to prepare and verify the mortality inputs.
#'
#'
#' BETA VERSION: To use coefficients for sub-Saharan Africa and south Asia, specify 'B' for the parameter 'model'. In the case scenario of having a single mortality input and no information on the value of k, the user can specify the parameter 'region' ('Eastern Africa', 'Middle Africa', 'Western Africa', or 'South Asia') in order to benefit from a regional prior of k instead of assuming k = 0 (see example 6 below)
#'
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
#'#4. Two inputs: M(0,1y) and z(28d) (Australia 1935)
#"input <-  format_data(
#'  lower_age = c(0,0),
#'  upper_age = c(365.25, 28),
#'  rate      = c(0.04196866, 0.68614291),
#'  type      = c("mx", "zx"),
#'  sex       = c("total", "total"))
#'
#'lagrange5q0(data = input)
#'
#'#5. 22 inputs + 1 match: q(x) (Finland 1933)
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
#'
#'#6. One input (Model B & regional k) : q(5y) (DHS DR Congo 2013-14)
#'input <- format_data(
#'  rate      = 0.10924,
#'  lower_age = 0,
#'  upper_age = 365.25*5,
#'  type      = "qx",
#'  sex       = "total")
#'
#'lagrange5q0(data = input, model = 'B', region = "Middle Africa")



lagrange5q0 <- function(data,k, model, region){

  if(F == ("format" %in% names(data)))                     stop('Inputs were no formated with "format_data".')
  data <- data$input


  if(missing(k) == F){
    if(!is.numeric(k))                                     stop('k must be numeric.')
  }
  if(missing(model) == F){
    if(!is.character(model))                               stop('"model" must be a character string.')
  }
  if(missing(region) == F){
    if(!is.character(region))                              stop('"region" must be a character string.')
  }



  if(missing(model) == F){
    if(F %in% (model   %in%   c("A","B")))                 stop('Wrong value for "model". Use "A" or "B".')
  }else{
    model <- "A"
  }

  if(model == 'A' & missing(region) == F){
    region <- 'Model A'
                                                           warning('"region" ignored with (default) model "A". Use "region" only with Model "B".')
  }
  if(model == 'A' & missing(region) == T){
    region <- 'Model A'
  }

  data(k_bounds)
  if(model == 'A'){
  bound_max <- k_bounds1$max[k_bounds1$region == 'Model A' & k_bounds1$sex == unique(data$sex)]
  bound_min <- k_bounds1$min[k_bounds1$region == 'Model A' & k_bounds1$sex == unique(data$sex)]
  q5 <- 0.150
  }else{
  bound_max <- k_bounds1$max[k_bounds1$region == 'Model B' & k_bounds1$sex == unique(data$sex)]
  bound_min <- k_bounds1$min[k_bounds1$region == 'Model B' & k_bounds1$sex == unique(data$sex)]
  q5 <- 0.330
  }


  if(model %in%   c("B")){
    data(coef_B)
    pred <- subset(coef, sex == unique(data$sex))
  }else{
  data(coef)
    pred <- subset(coef, sex == unique(data$sex))
  }


  data(k_bounds)


  if(nrow(data) > 1){
  if(sum(c(0,7,14,21) %in% data$lower_age) == 0)           warning('No mortality input starting at age 0. Select a single mortality input instead for stable results.')
  }


  k_r <- NULL

  if(nrow(data) > 1){
    if(missing(k) == F)                                    stop('k cannot be constrained with more than one input.')
    if(model == 'B' & missing(region) == F){
      region <- 'Model B'
                                                           warning('"region" ignored with more than one input.')
    }
    if(model == 'B'){region <- 'Model B'}
    par <- newton(data,pred,0)
  }else{
    if(missing(k) == F){

      if(k < bound_min | k > bound_max)                    warning(paste('Input value of k extrapolated. k <', round(bound_min,1), 'or k >', round(bound_max,1)))

      if(model == 'B' & missing(region) == F){
                                                           warning('"region" ignored when providing a value of k.')
      }
      if(model == 'B'){region <- 'Model B'}

      par <- newton(data,pred,k)

    }else{
      if(model == 'B' & missing(region) == F){
        if(F %in% (region   %in%   c("Eastern Africa",
                                     "Middle Africa",
                                     "Western Africa",
                                     "South Asia")))       stop('Wrong value for "region". Choose "Eastern Africa", "Middle Africa",  "Western Africa", or "South Asia".')
      }

      if(model == 'B' & missing(region) == T){
        region <- 'Model B'
                                                           warning('"region" not specified. Averaged value of model "B" used (k = 0).')
      }

      k_r   <- k_bounds2$p50[k_bounds2$region   == region & k_bounds2$sex == unique(data$sex)]
      par <- newton(data,pred,k_r)
    }
  }




  if(par[1] == "error")                                    stop("Model cannot find sensical solution.")

  if(par$k < bound_min | par$k > bound_max)                warning(paste('Predicted value of k extrapolated. k <', round(bound_min,1), 'or k >', round(bound_max,1)))
  if(exp(par$h) > q5)                                      warning(paste('Predicted value of q(5y) extrapolated. q(5y) >', q5))

  pred$p_qx <- logquad(pred,par$h, par$k)
  if(is.unsorted(pred$p_qx))                               stop('Model cannot converge to a solution.')

  pred$p_mx <- qx_to_mx(pred$p_qx)
  if(model == 'A'){
  if(is.unsorted(rev(pred$p_mx)))                          warning('Increase with age in the force of mortality (nMx). Prediction extrapolated.')
  }


  bound_max_ci <- k_bounds2$p97.5[k_bounds2$region == region & k_bounds2$sex == unique(data$sex)]
  bound_min_ci <- k_bounds2$p2.5[k_bounds2$region  == region & k_bounds2$sex == unique(data$sex)]


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
       is.unsorted(pred$upper_p_qx))                       warning('Model cannot converge to a solution for confidence interval.')
    if(model == 'A'){
    if(is.unsorted(rev(pred$lower_p_mx)) |
       is.unsorted(rev(pred$upper_p_mx)))                  warning('Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.')
    }
    if(is.unsorted(pred$lower_p_qx) |
       is.unsorted(pred$upper_p_qx)){
      pred <- pred[,c("lower_age", "lower_age_m", "upper_age", "p_qx", "p_mx")]
      names(pred)[1] <- "lower_age_q"
    }else{
      pred <- pred[,c("lower_age", "lower_age_m", "upper_age",
                      "p_qx", "lower_p_qx", "upper_p_qx",
                      "p_mx", "lower_p_mx", "upper_p_mx")]
      names(pred)[1] <- "lower_age_q"
    }
  }


  if(length(k_r) > 0 & nrow(data) == 1){
    check <- 2

    par1<- newton(data,pred,  bound_max_ci)
    par2<- newton(data,pred,  bound_min_ci)

    if(par1[1] == "error" | par2[1] == "error")            warning('Model cannot converge to a solution for confidence interval.')
    if(par1[1] == "error" | par2[1] == "error"){
      pred <- pred[,c("lower_age", "lower_age_m", "upper_age", "p_qx", "p_mx")]
      names(pred)[1] <- "lower_age_q"

    }else{

      pred$lower_p_qx  <- logquad(pred,par1$h,  bound_max_ci)
      pred$upper_p_qx  <- logquad(pred,par2$h,  bound_min_ci)
      pred$lower_p_mx  <- qx_to_mx(pred$lower_p_qx)
      pred$upper_p_mx  <- qx_to_mx(pred$upper_p_qx)

      if(is.unsorted(pred$p.qx_lowerCI) |
         is.unsorted(pred$p.qx_upperCI))                   warning('Model cannot converge to a solution for confidence interval.')
      if(model == 'A'){
      if(is.unsorted(rev(pred$p.mx_lowerCI)) |
         is.unsorted(rev(pred$p.mx_upperCI)))              warning('Increase with age in the force of mortality (nMx) in confidence interval. Values extrapolated.')
      }
      if(is.unsorted(pred$p.qx_lowerCI) |
         is.unsorted(pred$p.qx_upperCI)){
        pred <- pred[,c("lower_age", "lower_age_m", "upper_age", "p_qx", "p_mx")]
        names(pred)[1] <- "low_age_q"
      }else{
        pred <- pred[,c("lower_age", "lower_age_m", "upper_age",
                        "p_qx", "lower_p_qx", "upper_p_qx",
                        "p_mx", "lower_p_mx", "upper_p_mx")]
        names(pred)[1] <- "lower_age_q"
      }
    }
  }


  if(check == 0){
    pred <- pred[,c("lower_age", "lower_age_m", "upper_age", "p_qx", "p_mx")]
    names(pred)[1] <- "lower_age_q"
  }


  tmp     <- pred
  tmp$n   <- (pred$upper_age - pred$lower_age_m)/365.25
  tmp$lx  <- 1- pred$p_qx
  tmp$dx  <- c(1,tmp$lx[-length(tmp$lx)]) - tmp$lx
  tmp$Lx  <- tmp$dx/pred$p_mx

  tmp$ID  <- ifelse(pred$lower_age_m >= 0 & pred$upper_age <= 365.25, 1, NA )
  a0_1 <- round((sum(tmp$Lx*tmp$ID, na.rm = T) - 1*(1-tmp$p_qx[pred$upper_age == 365.25]))/tmp$p_qx[pred$upper_age == 365.25],4)

  tmp$ID  <- ifelse(pred$lower_age_m >= 365.25 & pred$upper_age <= 365.25*5, 1, NA )
  a1_4 <- round((sum(tmp$Lx*tmp$ID, na.rm = T) - 4*(1-tmp$p_qx[pred$upper_age == 365.25*5]))/(tmp$p_qx[pred$upper_age == 365.25*5]-tmp$p_qx[pred$upper_age == 365.25]),4)

  tmp$ID  <- ifelse(pred$lower_age_m >= 0 & pred$upper_age <= 365.25*5, 1, NA )
  a0_5 <- round((sum(tmp$Lx*tmp$ID, na.rm = T) - 5*(1-tmp$p_qx[pred$upper_age == 365.25*5]))/tmp$p_qx[pred$upper_age == 365.25*5],4)

  rownames(pred) <- 1:nrow(pred)

  return(list( "q(5y)"       = exp(par$h),
               "k"           = par$k,
               "predictions" = pred,
               "1a0"         = a0_1,
               "4a1"         = a1_4,
               "5a0"         = a0_5))

}
