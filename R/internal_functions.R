newton <- function(data, pred, k_input){

  d <- 0.0001

  if(sum(data$lower_age == 0 &  data$upper_age == 365.25*5 & data$type == "qx")>0){
    h <-log(unique(data$rate[data$lower_age == 0 &  data$upper_age == 365.25*5 & data$type == "qx"]))
  }else{
    h <- -3}

  k <- k_input

  l1 <- 1
  if(nrow(data) == 2){
    l2 <- 1}else{
    l2 <- NA
    }



  if(nrow(data) == 1){
    par  <-c("h","l1")
    par2 <-c("hh","hl1","hl1", "l1l1")
  }

  if(nrow(data) == 2){
    par  <- c("h","k","l1","l2")
    par2 <- c("hh", "hk","hl1","hl2",
              "hk","kk","kl1","kl2",
              "hl1","kl1","l1l1","l1l2",
              "hl2","kl2","l1l2","l2l2")
  }

  if(nrow(data) > 2){
    par  <-c("h","k","l1")
    par2 <-c("hh","hk","hl1",
             "hk","kk","kl1",
             "hl1", "kl1","l1l1")
  }


  repeat{

    v <- V(h,k,l1,l2,d)

    deriv1 <- NULL
    for(i in par){
      deriv1 <- c(deriv1, (L(data,pred,eval(parse(text=paste("v$",i,"[[1]]",sep = ""))))-
                           L(data,pred,eval(parse(text=paste("v$",i,"[[2]]",sep = "")))))/(2*d))
    }

    if(sum(abs(deriv1) < 1e-10) == length(deriv1)){break}
    if(sum(deriv1 == 0) == length(deriv1)){error <- 1
      break}


    deriv2 <- NULL
    for(i in par2){
      deriv2 <- c(deriv2, (L(data,pred,eval(parse(text=paste("v$",i,"[[1]]",sep = ""))))-
                           L(data,pred,eval(parse(text=paste("v$",i,"[[2]]",sep = "")))) -
                           L(data,pred,eval(parse(text=paste("v$",i,"[[3]]",sep = "")))) +
                           L(data,pred,eval(parse(text=paste("v$",i,"[[4]]",sep = "")))))/(4*d^2))
    }
    if(sum(deriv2 == 0) == length(deriv2)){error <- 1
    break}

    par_values <- NULL
    for(i in par){
      par_values <- c(par_values, eval(parse(text=i)))
    }

    M1      <- matrix(par_values,length(par),1)
    M2      <- matrix(deriv2,length(par),length(par))
    M3      <- matrix(deriv1,length(par),1)
    new_par <- as.vector(M1 - ginv(t(M2)) %*% M3)

    for(i in 1:length(par)){
      assign(paste(par[i]), new_par[i])
      }

    count <- 0
    repeat{
      count <-  count + 1
      if(count == 10){break}

      pred$p.qx <- logquad(pred,h,k)

      if(sum(pred$p.qx <= 0) > 1 |
         is.unsorted(pred$p.qx)  |
         h >= 0                  |
         k > 5 | k < -5){

        i <- as.vector(M1 - (ginv(t(M2)) %*% M3)*.3)
        h <- i[1]
        k <- i[2]
      }else{
        break}
    }

    pred$p.qx <- logquad(pred,h,k)

    if(sum(pred$p.qx < 0) > 1 |
       is.unsorted(pred$p.qx) |
       h >= 0                 |
       k > 5 | k < -5){
      error <- 1
      break
    }else{
      error <- 0}

  }

  if(error == 1){
    return("error")}else{
      return(list("h" = h, "k"=k ))}

}





V <- function(h,k,l1,l2,d){

  vect <- list(
    h    = list(c(h+d,k,l1,l2),c(h-d,k,l1,l2)),
    k    = list(c(h,k+d,l1,l2),c(h,k-d,l1,l2)),
    l1   = list(c(h,k,l1+d,l2),c(h,k,l1-d,l2)),
    l2   = list(c(h,k,l1,l2+d),c(h,k,l1,l2-d)),
    hh   = list(c(h+2*d,k,l1,l2),c(h,k,l1,l2),c(h,k,l1,l2),c(h-2*d,k,l1,l2)),
    kk   = list(c(h,k+2*d,l1,l2),c(h,k,l1,l2),c(h,k,l1,l2),c(h,k-2*d,l1,l2)),
    l1l1 = list(c(h,k,l1+2*d,l2),c(h,k,l1,l2),c(h,k,l1,l2),c(h,k,l1-2*d,l2)),
    l2l2 = list(c(h,k,l1,l2+2*d),c(h,k,l1,l2),c(h,k,l1,l2),c(h,k,l1,l2-2*d)),
    hk   = list(c(h+d,k+d,l1,l2),c(h-d,k+d,l1,l2),c(h+d,k-d,l1,l2),c(h-d,k-d,l1,l2)),
    hl1  = list(c(h+d,k,l1+d,l2),c(h-d,k,l1+d,l2),c(h+d,k,l1-d,l2),c(h-d,k,l1-d,l2)),
    hl2  = list(c(h+d,k,l1,l2+d),c(h-d,k,l1,l2+d),c(h+d,k,l1,l2-d),c(h-d,k,l1,l2-d)),
    kl1  = list(c(h,k+d,l1+d,l2),c(h,k-d,l1+d,l2),c(h,k+d,l1-d,l2),c(h,k-d,l1-d,l2)),
    kl2  = list(c(h,k+d,l1,l2+d),c(h,k-d,l1,l2+d),c(h,k+d,l1,l2-d),c(h,k-d,l1,l2-d)),
    l1l2 = list(c(h,k,l1+d,l2+d),c(h,k,l1-d,l2+d),c(h,k,l1+d,l2-d),c(h,k,l1-d,l2-d)))

  return(vect)}





L <- function(data, pred, vect){

   h <-vect[1]
   k <-vect[2]
  l1 <-vect[3]
  l2 <-vect[4]


  #logquad(pred,h,k)
  pred$p.qx <-  logquad(pred,h,k)
  pred$p.mx <-  qx_to_mx(pred$p.qx)

  if(nrow(data) > 2){
    L <- sum(f(data,pred), l1*g(data,pred)[1], l2*g(data,pred)[2], na.rm = T)
  }else{
    L <- sum(l1*g(data,pred)[1], l2*g(data,pred)[2], na.rm = T)
  }

  return(L)}





f <- function(data, pred){

  data_pred    <-  format_pred(data, pred)
  data_pred$wx <-  data_pred$weight
  rmse         <-  sum(data_pred$wx*(log(data_pred$rate) - log(data_pred$pred))^2, na.rm =T)

  return(rmse)}





g <- function(data, pred){

  tmp <- data[data$fit == "match",]

  l <- tmp$lower_age[1]
  u <- tmp$upper_age[1]
  e <- exp(1)
  pred$ID    <- ifelse(pred$lower_age_m >= l & pred$upper_age <= u, 1, NA )
  pred$ID12m <- ifelse(pred$upper_age <= 365.25, 1, NA )

  if(tmp$type[1] == "qx"){
    q    <-  1-e^-sum(pred$p.mx*(pred$n_d/365.25)*pred$ID,  na.rm = T)
    g1   <- log(tmp$rate[1]) - log(q)
  }

  if(tmp$type[1] == "mx"){
    p.lx      <- c(1, 1- pred$p.qx)
    pred$p.dx <- (p.lx[-length(p.lx)] - p.lx[-1])
    pred$p.Lx <- pred$p.dx/pred$p.mx
    m         <- sum(pred$p.dx*pred$ID, na.rm = T)/sum(pred$p.Lx*pred$ID, na.rm = T)
    g1        <- log(tmp$rate[1]) - log(m)
  }

  if(tmp$type[1] == "zx"){
    z          <- log(1-pred$p.qx[pred$upper_age == u])/log(1-pred$p.qx[pred$upper_age == 365.25])
    g1        <- log(tmp$rate[1]) - log(z)
  }





  if(nrow(tmp) == 1){
    g2 <- NA
    }

  if(nrow(tmp) == 2){
    l  <- tmp$lower_age[2]
    u  <- tmp$upper_age[2]
    e  <- exp(1)
    pred$ID    <- ifelse(pred$lower_age_m >= l & pred$upper_age <= u, 1, NA )
    pred$ID12m <- ifelse(pred$upper_age <= 365.25, 1, NA )


    if(tmp$type[2] == "qx"){
      q  <-  1-e^-sum(pred$p.mx*(pred$n_d/365.25)*pred$ID,  na.rm = T)
      g2   <- log(tmp$rate[2]) - log(q)
    }

    if(tmp$type[2] == "mx"){
      p.lx <- c(1, 1- pred$p.qx)
      pred$p.dx <- (p.lx[-length(p.lx)] - p.lx[-1])
      pred$p.Lx <- pred$p.dx/pred$p.mx
      m  <-  sum(pred$p.dx*pred$ID, na.rm = T)/sum(pred$p.Lx*pred$ID, na.rm = T)
      g2   <- log(tmp$rate[2]) - log(m)

    }

    if(tmp$type[2] == "zx"){
      z          <- log(1-pred$p.qx[pred$upper_age == u])/log(1-pred$p.qx[pred$upper_age == 365.25])
      g2         <- log(tmp$rate[2]) - log(z)
    }



  }

  rdif <- c(g1,g2)

  return(rdif)}





format_pred  <- function(data, pred){

  tmp <- data[data$fit == "min",]

  if(sum(tmp$type == "qx") > 0){

    tmp2        <- pred[, c("lower_age", "upper_age", "p.qx", "p.mx")]
    data_pred  <- merge(tmp,tmp2, c("lower_age", "upper_age"), all.x = T)
    data_pred  <- data_pred[order(data_pred$upper_age),]

    if(sum(is.na(data_pred$p.qx)) > 0){
      tmp  <- data_pred[which(is.na(data_pred$p.qx)),]
      e    <- exp(1)

      tmp2 <- lapply(1:nrow(tmp), FUN = function(i){
        pred$ID <- ifelse(pred$lower_age_m >= tmp$lower_age[i] &
                          pred$upper_age  <= tmp$upper_age[i], 1, NA)
        q       <-  1-e^-sum(pred$p.mx*(pred$n_d/365.25)*pred$ID,  na.rm = T)
        q
      })

      data_pred[is.na(data_pred$p.qx),]$p.qx <-  unlist(do.call(rbind.data.frame, tmp2))
    }

    data_pred$p.mx <- NULL
    names(data_pred)[which(names(data_pred) == "p.qx")] <- "pred"
    data_pred_qx   <- data_pred

  }else{
    data_pred_qx <- NULL}


  if(sum(tmp$type == "mx")>0){

    tmp2           <- pred[,c("lower_age_m", "upper_age", "p.qx", "p.mx")]
    names(pred)[1] <- "lower_age"
    data_pred     <- merge(tmp,tmp2, c("lower_age", "upper_age"), all.x = T)
    data_pred     <- data_pred[order(data_pred$upper_age),]


    if(sum(is.na(data_pred$p.mx)) > 0){
      tmp <- data_pred[which(is.na(data_pred$p.mx)),]

      tmp2 <- lapply(1:nrow(tmp), FUN = function(i){
        pred$ID   <- ifelse(pred$lower_age_m >= tmp$lower_age[i] &
                          pred$upper_age <= tmp$upper_age[i], 1, NA)

        p.lx      <- c(1, 1- pred$p.qx)
        pred$p.dx <- (p.lx[-length(p.lx)] - p.lx[-1])
        pred$p.Lx <- pred$p.dx/pred$p.mx
        m         <-  sum(pred$p.dx*pred$ID, na.rm = T)/sum(pred$p.Lx*pred$ID, na.rm = T)
        m
      })

      data_pred[is.na(data_pred$p.mx),]$p.mx <-  unlist(do.call(rbind.data.frame, tmp2))
    }

    data_pred$p.qx <-NULL
    names(data_pred)[which(names(data_pred) == "p.mx")] <- "pred"
    data_pred_mx   <- data_pred

  }else{
    data_pred_mx <- NULL}


  data_pred <- rbind(data_pred_qx, data_pred_mx)

  return(data_pred)
}





logquad <- function(pred,h, k){
  p.qx <- exp(pred$ax + h*pred$bx + h^2*pred$cx + pred$vx*k)
  return(p.qx)
}



qx_to_mx <- function(p.qx){
  tmp <-data.frame(
    "p.qx" = p.qx,
    "n_d"  = c(rep(7,4), 30.4375*2-28, rep(30.4375,10),  rep(30.4375*3,4), rep(365.25,3)))
  p.mx <- -log(1-c(tmp$p.qx[1],1-(1- tmp$p.qx[-1])/(1- tmp$p.qx[-22])))/(tmp$n_d/365.25)
  return(p.mx)
}


