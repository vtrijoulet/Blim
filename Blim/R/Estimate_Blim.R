##################################################################################################################
####          Simple function to estimate Blim according to the different stock types in ICES (2021):         ####
####          https://doi.org/10.17895/ices.advice.7891                                                       ####
####          Author: Vanessa Trijoulet, vtri@aqua.dtu.dk                                                     ####
##################################################################################################################



# Estimate hockey-stick SRR 
HS_fit <- function(y, x, years, doplot, log=FALSE, ...){
  
  fn <- function(par){
    i <- (x < exp(par[1]))
    slope <- exp(par[2]) / exp(par[1])
    sum( (log(y[i]) - log(slope * x[i]))^ 2 ) + sum( ( log(y[!i]) - par[2] )^ 2  )
  }
  n = 100
  gr <- expand.grid(x=seq(log(min(x)),log(max(x)),length=n), y=seq(log(min(y)),log(max(y)),length=n))
  m <- apply(gr,1,fn)
  # Get close to global optimum
  par = unlist(gr[which.min(m),]) # par are x,y coordinates for inflection point
  # Find the exact global optimum
  opt <- optim(par, fn)
  
  ssb <- seq(0.1, max(x)*1.5, length=1000)
  pred_rec <- ssb
  log.pred.y <- y
  slope <- exp(opt$par[2]) / exp(opt$par[1])
  i <- (ssb < exp(opt$par[1]))
  k <- (x < exp(opt$par[1]))
  if (length(pred_rec[i])>0) pred_rec[i] <- log(slope * ssb[i])
  if (length(pred_rec[!i])>0) pred_rec[!i] <-(opt$par[2])
  if (length(log.pred.y[k])>0) log.pred.y[k] <- log(slope * x[k])
  if (length(log.pred.y[!k])>0) log.pred.y[!k] <-(opt$par[2])
  
  
  sd.log.y <- sd(log(y)-log.pred.y)
  
  if (doplot){
    CI.low <- pred_rec-1.96*sd.log.y
    CI.up <- pred_rec+1.96*sd.log.y
    if (log){
      plot(y=log(y), x=log(x), type="l", col="red", xlab="log(SSB)", ylab="log(Recruitment)", ...)
      text(log(x), log(y), labels = names(x), cex = 0.7, col = "red")
      #plot(y=log(y), x=log(x), pch=16, type="p", xlab="log(SSB)", ylab="log(Recruitment)", ...)
      lines(pred_rec~log(ssb))
      polygon(x=c(sort(log(ssb)),rev(sort(log(ssb)))), y=c(sort(CI.low),rev(sort(CI.up))), col=rgb(0,0,0,alpha=0.2), border = NA)
    } else {
      plot(y=y, x=x, type="l", col="red", ylim=c(0, max(y)*1.1), xlim=c(0, max(x)*1.1), xlab="SSB", ylab="Recruitment", ...)
      text(x, y, labels = names(x), cex = 0.7, col = "red")
      #plot(y=y, x=x, pch=16, type="o", col="red", ylim=c(0, max(y)*1.1), xlim=c(0, max(x)*1.1), xlab="SSB", ylab="Recruitment", ...)
      lines(exp(pred_rec)~ssb)
      polygon(x=c(sort(ssb),rev(sort(ssb))), y=c(sort(exp(CI.low)),rev(sort(exp(CI.up)))), col=rgb(0,0,0,alpha=0.2), border = NA)
    }
    
  }
  tmp <- c(as.vector(opt$par), sd.log.y, opt$convergence)
  names(tmp) <- c("logBlim", "logRmax", "sd_logR", "convergence_code")
  
  return(tmp)
}


# Simple linear regression without constraint on intercept
lm_fn <- function(y, x, years, doplot, ...){
  opt <- lm(y~x)
  #opt <- lm(log(y)~log(x))
  if (doplot){
    plot(y=y, x=x, type="l", col="red", ylim=c(0, max(y)*1.1), xlim=c(0, max(x)*1.1), xlab="SSB", ylab="Recruitment", ...)
    text(x, y, labels = names(x), cex = 0.7, col = "red")
    #plot(y=y, x=x, pch=16, type="o", col="red", ylim=c(0, max(y)*1.1), xlim=c(0, max(x)*1.1), xlab="SSB", ylab="Recruitment", ...)
    sample <- seq(0, max(x)*1.1, length=1000)
    tmp=predict(opt,se.fit=TRUE,newdata=list(x=sample))
    CI.low <- tmp$fit-1.96*tmp$se.fit
    CI.up <- tmp$fit+1.96*tmp$se.fit
    polygon(x=c(sample,rev(sample)), y=c(CI.low,rev(CI.up)), col=rgb(0,0,0,alpha=0.2), border = NA)
    abline(opt)
  }
  return(opt)
}


# Function to determine stock type when recruitment is flat line (Types 5 or 6)
flat_rec <- function(x, narrow_threshold, doplot, ...){
  threshold <- max(x)/2
  tmp <- quantile(x,probs=seq(0.05,0.95,0.05))<threshold
  if (sum(tmp==TRUE)==0) prob_max <- "0%" else prob_max <- names(tmp[which(tmp==TRUE)][length(tmp[which(tmp==TRUE)])])
  if (as.numeric(gsub("%", "", prob_max))/100 < narrow_threshold){ # Type 6, narrow range, less than buffer of SSB below max(SSB)/2
    message <- "Type 6, No Blim from this data, only the PA reference point"
    cat(message)
    cat("\n")
    res <- NA
    attr(res, "code") <- message
    if (doplot) abline(v=threshold, lty=2, col="grey")
    return(res)
  } else { # Type 5, normal range
    message <- "Type 5, No evidence of impaired recruitment or no clear SRR: Blim is Bloss"
    cat(message)
    cat("\n")
    res <- min(x)
    names(res) <- "min(SSB)"
    attr(res, "code") <- message
    if (doplot) abline(v=res, lty=2, lwd=2)
    return(res)
  }
}


##' Function to estimate Blim given a time series of recruitment and SSB
##' @param data SAM fit, FLStock object, or data frame with columns "Rec" and "SSB" and with years as row names.
##' @param Rage age at recruitment to correctly identify SSB-Rec pairs, only needed if data is a data frame.
##' @param years the years to consider for the SSB-Rec pairs, default is the whole time series.
##' @param is.spasmodic \code{TRUE} if the stock is assumed to have spasmodic recruitment, default is \code{FALSE}.
##' @param types stock types to consider for Blim estimation, default are types 2 to 6 in ICES (2021). 
##' @param doplot if TRUE, the final stock-recruitment plot will be drawn during estimation, default is \code{FALSE}.
##' @param spasmodic_prob the probability to use to identify outliers in recruitment in the quantile function when \code{is.spasmodic=TRUE}, default is 0.9.
##' @param narrow_treshold threshold to determine if the range of SSB is narrow or not (difference between stock types 5 and 6), default is 0.25. The range of SSB is considered narrow if less than \code{narrow_threshold} of SSB is below max(SSB)/2.
##' @param conf_level_lm confidence level of the confidence interval of the slope parameter of the linear regression (types 3 to 6), default is 0.95.
##' @param ... extra arguments for internal functions.
##' @details The function estimates Blim according to the different stock types in ICES (2021, \url{https://doi.org/10.17895/ices.advice.7891}).
##' If \code{is.spasmodic=TRUE}, stock type 1, Blim is estimated as the lowest SSB, where large recruitment (probability \code{>=spasmodic_prob}) is observed, default is \code{FALSE}. Stock type 1 will only be considered if the argument is \code{TRUE}. \cr
##' \cr 
##' By default, the function will estimate Blim according to the stock types 2-6.
##' First, the function tries to estimate Blim as the inflection point of a hockey-stick stock-recruitment (SRR) curve (stock type 2). 
##' If the inflection point is outside the range of observed SSB, the stock type is considered to be either 5 or 6 if the inflection point 
##' is below the minimum SSB, or type 3 if above. 
##' If a consistent SRR cannot be estimated (non convergence), Blim will be estimated according to types 3-6 as explained below. \cr
##' \cr
##' To determine which type is the best between types 3-6, an unconstrained linear regression is fitted to the SSB-Rec pairs.
##' The slope indicates which type the stock is and Blim is estimated accordingly to ICES (2021).
##' If the slope of the linear regression is not different from zero (confidence interval given \code{conf_level_lm} contains zero), 
##' the SRR is considered to be flat (types 5 or 6) and the final type is chosen depending on the range of SSB values (narrow or not). If the slope is positive, the type 3 is assumed and the type 4 is assumed otherwise.                                         
##' \cr
##' It is possible to avoid fitting a hockey-stick SRR by using \code{types!=2}. In this case, any stock type between 3 and 6 can be considered to estimate Blim. 
##' @return Blim and stock type code as attribute.
##' @author Vanessa Trijoulet
##' @examples
##' s <- seq(100, 2000, length=50)
##' r <- exp((2+log(s)-0.001*s)+rnorm(length(s)))
##' data <- data.frame(Rec=r, SSB=s, row.names=1971:2020)
##' 
##' # The stock has spasmodic recruitment (type 1):
##' Blim_estim(data, Rage=0, is.spasmodic=TRUE, doplot=TRUE) # only type 1 is assumed
##' 
##' # The algorithm tests for ICES stock types 2-6 and chooses the best:
##' Blim_estim(data, Rage=0, doplot=TRUE) 
##' 
##' # Type 2 is not considered and the algorithm only tests for types 3-6:
##' Blim_estim(data, Rage=0, types=3:6, doplot=TRUE) 
##' @export
Blim_estim <- function (data, ...){
  UseMethod("Blim_estim")
}

##' @rdname Blim_estim
##' @method Blim_estim default
##' @export
Blim_estim.default <- function(data, Rage, years=as.numeric(rownames(data)), is.spasmodic=FALSE, types=2:6, doplot=FALSE, spasmodic_prob=0.9, narrow_threshold=0.25, conf_level_lm=0.95, ...){
  
  if (sum(colnames(data) %in% c("Rec", "SSB")) != 2) stop("The column names of data should be 'Rec' and 'SSB'")

  R <- as.numeric(data[, "Rec"])
  S <- as.numeric(data[, "SSB"])
  names(R) <- as.numeric(rownames(data))
  names(S) <- as.numeric(rownames(data))
  
  if (missing(Rage)) stop("The argument Rage (recruitment age) is missing")
  
  if (sum(types==1)>0) cat("Type 1 can only be considered by setting is.spasmodic=TRUE. Only the types 3-6 will be considered", sep="\n")
  
  if (Rage>0) {
    y=R[as.character(years[-(1:Rage)])]
    x=S[as.character(years[-((length(years)-(Rage-1)):length(years))])]
  } else {
    y=R[as.character(years)]
    x=S[as.character(years)]
  }
  
  if (is.spasmodic){ # Type 1
    threshold <- quantile(y, probs=spasmodic_prob)
    outlier_idx <- which(y>=threshold)
    res <- min(x[outlier_idx]) # could also use S instead of x!!!!
    names(res) <- "Blim"
    if (doplot){
      plot(y=y, x=x, type="l", col="red", ylim=c(0, max(y)*1.1), xlim=c(0, max(x)*1.1), xlab="SSB", ylab="Recruitment", ...)
      text(x, y, labels = names(x), cex = 0.7, col = "red")
      #plot(y=y, x=x, pch=16, type="o", col="red", ylim=c(0, max(y)*1.1), xlim=c(0, max(x)*1.1), xlab="SSB", ylab="Recruitment", ...)
      abline(h=threshold, lty=2, col="grey")
      abline(v=res, lty=2, lwd=2)
    }
    message <- "Type 1, Blim is based on the lowest SSB, where large recruitment is observed"
    attr(res, "code") <- message
    cat(message)
    cat("\n")
    return (res)
  } else {
    if (2 %in% types){  # Only if we want to fit a HS SRR
      # Try fit HS SRR
      SRR <- HS_fit(y, x, years, doplot, ...)
      if (SRR["convergence_code"]==0){
        if (exp(SRR["logBlim"])>min(x)){ 
          if (exp(SRR["logBlim"])>max(x)){ # Type 3 
            message <- "Type 3, Blim may be close to the highest SSB observed"
            cat(message)
            cat("\n")
            res <- max(x) # could also use S instead of x!!!!
            names(res) <- "max(SSB)"
            attr(res, "code") <- message
            if (doplot) abline(v=res, lty=2, lwd=2)
            return(res)
          } else { # Type 2 # could also use S instead of x!!!!
            message <- "Type 2, Blim is the inflection point of the hockey-stick curve"
            cat(message)
            cat("\n")
            res <- exp(SRR["logBlim"])
            names(res) <- "Blim"
            attr(res, "code") <- message
            if(doplot) abline(v=res, lty=2, lwd=2)
            return (res)
          }
        } else { # Recruitment is flat, Type 5 or 6
          flat_rec(x, narrow_threshold, doplot, ...) # could also use S instead of x!!!!
        }
      } else {
        # message <- "It is not possible to fit a consistent hockey-stick stock recruitment relationship to the data, try the other types"
        # cat(message)
        # res <- NA
        # attr(res, "code") <- message
        # return(res)
        
        # Fit a linear regression without constraint on intercept (for Types 3-6)
        res_lm <- lm_fn(y, x, years, doplot, ...)
        CI_slope <- confint(res_lm, level = conf_level_lm)["x",]
        if (CI_slope[1]<=0 & CI_slope[2]>=0){ # Recruitment is flat, Type 5 or 6
          flat_rec(x, narrow_threshold, doplot, ...) # could also use S instead of x!!!!
        } else {
          if (CI_slope[1]>=0){ # Type 3
            message <- "Type 3, Blim may be close to the highest SSB observed"
            cat(message)
            cat("\n")
            res <- max(x) # could also use S instead of x!!!!
            names(res) <- "max(SSB)"
            attr(res, "code") <- message
            if (doplot) abline(v=res, lty=2, lwd=2)
            return(res)
          } else { # Type 4
            message <- "Type 4, No Blim from this data, only the PA reference point"
            cat(message)
            cat("\n")
            res <- NA
            attr(res, "code") <- message
            return(res)
          }
        }
      }
    } else { # Fit a linear regression without constraint on intercept (for Types 3-6)
      res_lm <- lm_fn(y, x, years, doplot, ...)
      CI_slope <- confint(res_lm, level = conf_level_lm)["x",]
      if (CI_slope[1]<=0 & CI_slope[2]>=0){ # Recruitment is flat, Type 5 or 6
        flat_rec(x, narrow_threshold, doplot, ...) # could also use S instead of x!!!!
      } else {
        if (CI_slope[1]>=0){ # Type 3
          message <- "Type 3, Blim may be close to the highest SSB observed"
          cat(message)
          cat("\n")
          res <- max(x) # could also use S instead of x!!!!
          names(res) <- "max(SSB)"
          attr(res, "code") <- message
          if (doplot) abline(v=res, lty=2, lwd=2)
          return(res)
        } else { # Type 4
          message <- "Type 4, No Blim from this data, only the PA reference point"
          cat(message)
          cat("\n")
          res <- NA
          attr(res, "code") <- message
          return(res)
        }
      }
    }
  }
}



##' @rdname Blim_estim
##' @method Blim_estim sam
##' @export
Blim_estim.sam <- function(data, years=as.numeric(rownames(rectable(data))), is.spasmodic=FALSE, types=2:6, doplot=FALSE, spasmodic_prob=0.9, narrow_threshold=0.25, conf_level_lm=0.95, ...){
  
  R <- rectable(data)[, "Estimate"]
  S <- ssbtable(data)[, "Estimate"]
  Rage <- data$conf$minAge

  if (sum(types==1)>0) cat("Type 1 can only be considered by setting is.spasmodic=TRUE. Only the types 3-6 will be considered", sep="\n")
  
  if (Rage>0) {
    y=R[as.character(years[-(1:Rage)])]
    x=S[as.character(years[-((length(years)-(Rage-1)):length(years))])]
  } else {
    y=R[as.character(years)]
    x=S[as.character(years)]
  }
  
  if (is.spasmodic){ # Type 1
    threshold <- quantile(y, probs=spasmodic_prob)
    outlier_idx <- which(y>=threshold)
    res <- min(x[outlier_idx]) # could also use S instead of x!!!!
    names(res) <- "Blim"
    if (doplot){
      plot(y=y, x=x, type="l", col="red", ylim=c(0, max(y)*1.1), xlim=c(0, max(x)*1.1), xlab="SSB", ylab="Recruitment", ...)
      text(x, y, labels = names(x), cex = 0.7, col = "red")
      #plot(y=y, x=x, pch=16, type="o", col="red", ylim=c(0, max(y)*1.1), xlim=c(0, max(x)*1.1), xlab="SSB", ylab="Recruitment", ...)
      abline(h=threshold, lty=2, col="grey")
      abline(v=res, lty=2, lwd=2)
    }
    message <- "Type 1, Blim is based on the lowest SSB, where large recruitment is observed"
    attr(res, "code") <- message
    cat(message)
    cat("\n")
    return (res)
  } else {
    if (2 %in% types){  # Only if we want to fit a HS SRR
      # Try fit HS SRR
      SRR <- HS_fit(y, x, years, doplot, ...)
      if (SRR["convergence_code"]==0){
        if (exp(SRR["logBlim"])>min(x)){ 
          if (exp(SRR["logBlim"])>max(x)){ # Type 3 
            message <- "Type 3, Blim may be close to the highest SSB observed"
            cat(message)
            cat("\n")
            res <- max(x) # could also use S instead of x!!!!
            names(res) <- "max(SSB)"
            attr(res, "code") <- message
            if (doplot) abline(v=res, lty=2, lwd=2)
            return(res)
          } else { # Type 2 # could also use S instead of x!!!!
            message <- "Type 2, Blim is the inflection point of the hockey-stick curve"
            cat(message)
            cat("\n")
            res <- exp(SRR["logBlim"])
            names(res) <- "Blim"
            attr(res, "code") <- message
            if(doplot) abline(v=res, lty=2, lwd=2)
            return (res)
          }
        } else { # Recruitment is flat, Type 5 or 6
          flat_rec(x, narrow_threshold, doplot, ...) # could also use S instead of x!!!!
        }
      } else {
        # message <- "It is not possible to fit a consistent hockey-stick stock recruitment relationship to the data, try the other types"
        # cat(message)
        # res <- NA
        # attr(res, "code") <- message
        # return(res)
        
        # Fit a linear regression without constraint on intercept (for Types 3-6)
        res_lm <- lm_fn(y, x, years, doplot, ...)
        CI_slope <- confint(res_lm, level = conf_level_lm)["x",]
        if (CI_slope[1]<=0 & CI_slope[2]>=0){ # Recruitment is flat, Type 5 or 6
          flat_rec(x, narrow_threshold, doplot, ...) # could also use S instead of x!!!!
        } else {
          if (CI_slope[1]>=0){ # Type 3
            message <- "Type 3, Blim may be close to the highest SSB observed"
            cat(message)
            cat("\n")
            res <- max(x) # could also use S instead of x!!!!
            names(res) <- "max(SSB)"
            attr(res, "code") <- message
            if (doplot) abline(v=res, lty=2, lwd=2)
            return(res)
          } else { # Type 4
            message <- "Type 4, No Blim from this data, only the PA reference point"
            cat(message)
            cat("\n")
            res <- NA
            attr(res, "code") <- message
            return(res)
          }
        }
      }
    } else { # Fit a linear regression without constraint on intercept (for Types 3-6)
      res_lm <- lm_fn(y, x, years, doplot, ...)
      CI_slope <- confint(res_lm, level = conf_level_lm)["x",]
      if (CI_slope[1]<=0 & CI_slope[2]>=0){ # Recruitment is flat, Type 5 or 6
        flat_rec(x, narrow_threshold, doplot, ...) # could also use S instead of x!!!!
      } else {
        if (CI_slope[1]>=0){ # Type 3
          message <- "Type 3, Blim may be close to the highest SSB observed"
          cat(message)
          cat("\n")
          res <- max(x) # could also use S instead of x!!!!
          names(res) <- "max(SSB)"
          attr(res, "code") <- message
          if (doplot) abline(v=res, lty=2, lwd=2)
          return(res)
        } else { # Type 4
          message <- "Type 4, No Blim from this data, only the PA reference point"
          cat(message)
          cat("\n")
          res <- NA
          attr(res, "code") <- message
          return(res)
        }
      }
    }
  }
}

##' @rdname Blim_estim
##' @method Blim_estim FLStock
##' @export
Blim_estim.FLStock <- function(data, years=as.numeric(names(rec(data)[drop=TRUE])), is.spasmodic=FALSE, types=2:6, doplot=FALSE, spasmodic_prob=0.9, narrow_threshold=0.25, conf_level_lm=0.95, ...){

  R <- rec(data)[drop=TRUE]
  S <- ssb(data)[drop=TRUE]
  Rage <- as.numeric(rownames(rec(data)))
  
  if (sum(types==1)>0) cat("Type 1 can only be considered by setting is.spasmodic=TRUE. Only the types 3-6 will be considered", sep="\n")
  
  if (Rage>0) {
    y=R[as.character(years[-(1:Rage)])]
    x=S[as.character(years[-((length(years)-(Rage-1)):length(years))])]
  } else {
    y=R[as.character(years)]
    x=S[as.character(years)]
  }
  
  if (is.spasmodic){ # Type 1
    threshold <- quantile(y, probs=spasmodic_prob)
    outlier_idx <- which(y>=threshold)
    res <- min(x[outlier_idx]) # could also use S instead of x!!!!
    names(res) <- "Blim"
    if (doplot){
      plot(y=y, x=x, type="l", col="red", ylim=c(0, max(y)*1.1), xlim=c(0, max(x)*1.1), xlab="SSB", ylab="Recruitment", ...)
      text(x, y, labels = names(x), cex = 0.7, col = "red")
      #plot(y=y, x=x, pch=16, type="o", col="red", ylim=c(0, max(y)*1.1), xlim=c(0, max(x)*1.1), xlab="SSB", ylab="Recruitment", ...)
      abline(h=threshold, lty=2, col="grey")
      abline(v=res, lty=2, lwd=2)
    }
    message <- "Type 1, Blim is based on the lowest SSB, where large recruitment is observed"
    attr(res, "code") <- message
    cat(message)
    cat("\n")
    return (res)
  } else {
    if (2 %in% types){  # Only if we want to fit a HS SRR
      # Try fit HS SRR
      SRR <- HS_fit(y, x, years, doplot, ...)
      if (SRR["convergence_code"]==0){
        if (exp(SRR["logBlim"])>min(x)){ 
          if (exp(SRR["logBlim"])>max(x)){ # Type 3 
            message <- "Type 3, Blim may be close to the highest SSB observed"
            cat(message)
            cat("\n")
            res <- max(x) # could also use S instead of x!!!!
            names(res) <- "max(SSB)"
            attr(res, "code") <- message
            if (doplot) abline(v=res, lty=2, lwd=2)
            return(res)
          } else { # Type 2 # could also use S instead of x!!!!
            message <- "Type 2, Blim is the inflection point of the hockey-stick curve"
            cat(message)
            cat("\n")
            res <- exp(SRR["logBlim"])
            names(res) <- "Blim"
            attr(res, "code") <- message
            if(doplot) abline(v=res, lty=2, lwd=2)
            return (res)
          }
        } else { # Recruitment is flat, Type 5 or 6
          flat_rec(x, narrow_threshold, doplot, ...) # could also use S instead of x!!!!
        }
      } else {
        # message <- "It is not possible to fit a consistent hockey-stick stock recruitment relationship to the data, try the other types"
        # cat(message)
        # res <- NA
        # attr(res, "code") <- message
        # return(res)
        
        # Fit a linear regression without constraint on intercept (for Types 3-6)
        res_lm <- lm_fn(y, x, years, doplot, ...)
        CI_slope <- confint(res_lm, level = conf_level_lm)["x",]
        if (CI_slope[1]<=0 & CI_slope[2]>=0){ # Recruitment is flat, Type 5 or 6
          flat_rec(x, narrow_threshold, doplot, ...) # could also use S instead of x!!!!
        } else {
          if (CI_slope[1]>=0){ # Type 3
            message <- "Type 3, Blim may be close to the highest SSB observed"
            cat(message)
            cat("\n")
            res <- max(x) # could also use S instead of x!!!!
            names(res) <- "max(SSB)"
            attr(res, "code") <- message
            if (doplot) abline(v=res, lty=2, lwd=2)
            return(res)
          } else { # Type 4
            message <- "Type 4, No Blim from this data, only the PA reference point"
            cat(message)
            cat("\n")
            res <- NA
            attr(res, "code") <- message
            return(res)
          }
        }
      }
    } else { # Fit a linear regression without constraint on intercept (for Types 3-6)
      res_lm <- lm_fn(y, x, years, doplot, ...)
      CI_slope <- confint(res_lm, level = conf_level_lm)["x",]
      if (CI_slope[1]<=0 & CI_slope[2]>=0){ # Recruitment is flat, Type 5 or 6
        flat_rec(x, narrow_threshold, doplot, ...) # could also use S instead of x!!!!
      } else {
        if (CI_slope[1]>=0){ # Type 3
          message <- "Type 3, Blim may be close to the highest SSB observed"
          cat(message)
          cat("\n")
          res <- max(x) # could also use S instead of x!!!!
          names(res) <- "max(SSB)"
          attr(res, "code") <- message
          if (doplot) abline(v=res, lty=2, lwd=2)
          return(res)
        } else { # Type 4
          message <- "Type 4, No Blim from this data, only the PA reference point"
          cat(message)
          cat("\n")
          res <- NA
          attr(res, "code") <- message
          return(res)
        }
      }
    }
  }
}

