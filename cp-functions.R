#-----------------------------------------------------------------------------------------------------#
#   Estimation of Conditional Power for Cluster-Randomized Trials with Interval-Censored Endpoints    #
#-----------------------------------------------------------------------------------------------------#

library(icenReg)
library(rootSolve)
library(coxme)
library(survival)
library(frailtypack)


#-----------------------------------------------------------------------------------------------------
# Data Generation Functions
#-----------------------------------------------------------------------------------------------------

### Reparametrizes the Weibull distribution
  # n             numeric         number of observations
  # shape         numeric         shape parameter for the Weibull distribution
  # scale         numeric         scale parameter for the Weibull distribution
  # beta          vector          log hazard ratio
  # x             vector          covariate data
  # eta           numeric         frailty term
rweibull_reparam <- function(n, shape, scale, beta, x, eta){
  
  u <- runif(n, 0, 1)
  s <- scale*(-log(u)/(exp(t(beta) %*% x)*eta))^(1/shape)
  
  return(as.vector(s))
}

### Generates clustered and interval-censored failure times
  # M             numeric         number of clusters
  # ni            vector          minimum and maximum permissable cluster sizes
  # hr            numeric         hazard ratio
  # t.dist        character       failure time distribution ("exponential"/"weibull")
  # t.theta       vector          failure time distribution parameters
  # f.dist        character       frailty distribution ("lognormal"/"gamma")
  # f.theta       vector          frailty variance parameter
  # start.times   vector          (calendar) times of cluster enrollment in the study 
  # start.delta   numeric         individual study entry times ~ Unif(start.times, start.times + start.delta)
  # visit.times   vector          (internal) time of all planned monitoring visits
  # visit.delta   numeric         individual monitoring times ~ Unif(visit.times - visit.delta, visit.times + visit.delta)
  # cens          numeric         rate of loss to follow-up and drop-out
  # paired        boolean         should the clusters be pair-matched?
gen_ic_data <- function(M, ni, hr, t.dist, t.theta, f.dist, f.theta, start.times, start.delta, 
                        visit.times, visit.delta, cens, paired){
  
  # Checking distributional inputs
  if (t.dist != "exponential" & t.dist != "weibull"){
    stop("Please choose either the 'exponential' or 'weibull' survival distribution.")
  }
  if ((t.dist=="exponential" & length(t.theta)!= 1) | (t.dist=="weibull" & length(t.theta)!= 2)){
    stop("Please enter the correct number of parameters for the survival distribution.")
  }
  if (f.dist != "lognormal" & f.dist != "gamma"){
    stop("Please choose either the 'lognormal' or 'gamma' distribution for the frailty.")
  }
  
  # Generating treatment vector for all observations
  clus.sizes <- sample(ni[1]:ni[2], M, replace=TRUE)
  data <- data.frame(group = rep(seq(1:M), clus.sizes))
  if (!paired){
    data$pair.id <- (data$group + data$group %% 2)/2
    data$treat <- ifelse(data$group %% 2 == 0, 1, 0)
    switch(f.dist, "gamma"={eta <- rgamma(M, shape=f.theta, rate=f.theta)},
           "lognormal"={eta <- rlnorm(M, meanlog=0, sdlog=sqrt(f.theta))})
  } else {
    switch(f.dist, "gamma"={eta <- rgamma(M, shape=f.theta, rate=f.theta)},
           "lognormal"={eta <- rlnorm(M, meanlog=0, sdlog=sqrt(f.theta))})
    if (t.dist=="exponential"){
      l.vec <- t.theta[1]*eta
    } else if (t.dist=="weibull"){
      l.vec <- (t.theta[1]/t.theta[2])*eta
    }
    pair.assignment <- order(l.vec)
    data$pair.id <- rep(NA, nrow(data))
    data$treat <- rep(NA, nrow(data))
    for (i in 1:(M/2)){
      members <- pair.assignment[c(i*2-1, i*2)]
      tx.assign <- sample(c(0, 1), 2)
      for (j in 1:2){
        data$pair.id[which(data$group==members[j])] <- i
        data$treat[which(data$group==members[j])] <- tx.assign[j]
      }
    }
  }
  
  # Generating (calendar) times of study entry 
  data$start.date <- do.call(c, apply(cbind(clus.sizes, start.times), 1, FUN=function(x)runif(x[1], x[2], x[2]+start.delta)))
  
  # Generating (internal) inspection times
  N <- nrow(data)
  visits <- matrix(NA, nrow=N, ncol=length(visit.times))
  for (i in 1:length(visit.times)){
    visits[,i] <- runif(N, max(0, visit.times[i]-visit.delta), visit.times[i]+visit.delta)
  }
  
  # Generating survival outcomes
  event <- matrix(NA, nrow=nrow(visits), ncol=ncol(visits))
  beta <- log(hr)
  for (i in 1:M){
    ids <- which(data[, 1] == unique(data[, 1])[i])
    ni <- length(ids)
    # Drawing (internal) survival and censoring times
    switch(t.dist, "exponential"={s.time <- rexp(ni, t.theta[1]*exp(t(beta) %*% data[ids, -c(1, 2, ncol(data))])*eta[i])},
           "weibull"={s.time <- rweibull_reparam(ni, shape=t.theta[1], scale=t.theta[2], beta=beta, x=data[ids, -c(1, 2, ncol(data))], eta=eta[i])})
    if (cens==0){
      c.time <- visits[ids, ncol(visits)] + 1
    } else{
      c.time <- rexp(ni, cens)
    }
    # Determining which visits were attended
    visits[ids, ] <- t(apply(cbind(visits[ids, ], c.time), 1, FUN=function(x){
      obs <- x[which(x[-length(x)] < x[length(x)])]
      out <- c(obs, rep(NA, length(x) - length(obs)-1))
      return(out)
    }))
    # Recording event indicator at each attended visit
    event[ids, ] <- t(apply(cbind(visits[ids, ], s.time), 1, FUN=function(x) return(as.numeric(x[-length(x)] > x[length(x)]))))
  }
  
  # Appending visit results to the complete-trial dataset
  visits <- visits + data$start.date
  data <- cbind(cbind(data, visits), event)
  visit.c.names <- sapply(seq(1:ncol(visits)), FUN=function(x)paste0("visit.date.", x))
  event.c.names <- sapply(seq(1:ncol(visits)), FUN=function(x)paste0("event.result.", x))
  names(data)[ncol(data) - ((2*ncol(visits)-1):0)] <- c(visit.c.names, event.c.names)
  
  # Retain pair designation?
  if (!paired){
    data$pair.id <- NULL
  }
  
  return(list(data, eta))
}  

### Transforms the wide-form (multinomial) records into interval-censored internal (study) failure times
  # data          dataframe       trial data with multinomial representation of interval censoring 
format_data <- function(data){
  
  # Working complete trial, visit, and event indicator datasets
  out.data <- data[, -c(grep("visit.date", names(data)), 
                        grep("event.result", names(data)))]
  visits <- data[, grep("visit.date", names(data))]
  event <- data[, grep("event.result", names(data))]
  
  # Deriving the observed intervals
  v.col <- seq(1:ncol(visits)); e.col <- ncol(visits) + seq(1:ncol(event))
  intervals <- t(apply(cbind(visits, event, data$start.date), 1, FUN=function(z){
    x <- z[v.col]; y <- z[e.col]; e <- z[length(z)]
    if (sum(is.na(x))==length(x)){
      return(c(0, Inf))
    }
    event.status <- max(y, na.rm=TRUE)
    if (event.status==1){
      first.event.time <- min(which(y==1))
      left <- ifelse(first.event.time==1, 0, 
                     ifelse(sum(!is.na(x[1:(first.event.time-1)]))==0, 0,
                            max(x[1:(first.event.time-1)], na.rm=TRUE)-e))
      right <- x[first.event.time]-e
    } else{
      left <- max(x, na.rm=TRUE)-e
      right <- Inf
    }
    return(c(left, right))
  }))
  
  # Creating the cleaned interim dataset
  colnames(intervals) <- c("left", "right")
  out.data <- cbind(out.data, intervals)
  
  return(out.data)
}

### Returns the interim analysis dataset corresponding to a given completed trial
  # data          dataframe       trial data, with multinomial representation of interval censoring 
  # interim       numeric         (calendar) time of the interim analysis
return_int <- function(data, interim){
  
  # Working interim, visit, and event indicator datasets
  int.data <- data[, -c(grep("visit.date", names(data)),
                        grep("event.result", names(data)))]
  visits <- data[, grep("visit.date", names(data))]
  event <- data[, grep("event.result", names(data))]
  
  # Deriving the observed intervals at interim
  v.col <- seq(1:ncol(visits)); e.col <- ncol(visits) + seq(1:ncol(event))
  interim.intervals <- t(apply(cbind(visits, event, data$start.date), 1, FUN=function(z){
    x <- z[v.col]; y <- z[e.col]; e <- z[length(z)]
    if (length(which(x <= interim))==0){
      return(c(0, Inf))
    }
    event.status <- max(y[which(x <= interim)], na.rm=TRUE)
    if (event.status==1){
      first.event.time <- min(which(y==1))
      left <- ifelse(first.event.time==1, 0, 
                     ifelse(sum(!is.na(x[1:(first.event.time-1)]))==0, 0,
                            max(x[1:(first.event.time-1)], na.rm=TRUE)-e))
      right <- x[first.event.time]-e
    } else{
      left <- max(x[which(x <= interim)], na.rm=TRUE)-e
      right <- Inf
    }
    return(c(left, right))
  }))
  
  # Creating the cleaned interim dataset
  colnames(interim.intervals) <- c("left", "right")
  int.data <- cbind(int.data, interim.intervals)
  
  return(int.data)
}


#-----------------------------------------------------------------------------------------------------
# Conditional Power Estimation Functions
#-----------------------------------------------------------------------------------------------------

### Estimates frailty terms and variance component based on the observed interim data
  # data          dataframe       interval-censored interim data
  # f.dist        character       frailty distribution ("lognormal"/"gamma")
est_frail <- function(data, f.dist){
  
  # Midpoint imputation
  rc <- data.frame(group=data$group, treat=data$treat)
  rc$delta <- ifelse(data$right==Inf, 0, 1)
  rc$time <- ifelse(data$right==Inf, data$left, (data$left + data$right)/2)
  
  if (f.dist=="lognormal"){
    # Fitting a Cox model with cluster-specific frailties
    rc.cox <- coxme(Surv(time, delta, type="right")~treat+(1 | group), data=rc)
    eta <-exp(ranef(rc.cox)$group)
    f.theta <- VarCorr(rc.cox)$group
  } else if (f.dist=="gamma"){
    # Fitting a gamma frailty model
    rc.cox <- frailtyPenal(Surv(time, delta, type="right")~treat + cluster(group), data=rc,
                           hazard="Piecewise-per", RandDist="Gamma", nb.int=10)
    eta <- rc.cox$frailty.pred
    f.theta <- 1/rc.cox$theta
  }
  
  return(list('eta'=eta, 'theta'=f.theta))
}

### Recovers the cluster-conditional survival function from the marginal survival function at time t
  # s.bar         numeric         marginal survival at time t
  # f.dist        character       frailty distribution ("lognormal"/"gamma")
  # f.theta       numeric         frailty variance parameter
recover_cond_curve <- function(s.bar, f.dist, f.theta){
  
  if (f.dist=="lognormal"){
    wrapper.surv <- function(s.bar, f.theta){
      return(fun <- function(x){
        x*(1-f.theta/2*log(x)*(-log(x)-1)) - s.bar
      })
    }
    y <- uniroot.all(wrapper.surv(s.bar, f.theta), interval=c(0, 1))
    return(ifelse(length(y)==0, 0.001, y))
    
  } else if (f.dist=="gamma"){
    return(exp(-f.theta*(s.bar^(-1/f.theta)-1)))
  }
}

### Estimates the survival function in each arm based on the observed interim data
  # data          dataframe       interval-censored interim data
  # f.dist        character       frailty distribution ("lognormal"/"gamma")
  # f.theta       numeric         frailty variance parameter
  # marginal      boolean         should the marginal or cluster-conditional survival curve be returned?
est_curve <- function(data, f.dist, f.theta, marginal=FALSE){
  
  # Turning treatment into a factor variable
  data$treat <- as.factor(data$treat)
  
  # Fitting the Turnbull curves
  turn.fit <- ic_np(Surv(left, right, type="interval2")~treat, data=data)
  cntrl.fit <- matrix(NA, nrow=nrow(turn.fit$scurves$'0'[[1]])+1, ncol=3)
  cntrl.fit[, 1] <- c(0, turn.fit$scurves$'0'[[1]][, 1])
  cntrl.fit[, 2] <- c(0, turn.fit$scurves$'0'[[1]][, 2])
  cntrl.fit[, 3] <- c(1, turn.fit$scurves$'0'[[2]]$baseline)
  tx.fit <- matrix(NA, nrow=nrow(turn.fit$scurves$'1'[[1]])+1, ncol=3)
  tx.fit[, 1] <- c(0, turn.fit$scurves$'1'[[1]][, 1])
  tx.fit[, 2] <- c(0, turn.fit$scurves$'1'[[1]][, 2])
  tx.fit[, 3] <- c(1, turn.fit$scurves$'1'[[2]]$baseline)
  
  # Ensuring the curves are estimated with sufficient information/precision
  q0 <- min(0.95, 1-25/nrow(data[data$treat==0,])); q1 <- min(0.95, 1-25/nrow(data[data$treat==1,]))
  cntrl.max <- quantile(data$left[data$treat==0], prob=q0); tx.max <- quantile(data$left[data$treat==1], prob=q1)
  cntrl.fit <- cntrl.fit[cntrl.fit[,2] < cntrl.max, ]
  tx.fit <- tx.fit[tx.fit[,2] < tx.max, ]
  
  # Ensuring that the curves span the positive real line
  if (cntrl.fit[nrow(cntrl.fit), 2] != Inf){
    cntrl.fit <- rbind(cntrl.fit, c(cntrl.fit[nrow(cntrl.fit), 2], Inf, cntrl.fit[nrow(cntrl.fit), 3]))
  }
  if (tx.fit[nrow(tx.fit), 2] != Inf){
    tx.fit <- rbind(tx.fit, c(tx.fit[nrow(tx.fit), 2], Inf, tx.fit[nrow(tx.fit), 3]))
  }
  
  if (!marginal){
    # Corrected control curve
    cntrl.surv <- sapply(cntrl.fit[, 3], FUN=function(x)recover_cond_curve(x, f.dist, f.theta))
    if (class(cntrl.surv)=="matrix"){
      warning("The non-marginal control survival function is not unique: we have selected the survival estimates that minimize the size of the correction needed.")
      cntrl.surv.2 <- apply(rbind(cntrl.surv, cntrl.fit[, 3]), 2, FUN=function(x)x[ifelse(abs(x[1]-x[3]) > abs(x[2]-x[3]), 2, 1)])
      cntrl.surv <- cntrl.surv.2
    }
    else if (class(cntrl.surv)=="list"){
      warning("The non-marginal control survival function is not unique: we have selected the survival estimates that minimize the size of the correction needed.")
      cntrl.surv.2 <- mapply(function(r, s){
        d <- sapply(r, FUN=function(x)abs(x-s))
        return(r[d==min(d)])
      }, r=cntrl.surv, s=cntrl.fit[, 3])
      cntrl.surv <- cntrl.surv.2
    }
    cntrl.fit[, 3] <- cntrl.surv
    
    # Corrected treatment curve
    tx.surv <- sapply(tx.fit[, 3], FUN=function(x)recover_cond_curve(x, f.dist, f.theta))
    if (class(tx.surv)=="matrix"){
      warning("The non-marginal intervention survival function is not unique: we have selected the survival estimates that minimize the size of the correction needed.")
      tx.surv.2 <- apply(rbind(tx.surv, tx.fit[, 3]), 2, FUN=function(x)x[ifelse(abs(x[1]-x[3]) > abs(x[2]-x[3]), 2, 1)])
      tx.surv <- tx.surv.2
    }
    else if (class(tx.surv)=="list"){
      warning("The non-marginal intervention survival function is not unique: we have selected the survival estimates that minimize the size of the correction needed.")
      tx.surv.2 <- mapply(function(r, s){
        d <- sapply(r, FUN=function(x)abs(x-s))
        return(r[d==min(d)])
      }, r=tx.surv, s=tx.fit$surv)
      tx.surv <- tx.surv.2
    }
    tx.fit[, 3] <- tx.surv
  }
  
  return(list('control'=cntrl.fit, 'treat'=tx.fit))
}

### Piecewise approximates (and averages) the conditional hazard function based on the observed interim data
  # s.curve       dataframe       estimated conditional survival function
  # end.point     numeric         maximum knot point for the piecewise approximation
est_hazard <- function(s.curve, end.point = NULL){
  
  get_surv <- function(t){
    v.lower <- s.curve[length(which(s.curve[, 1] <= t)), ]
    v.upper <- s.curve[which(s.curve[, 2] > t)[1], ]
    if (v.lower[3]==v.upper[3]){
      s <- v.lower[3]
    } else {
      s <- ((v.upper[3]-v.lower[3])/(v.upper[1]-v.lower[2]))*t + (v.lower[3]*v.upper[1] - v.upper[3]*v.lower[2])/(v.upper[1]-v.lower[2])
    }
  }
  
  # Determining the estimated survival at each knot point
  if (is.null(end.point)){
    times <- seq(0, as.integer(s.curve[nrow(s.curve)-1, 2]))
  } else {
    times <- seq(0, as.integer(end.point))
  }
  surv <- vapply(times, get_surv, numeric(1))
  surv.grid <- matrix(c(times, surv), ncol=2)
  
  # Approximating the associated conditional hazard function 
  haz.list <- rep(NA, nrow(surv.grid)-1)
  time.diff <- rep(NA, nrow(surv.grid)-1)
  for (i in 2:nrow(surv.grid)){
    haz.list[i-1] <- -log(surv.grid[i, 2]/surv.grid[i-1, 2])/(surv.grid[i, 1]-surv.grid[i-1, 1])
    time.diff[i-1] <- surv.grid[i, 1]-surv.grid[i-1, 1]
  }
  
  # Returning the mean of these piecewise components
  return(weighted.mean(haz.list, time.diff))
}

### Simulates the complete-trial record for one study participant
  # left          numeric         (internal) time of last study visit
  # interim       numeric         (internal) time of the interim analysis
  # visit.times   vector          (internal) times of planned future monitoring visits
  # visit.delta   numeric         individual monitoring times ~ Unif(visit.times - visit.delta, visit.times + visit.delta)
  # s.curve       dataframe       fitted survival curve within community i's trial arm
  # hazard        numeric         projected conditional hazard over the remainder of the study
  # eta           numeric         estimated frailty term for community i
  # cens          numeric         rate of loss to follow-up and drop-out
proj_observation <- function(left, interim, visit.times, visit.delta, s.curve, hazard, eta, cens){
  
  # Community-specific survival curve
  s.curve[, 3] <- s.curve[, 3]^eta
  
  get_surv <- function(t){
    t.inf <- min(left, s.curve[nrow(s.curve), 1])
    if (t <= t.inf){
      v.lower <- s.curve[length(which(s.curve[, 1] <= t)), ]
      v.upper <- s.curve[which(s.curve[, 2] > t)[1], ]
      if (v.lower[3]==v.upper[3]){
        s <- v.lower[3]
      } else {
        s <- ((v.upper[3]-v.lower[3])/(v.upper[1]-v.lower[2]))*t + (v.lower[3]*v.upper[1] - v.upper[3]*v.lower[2])/(v.upper[1]-v.lower[2])
      }
    } else {
      inf.lower <- s.curve[length(which(s.curve[, 1] <= t.inf)), ]
      inf.upper <- s.curve[which(s.curve[, 2] > t.inf)[1], ]
      if (inf.lower[3]==inf.upper[3]){
        s.inf <- inf.lower[3]
      } else {
        s.inf <- ((inf.upper[3]-inf.lower[3])/(inf.upper[1]-inf.lower[2]))*t.inf + (inf.lower[3]*inf.upper[1] - inf.upper[3]*inf.lower[2])/(inf.upper[1]-inf.lower[2])
      }
      s <- s.inf*exp(-hazard*eta*(t-t.inf))
    }
    return(s)
  }
  
  # Constructing the inverse survival function
  get_inv_surv <- function(s){
    t.inf <- min(left, s.curve[nrow(s.curve), 1])
    s.inf <- get_surv(t.inf)
    if (s >= s.inf){
      v.lower <- s.curve[length(which(s.curve[, 3] > s)), ]
      v.upper <- s.curve[which(s.curve[, 3] <= s)[1], ]
      t <- v.upper[1] + (s - v.upper[3])*(v.lower[2] - v.upper[1])/(v.lower[3]-v.upper[3])
    } else {
      t <- -log(s/s.inf)/(hazard*eta) + t.inf
    }
    return(max(t, left))
  }
  
  # Generating visit times
  visits <- runif(length(visit.times), vapply(visit.times - visit.delta, FUN=function(x)max(x, interim), numeric(1)), visit.times + visit.delta)
  
  # Generating survival and censoring times
  u.t <- runif(1, 0, get_surv(left))
  s.time <- get_inv_surv(u.t)
  if (cens==0){
    c.time <- visits[length(visits)] + 1
  } else{
    u.c <- runif(1, 0, exp(-cens*left))
    c.time <- -log(u.c)/cens
  }
  
  # Deriving final interval censored observation and return (left, right)
  visits[visits > c.time] <- NA
  event <- as.numeric(visits > s.time)
  if (sum(is.na(event))==length(event)){
    return(c(left, Inf))
  }
  event.status <- max(event, na.rm=TRUE)
  if (event.status==1){
    first.event.time <- min(which(event==1))
    l <- ifelse(first.event.time==1, left, 
                ifelse(sum(!is.na(visits[1:(first.event.time-1)]))==0, left,
                       max(visits[1:(first.event.time-1)], na.rm=TRUE)))
    r <- visits[first.event.time]
  } else{
    l <- max(visits, na.rm=TRUE)
    r <- Inf
  }
  
  return(c(l, r))
}

### Simulates one complete-trial dataset corresponding to a given interim dataset
  # int.data      dataframe       interval-censored interim data
  # interim       numeric         (calendar) time of the interim analysis
  # visit.times   numeric         (internal) time of all planned monitoring visits
  # visit.delta   numeric         individual monitoring times ~ Unif(visit.times - visit.delta, visit.times + visit.delta)
  # s.curves      list            estimated survival curves in each trial arm
  # hazards       list            projected conditional hazard function in each trial arm
  # eta           numeric         estimated community-specific frailty terms
  # cens          numeric         projected rate of loss to follow-up and drop-out
proj_full_data <- function(int.data, interim, visit.times, visit.delta, s.curves, hazards, eta, cens){
  
  # Identifying individuals still in the risk set at interim
  risk.set <- which(int.data$proj==1)
  
  # Identifying column positions of relevant components
  tx.col <- which(colnames(int.data)=="treat")
  e.col <- which(colnames(int.data)=="start.date")
  i.col <- which(colnames(int.data)=="group")
  l.col <- which(colnames(int.data)=="left")
  
  # Obtaining projected end-of-trial observations for those in the risk set
  projected.intervals <- t(apply(int.data[risk.set, ], 1, FUN=function(z){
    tx <- z[tx.col]; s.curve <- s.curves[[tx+1]]; hazard <- hazards[[tx+1]]
    i <- z[i.col]; eta.i <- eta[i]
    e <- z[e.col]; int.ij <- interim - e
    l.ij <- z[l.col]
    visits.ij <- visit.times[((visit.times + visit.delta) > int.ij) & ((visit.times - visit.delta) > l.ij)]
    return(proj_observation(l.ij, int.ij, visits.ij, visit.delta, s.curve, hazard, eta.i, cens))
  }))
  
  # Constructing the projected dataset
  proj.data <- int.data
  proj.data[risk.set, c('left', 'right')] <- projected.intervals
  
  return(proj.data)
}

### Estimates study conditional power based on available interim data
  # int.data       dataframe     interval-censored interim data
  # interim        numeric       (calendar) time of the interim analysis
  # visit.times    numeric       (internal) time of all planned monitoring visits
  # visit.delta    numeric       individual monitoring times ~ Unif(visit.times - visit.delta, visit.times + visit.delta)
  # l0.delta       numeric       multiplicative change in the hazard within the control arm
  # l1.delta       numeric       multiplicative change in the hazard within the intervention arm
  # p_value        function      planned final analysis and p-value calculation
  # alpha          numeric       prespecified significance level for the final analysis
  # R              numeric       number of projected datasets with which to estimate the conditional power
  # f.dist         character     frailty distribution ("lognormal"/"gamma")
  # eta            numeric       (optional) vector of user-specified community-specific frailty terms
  # f.theta        numeric       (optional) user-specified frailty variance parameter
  # l0             numeric       (optional) user-specified conditional hazard within the control arm at interim
  # l1             numeric       (optional) user-specified conditional hazard within the intervention arm at interim
  # cens           numeric       (optional) user-specified rate of loss to follow-up and drop-out
  # q              numeric       (optional) user-specified quantile of the maximum knot point for the hazard function approximation
  # trunc.time     numeric       (optional) user-specified maximum knot point for the hazard function approximation
est_cp <- function(int.data, interim, visit.times, visit.delta, l0.delta, l1.delta, p_value, alpha, R, 
                   f.dist="lognormal", eta=NULL, f.theta=NULL, l0=NULL, l1=NULL, cens=NULL, q=NULL, trunc.time=NULL){
  
  call <- match.call()
  
  # Initializing vectors
  p.values <- rep(NA, R)
  full.l0 <- rep(NA, R); full.l1 <- rep(NA, R); full.hr <- rep(NA, R); full.cens <- rep(NA, R)
  full.eta <- matrix(NA, nrow=R, ncol=length(unique(int.data$group))); full.f.theta <- rep(NA, R)
  full.events <- data.frame(control=rep(NA, R), treat=rep(NA, R))
  full.p.time <- data.frame(control=rep(NA, R), treat=rep(NA, R))
  
  # Calculating quantities at interim
  if (is.null(eta)){
    f <- est_frail(int.data, f.dist)
  } else{
    f <- list('eta'=eta, 'theta'=f.theta)
  }
  s.curves <- est_curve(int.data, f.dist, f$theta)
  if (is.null(q) & is.null(trunc.time)){
    q0 <- q1 <- visit.times[sum(visit.times < (min(interim-int.data$start.date) + visit.delta))]
  } else if (is.null(q)){
    q0 <- q1 <- trunc.time[1]
  } else {
    q0 <- quantile(ifelse(int.data$right[int.data$treat==0]==Inf, int.data$left[int.data$treat==0], int.data$right[int.data$treat==0]), probs=q)
    q1 <- quantile(ifelse(int.data$right[int.data$treat==1]==Inf, int.data$left[int.data$treat==1], int.data$right[int.data$treat==1]), probs=q)
  }
  if (is.null(l0)){
    l0 <- est_hazard(s.curves$control, end.point = q0)
  }
  if (is.null(l1)){
    l1 <- est_hazard(s.curves$treat, end.point = q1)
  }
  if (is.null(cens)){
    cens.data <- int.data
    cens.data[which(int.data$right != Inf), 'left'] <- int.data$right[which(int.data$right != Inf)]
    cens.data[which(int.data$right != Inf), 'right'] <- Inf
    cens.data[which(int.data$right == Inf & int.data$proj == 0), 'left'] <- int.data$left[which(int.data$right == Inf & int.data$proj == 0)]
    cens.data[which(int.data$right == Inf & int.data$proj == 0), 'right'] <- sapply(int.data$left[which(int.data$right == Inf & int.data$proj == 0)], 
                                                                                    FUN=function(x){
                                                                                      ifelse(x + visit.delta > max(visit.times), Inf, 
                                                                                             visit.times[min(which(visit.times > x + visit.delta))])
                                                                                    })
    cens.data$left[which(cens.data$left==0)] <- 0.001
    cens <- exp(-survreg(Surv(left, right, type="interval2") ~ 1, data=cens.data, dist="exponential")$coefficients)
  }
  hr <- l1/l0
  
  # Calculating projected hazards
  hazards <- list('control'=l0*l0.delta, 'treat'=l1*l1.delta)
  
  for (i in 1:R){
    # Generate end-of-trial data under desired assumptions
    full.data <- proj_full_data(int.data, interim, visit.times, visit.delta, s.curves, hazards, f[[1]], cens)
    
    # Storing attributes of full data
    full.f <- est_frail(full.data, f.dist)
    full.eta[i, ] <- full.f$eta
    full.f.theta[i] <- full.f$theta
    full.curves <- est_curve(full.data, f.dist, full.f$theta)
    if ((is.null(q) & is.null(trunc.time)) | (is.null(q) & length(trunc.time)==1)){
      q0 <- q1 <- visit.times[length(visit.times)]
    } else if (is.null(q)){
      q0 <- q1 <- trunc.time[2]
    } else {
      q0 <- quantile(ifelse(full.data$right[full.data$treat==0]==Inf, full.data$left[full.data$treat==0], full.data$right[full.data$treat==0]), probs=q)
      q1 <- quantile(ifelse(full.data$right[full.data$treat==1]==Inf, full.data$left[full.data$treat==1], full.data$right[full.data$treat==1]), probs=q)
    }
    full.l0[i] <- est_hazard(full.curves$control, end.point = q0)
    full.l1[i] <- est_hazard(full.curves$treat, end.point = q1)
    full.hr[i] <- full.l1[i]/full.l0[i]
    cens.data <- full.data
    cens.data[which(full.data$right != Inf), 'left'] <- full.data$right[which(full.data$right != Inf)]
    cens.data[which(full.data$right != Inf), 'right'] <- Inf
    cens.data[which(full.data$right == Inf), 'left'] <- full.data$left[which(full.data$right == Inf)]
    cens.data[which(full.data$right == Inf), 'right'] <- sapply(full.data$left[which(full.data$right == Inf)],
                                                                FUN=function(x){ifelse(x + visit.delta > max(visit.times), Inf,
                                                                                       visit.times[min(which(visit.times > x + visit.delta))])
                                                                })
    cens.data$left[which(cens.data$left==0)] <- 0.001
    full.cens[i] <- exp(-survreg(Surv(left, right, type="interval2") ~ 1, data=cens.data, dist="exponential")$coefficients)
    full.events[i, ] <- by(full.data, full.data$treat, FUN=function(x) sum(x$right != Inf))
    full.p.time[i, ] <- by(full.data, full.data$treat, FUN=function(x) sum(ifelse(x$right==Inf, x$left, (x$left + x$right)/2)))
    
    # Calculating the corresponding p-value
    p.values[i] <- p_value(full.data)
  }
  
  return(list('call'=call, 'interim'=interim, 'cp'=mean(p.values <= alpha), 'p.values'=p.values, 
              'int.hazard'=list('control'=l0, 'treat'=l1), 'int.hr'=hr, 
              'full.hazard'=list('control'=full.l0, 'treat'=full.l1), 'full.hr'=full.hr, 
              'f.dist'=f.dist, 'int.eta'=f$eta, 'int.theta'=f$theta, 'full.eta'=full.eta, 
              'full.theta'=full.f.theta, 'int.cens'=cens, 'full.cens'=full.cens, 
              'int.events'=by(int.data, int.data$treat, FUN=function(x)sum(x$right!=Inf)),
              'int.person.time'=by(int.data, int.data$treat, FUN=function(x) sum(ifelse(x$right==Inf, x$left, (x$left + x$right)/2))), 
              'full.events'=full.events, 'full.person.time'=full.p.time))
}


#-----------------------------------------------------------------------------------------------------
# Power Estimation Function
#-----------------------------------------------------------------------------------------------------

### Estimates study power
  # M              numeric      number of clusters
  # ni             vector       two-dimensional vector of minimum and maximum cluster sizes
  # hr             numeric      hazard ratio
  # t.dist         character    failure time distribution ("exponential"/"weibull")
  # t.theta        vector       parameters of the failure time distribution (rate for exponential; shape and scale for weibull)
  # f.dist         character    frailty distribution ("lognormal"/"gamma")
  # f.theta        vector       frailty variance parameter
  # start.times    vector       calendar times corresponding to each cluster's enrollment into the study 
  # start.delta    numeric      individual study entry dates ~ Unif(start.times, start.times + start.delta)
  # visit.times    vector       scheduled internal times of the monitoring visits
  # visit.delta    numeric      individual monitoring visits ~ Unif(visit.times - delta, visit.times + delta)
  # cens           numeric      rate of loss to follow-up and drop-out
  # paired         boolean      should the clusters be pair-matched?
  # p_value        function     planned final analysis and p-value calculation
  # R              numeric      number of simulations from which to calculate the power
  # features       logical      should a summary of the simulated data features be returned?
  # q              numeric      (optional) user-specified quantile of the maximum knot point for the hazard function approximation
  # trunc.time     numeric      (optional) user-specified maximum knot point for the hazard function approximation
est_power <- function(M, ni, hr, t.dist, t.theta, f.dist, f.theta, start.times, start.delta, visit.times, 
                      visit.delta, cens, paired, p_value, R, features=FALSE, q=NULL){
  
  call <- match.call()
  
  # Check inputs for covariate generation
  if (M %% 2 != 0){
    stop("Must specify even number of clusters")
  }
  if (ni[1] > ni[2]){
    stop("Please input cluster size specifications as c(lower, upper)")
  }
  
  if (features){
    full.l0 <- rep(NA, R); full.l1 <- rep(NA, R); full.hr <- rep(NA, R); full.cens <- rep(NA, R)
    full.eta <- matrix(NA, nrow=R, ncol=length(unique(int.data$group))); full.f.theta <- rep(NA, R)
    full.events <- data.frame(control=rep(NA, R), treat=rep(NA, R))
    full.p.time <- data.frame(control=rep(NA, R), treat=rep(NA, R))
  }
  p.values <- rep(NA, R)
  
  for (i in 1:R){
    # Generating the complete data
    full.data <- gen_ic_data(M, ni, hr, t.dist, t.theta, f.dist, f.theta, start.times, start.delta,
                             visit.times, visit.delta, cens, paired)[[1]]
    full.data <- format_data(full.data)
    
    if (features){
      full.f <- est_frail(full.data, f.dist)
      full.eta[i, ] <- full.f$eta
      full.f.theta[i] <- full.f$theta
      full.curves <- est_curve(full.data, f.dist, full.f$theta)
      if (is.null(q) & is.null(trunc.time)){
        q0 <- q1 <- NULL
      } else if (is.null(q)){
        q0 <- q1 <- trunc.time[1]
      } else {
        q0 <- quantile(ifelse(full.data$right[full.data$treat==0]==Inf, full.data$left[full.data$treat==0], full.data$right[full.data$treat==0]), probs=q)
        q1 <- quantile(ifelse(full.data$right[full.data$treat==1]==Inf, full.data$left[full.data$treat==1], full.data$right[full.data$treat==1]), probs=q)
      }
      full.l0[i] <- est_hazard(full.curves$control, end.point = q0)
      full.l1[i] <- est_hazard(full.curves$treat, end.point = q1)
      full.hr[i] <- full.l1[i]/full.l0[i]
      cens.data <- full.data
      cens.data[which(full.data$right != Inf), 'left'] <- full.data$right[which(full.data$right != Inf)]
      cens.data[which(full.data$right != Inf), 'right'] <- Inf
      cens.data[which(full.data$right == Inf), 'left'] <- full.data$left[which(full.data$right == Inf)]
      cens.data[which(full.data$right == Inf), 'right'] <- sapply(full.data$left[which(full.data$right == Inf)],
                                                                  FUN=function(x){ifelse(x + visit.delta > max(visit.times), Inf,
                                                                                         visit.times[min(which(visit.times > x + visit.delta))])
                                                                  })
      cens.data$left[which(cens.data$left==0)] <- 0.001
      full.cens[i] <- exp(-survreg(Surv(left, right, type="interval2") ~ 1, data=cens.data, dist="exponential")$coefficients)
      full.events[i, ] <- by(full.data, full.data$treat, FUN=function(x) sum(x$right != Inf))
      full.p.time[i, ] <- by(full.data, full.data$treat, FUN=function(x) sum(ifelse(x$right==Inf, x$left, (x$left + x$right)/2)))
    }
    
    # Calculating the corresponding p-value
    p.values[i] <- p_value(full.data)
  }
  
  out <- list('call'=call, 'power'=mean(p.values <= 0.05, 'p_values'=p.values))
  if (features){
    out <- append(out, list('full.hazard'=list('control'=full.l0, 'treat'=full.l1), 'full.hr'=full.hr, 
                            'f.dist'=f.dist, 'full.eta'=full.eta, 'full.theta'=full.f.theta, 'full.cens'=full.cens, 
                            'full.events'=full.events, 'full.person.time'=full.p.time))
  }
  
  return(out)
}


#-----------------------------------------------------------------------------------------------------
# Final Analysis Functions 
#-----------------------------------------------------------------------------------------------------

### Performs permutation test from Wang and De Gruttola (2017) under a matched CRT design
  # data           dataframe    trial data for final analysis
matched_analysis <- function(data){
  
  # Calculating the within-pair incidence difference
  inc.diff <- c()
  for (pair in unique(data$pair.id)){
    pair.data <- data[which(data$pair.id==pair), ]
    incidence <- by(pair.data, pair.data$treat, FUN=function(x) sum(x$right !=Inf)/
                      sum(ifelse(x$right==Inf, x$left, (x$left + x$right)/2)))
    inc.diff <- c(inc.diff, incidence[2]-incidence[1])
  }
  
  # Calculating the test statistic of interest
  test.stat <- sum(inc.diff)
  
  # Determining the associated p-value
  perm.stat <- rep(NA, 1000)
  for (r in 1:1000){
    coefs <- sample(c(-1, 1), length(inc.diff), replace=TRUE, prob=c(0.5, 0.5))
    perm.stat[r] <- sum(inc.diff*coefs)
  }
  
  return((length(which(abs(perm.stat) >= abs(test.stat))) + 1)/(1000+1))
}

### Performs permutation test from Wang and De Gruttola (2017) under an unmatched CRT design
  # data           dataframe    trial data for final analysis
unmatched_analysis <- function(data){
  
  M <- length(unique(data$group))
  n.tx <- length(unique(data$group[data$treat==1]))
  n.cntrl <- length(unique(data$group[data$treat==0]))
  
  # Calculating incidence within each community
  inc <- as.numeric(by(data, data$group, FUN=function(x)sum(x$right != Inf)/sum(ifelse(x$right==Inf, x$left, (x$left + x$right)/2))))
  treat <- as.numeric(by(data, data$group, FUN=function(x)unique(x$treat)))
  
  # Finding the test statistic within the observed data
  test.mat <- matrix(c(rep(-inc[treat==0], each=n.tx), rep(inc[treat==1], n.cntrl)), ncol=2)
  test.stat <- sum(rowSums(test.mat))
  
  # Finding the permutation distribution
  perm.stat <- rep(NA, 1000)
  for (r in 1:1000){
    id.tx <- sample(1:M, n.tx, replace=FALSE)
    id.cntrl <- setdiff(1:M, id.tx)
    perm.mat <- matrix(c(rep(-inc[id.cntrl], each=n.tx), rep(inc[id.tx], n.cntrl)), ncol=2)
    perm.stat[r] <- sum(rowSums(perm.mat))
  }
  
  return((length(which(abs(perm.stat) >= abs(test.stat))) + 1)/(1000+1))
}
