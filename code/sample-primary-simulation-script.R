#-----------------------------------------------------------------------------------------------------#
#   Simulation Studies: Sample Simulation Script                                                      #
#-----------------------------------------------------------------------------------------------------#

### Note:
  # The following represents a sample simulation script, which returns a single .csv file summarizing:
  #     1) the estimated power of the study design;
  #     2) the estimated conditional power, assuming that the true data-generating process is used to 
  #        construct the projected survival functions; and
  #     3) the estimated conditional power, assuming that the components of the projected survival 
  #        functions are estimated from/informed by the interim data
  # across 8 simulated interim datasets with the following time-to-event and study design parameters:
  #     a) low baseline incidence (\lambda_0 = 0.001),
  #     b) a mild clustering effect (\sigma^2 = 0.06),
  #     c) a null intervention effect (\beta = 0), and
  #     d) an unmatched CRT design.
  # To obtain all simulation replicates for this single simulation setting, this script must be run 
  # 125 times. Wherever possible, I have marked where the variables in the sample simulation script may 
  # be changed in order to obtain the full range of data-generating mechanisms and CRT designs considered 
  # in the primary simulation study.


source("code/cp-functions.R")
library(foreach)
library(doParallel)


#-----------------------------------------------------------------------------------------------------
# Function to Execute a Single Simulation Replicate
#-----------------------------------------------------------------------------------------------------

### Conducts one replicate of the primary simulation study
  # M               numeric         number of clusters
  # ni              vector          minimum and maximum permissable cluster sizes
  # hr              numeric         hazard ratio
  # t.dist          character       failure time distribution ("exponential"/"weibull")
  # t.theta         vector          failure time distribution parameters
  # f.dist          character       frailty distribution ("lognormal"/"gamma")
  # f.theta         vector          frailty variance parameter
  # start.times     vector          (calendar) times of cluster enrollment in the study 
  # start.delta     numeric         individual study entry times ~ Unif(start.times, start.times + start.delta)
  # visit.times     vector          (internal) time of all planned monitoring visits
  # visit.delta     numeric         individual monitoring times ~ Unif(visit.times - delta, visit.times + delta)
  # cens            numeric         rate of loss to follow-up and drop-out
  # interim.visit   numeric         monitoring visit after which the interim analysis is scheduled
  # paired          boolean         should the clusters be pair-matched?
  # l0.delta        numeric         multiplicative change in the hazard within the control arm
  # l1.delta        numeric         multiplicative change in the hazard within the intervention arm
one_run <- function(M, ni, hr, t.dist, t.theta, f.dist, f.theta, start.times, start.delta, visit.times, visit.delta,
                    cens, interim.visit, paired, l0.delta, l1.delta){
  
  # Determining the final analysis method
  if (paired){
    p_value <- matched_analysis
  } else {
    p_value <- unmatched_analysis
  }
  
  # Generating the full-trial data 
  full.out <- gen_ic_data(M, ni, hr, t.dist, t.theta, f.dist, f.theta, start.times, start.delta, visit.times, 
                          visit.delta, cens, paired)
  full.data <- format_data(full.out[[1]]); true.eta <- full.out[[2]]
  
  # Summarizing the generated full-trial data
  actual.events <- by(full.data, full.data$treat, FUN=function(x) sum(x$right != Inf))
  actual.pt <- by(full.data, full.data$treat, FUN=function(x) sum(ifelse(x$right==Inf, x$left, (x$left + x$right)/2)))
  actual.f <- est_frail(full.data, f.dist)
  actual.curves <- est_curve(full.data, f.dist, actual.f$theta)
  q0 <- quantile(ifelse(full.data$right[full.data$treat==0]==Inf, full.data$left[full.data$treat==0], full.data$right[full.data$treat==0]), probs=q)
  q1 <- quantile(ifelse(full.data$right[full.data$treat==1]==Inf, full.data$left[full.data$treat==1], full.data$right[full.data$treat==1]), probs=q)
  actual.hazard.0 <- est_hazard(actual.curves$control, end.point = q0)
  actual.hazard.1 <- est_hazard(actual.curves$treat, end.point = q1)
  actual.hr <- actual.hazard.1/actual.hazard.0
  cens.data <- full.data
  cens.data[which(full.data$right != Inf), 'left'] <- full.data$right[which(full.data$right != Inf)]
  cens.data[which(full.data$right != Inf), 'right'] <- Inf
  cens.data[which(full.data$right == Inf), 'left'] <- full.data$left[which(full.data$right == Inf)]
  cens.data[which(full.data$right == Inf), 'right'] <- sapply(full.data$left[which(full.data$right == Inf)],
                                                              FUN=function(x){ifelse(x + visit.delta > max(visit.times), Inf,
                                                                                     visit.times[min(which(visit.times > x + visit.delta))])
                                                              })
  cens.data$left[which(cens.data$left==0)] <- 0.001
  actual.cens <- exp(-survreg(Surv(left, right, type="interval2") ~ 1, data=cens.data, dist="exponential")$coefficients)
  
  
  # Creating the interim analysis dataset
  interim <- max(full.out[[1]][, grep(paste0("visit.date.", interim.visit), colnames(full.out[[1]]))], na.rm=TRUE)
  int.data <- return_int(full.out[[1]], interim)
  last.visit <- sapply(int.data$start.date, FUN=function(x) visit.times[max(which(visit.times < (interim - x - visit.delta)))])
  int.data$proj <- ifelse(last.visit==max(visit.times), 0, ifelse(int.data$right != Inf, 0, 
                                                                  ifelse(int.data$left < last.visit - visit.delta, 0, 1)))
  
  # Calculating censoring rate at interim
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
  int.cens <- exp(-survreg(Surv(left, right, type="interval2") ~ 1, data=cens.data, dist="exponential")$coefficients)
  
  # Calculating unconditional power
  est.power <- est_power(M, ni, hr, t.dist, t.theta, f.dist, f.theta, start.times, start.delta, visit.times, 
                         visit.delta, cens, paired, p_value=p_value, R=500)
  
  # Calculating conditional power: point-estimated hazard, HR, and frailties
  time.start <- Sys.time()
  est.cp <- est_cp(int.data, interim, visit.times, visit.delta, l0.delta, l1.delta, p_value, 0.05, 500, 
                   f.dist=f.dist, trunc.time=c(visit.times[interim.visit], visit.times[length(visit.times)]))
  time.end <- Sys.time()
  cp.time <- as.numeric(difftime(time.end, time.start, units="mins"))
  
  # Calculating conditional power: known hazards and frailties
  known.cp <- est_cp(int.data, interim, visit.times, visit.delta, l0.delta, l1.delta, p_value, 0.05, 500, 
                     f.dist=f.dist, eta=true.eta, f.theta=f.theta, l0=t.theta, l1=t.theta*hr, cens=cens,
                     trunc.time=c(visit.times[interim.visit], visit.times[length(visit.times)]))
  
  # Store simulation results
  res <- data.frame(M = M, t.dist = t.dist, hazard = t.theta, hr=hr, f.dist=f.dist, f.theta = f.theta, cens=cens,
                    actual.f.theta = actual.f[[2]], actual.cens = actual.cens, 
                    actual.events.0 = actual.events[1], actual.pt.0 = actual.pt[1], actual.hazard.0 = actual.hazard.0, 
                    actual.events.1 = actual.events[2], actual.pt.1 = actual.pt[2], actual.hazard.1 = actual.hazard.1, 
                    actual.hr = actual.hr, actual.p.value = p_value(full.data),
                    
                    int.f.theta = est.cp$int.theta, int.cens = int.cens, int.events.0 = est.cp$int.events[1], 
                    int.pt.0 = est.cp$int.person.time[1], int.hazard.0 = est.cp$int.hazard$control, 
                    int.events.1 = est.cp$int.events[2], int.pt.1 = est.cp$int.person.time[2], 
                    int.hazard.1 = est.cp$int.hazard$treat, int.hr=est.cp$int.hr, 
                    
                    l0.delta = l0.delta, l1.delta = l1.delta, power = est.power$power, 
                    
                    cp = est.cp$cp, time.cp=cp.time, proj.mean.f.theta = mean(est.cp$full.theta), 
                    proj.min.f.theta = min(est.cp$full.theta), proj.max.f.theta = max(est.cp$full.theta), 
                    proj.var.f.theta = var(est.cp$full.theta), proj.mse.f.theta.int = mean((est.cp$full.theta-est.cp$int.theta)^2), 
                    proj.mse.f.theta.actual = mean((est.cp$full.theta-actual.f[[2]])^2), proj.mean.cens = mean(est.cp$full.cens),
                    proj.mean.events.0 = mean(est.cp$full.events[,1]), proj.min.events.0 = min(est.cp$full.events[,1]),
                    proj.max.events.0 = max(est.cp$full.events[,1]), proj.var.events.0 = var(est.cp$full.events[,1]),
                    proj.mean.pt.0 = mean(est.cp$full.person.time[,1]), proj.min.pt.0 = min(est.cp$full.person.time[,1]),
                    proj.max.pt.0 = max(est.cp$full.person.time[,1]), proj.var.pt.0 = var(est.cp$full.person.time[,1]),
                    proj.mean.hazard.0 = mean(est.cp$full.hazard$control), proj.min.hazard.0 = min(est.cp$full.hazard$control),
                    proj.max.hazard.0 = max(est.cp$full.hazard$control), proj.var.hazard.0 = var(est.cp$full.hazard$control),
                    proj.mse.hazard.0.int = mean((est.cp$full.hazard$control - est.cp$int.hazard$control)^2), 
                    proj.mse.hazard.0.actual = mean((est.cp$full.hazard$control-actual.hazard.0)^2),
                    proj.mean.events.1 = mean(est.cp$full.events[,2]), proj.min.events.1 = min(est.cp$full.events[,2]),
                    proj.max.events.1 = max(est.cp$full.events[,2]), proj.var.events.1 = var(est.cp$full.events[,2]),
                    proj.mean.pt.1 = mean(est.cp$full.person.time[,2]), proj.min.pt.1 = min(est.cp$full.person.time[,2]),
                    proj.max.pt.1 = max(est.cp$full.person.time[,2]), proj.var.pt.1 = var(est.cp$full.person.time[,2]),
                    proj.mean.hazard.1 = mean(est.cp$full.hazard$treat), proj.min.hazard.1 = min(est.cp$full.hazard$treat),
                    proj.max.hazard.1 = max(est.cp$full.hazard$treat), proj.var.hazard.1 = var(est.cp$full.hazard$treat),
                    proj.mse.hazard.1.int = mean((est.cp$full.hazard$treat - est.cp$int.hazard$treat)^2), 
                    proj.mse.hazard.1.actual = mean((est.cp$full.hazard$treat-actual.hazard.1)^2),
                    proj.mean.hr = mean(est.cp$full.hr), proj.geom.mean.hr = exp(mean(log(est.cp$full.hr))),
                    
                    known.cp = known.cp$cp, known.mean.f.theta = mean(known.cp$full.theta), 
                    known.min.f.theta = min(known.cp$full.theta), known.max.f.theta = max(known.cp$full.theta), 
                    known.var.f.theta = var(known.cp$full.theta), known.mse.f.theta.int = mean((known.cp$full.theta-est.cp$int.theta)^2), 
                    known.mse.f.theta.actual = mean((known.cp$full.theta-actual.f[[2]])^2), known.mean.cens = mean(known.cp$full.cens),
                    known.mean.events.0 = mean(known.cp$full.events[,1]), known.min.events.0 = min(known.cp$full.events[,1]),
                    known.max.events.0 = max(known.cp$full.events[,1]), known.var.events.0 = var(known.cp$full.events[,1]),
                    known.mean.pt.0 = mean(known.cp$full.person.time[,1]), known.min.pt.0 = min(known.cp$full.person.time[,1]),
                    known.max.pt.0 = max(known.cp$full.person.time[,1]), known.var.pt.0 = var(known.cp$full.person.time[,1]),
                    known.mean.hazard.0 = mean(known.cp$full.hazard$control), known.min.hazard.0 = min(known.cp$full.hazard$control),
                    known.max.hazard.0 = max(known.cp$full.hazard$control), known.var.hazard.0 = var(known.cp$full.hazard$control),
                    known.mse.hazard.0.int = mean((known.cp$full.hazard$control - est.cp$int.hazard$control)^2), 
                    known.mse.hazard.0.actual = mean((known.cp$full.hazard$control-actual.hazard.0)^2),
                    known.mean.events.1 = mean(known.cp$full.events[,2]), known.min.events.1 = min(known.cp$full.events[,2]),
                    known.max.events.1 = max(known.cp$full.events[,2]), known.var.events.1 = var(known.cp$full.events[,2]),
                    known.mean.pt.1 = mean(known.cp$full.person.time[,2]), known.min.pt.1 = min(known.cp$full.person.time[,2]),
                    known.max.pt.1 = max(known.cp$full.person.time[,2]), known.var.pt.1 = var(known.cp$full.person.time[,2]),
                    known.mean.hazard.1 = mean(known.cp$full.hazard$treat), known.min.hazard.1 = min(known.cp$full.hazard$treat),
                    known.max.hazard.1 = max(known.cp$full.hazard$treat), known.var.hazard.1 = var(known.cp$full.hazard$treat),
                    known.mse.hazard.1.int = mean((known.cp$full.hazard$treat - est.cp$int.hazard$treat)^2), 
                    known.mse.hazard.1.actual = mean((known.cp$full.hazard$treat-actual.hazard.1)^2),
                    known.mean.hr = mean(known.cp$full.hr), known.geom.mean.hr=exp(mean(log(known.cp$full.hr))))
  
  return(res)
  
}


#-----------------------------------------------------------------------------------------------------
# (Sample) Simulation Parameter Values
#-----------------------------------------------------------------------------------------------------

# Number of clusters
M <- 30       

# Cluster size
ni <- c(250, 350)   

# Changes to the baseline hazard and hazard ratio
l0.delta <- 1
l1.delta <- 1

# Time to event and frailty distributions
t.dist <- c("exponential") 
t.theta <- 0.001            ## High Incidence Simulations: t.theta <- 0.01
f.dist <- c("lognormal")
f.theta <- 0.06             ## No Clustering: f.theta <- 0; Moderate Clustering: f.theta <- 0.22
hr <- exp(0)                ## Full Simulation Study: hr %in% exp(seq(-1, 1, by=0.1))

# Study entry and visit times
paired <- FALSE             ## Matched CRT Design: paired <- TRUE
if (paired){
  start.times <- rep(seq(0, 4*(M/2-1), by=4), each=2)
} else {
  start.times <- seq(0, 2*(M-1), by=2)
}
start.delta <- 4
visit.times <- c(52, 104, 156, 208)
visit.delta <- 4
interim.visit <- 2

# Censoring rate
cens <- 0.002

# Index
index <- 1                  ## Full Simulation Study: index %in% 1:125


#-----------------------------------------------------------------------------------------------------
# 8 Simulation Replicates, in Parallel
#-----------------------------------------------------------------------------------------------------

# Running the simulation 8 times in parallel
registerDoParallel(cores=8)
out <- foreach(n = 1:8, .combine=rbind) %dopar% one_run(M, ni, hr, t.dist, t.theta, f.dist, f.theta, start.times, start.delta, visit.times, visit.delta, cens, interim.visit, paired, l0.delta, l1.delta)
registerDoSEQ()

# Creating and writing output dataframe
out.name <- paste0(paste("est_cp", t.theta, f.theta, log(hr), index, sep="_"), ".csv")
write.csv(out, out.name, row.names=FALSE)

