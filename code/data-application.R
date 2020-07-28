#-----------------------------------------------------------------------------------------------------#
#   Applied Example: Simulated Botswana Combination Prevention Project (BCPP) Data                    #
#-----------------------------------------------------------------------------------------------------#

source("code/cp-functions.R")


#-----------------------------------------------------------------------------------------------------
# Data and Conditional Power Summary Functions
#-----------------------------------------------------------------------------------------------------

### Summarizes the number of events, contibuted person-time, and estimated hazard within each trial arm
  # data           dataframe    interval-censored interim data
  # trunc.time     numeric      (optional) user-specified maximum knot point for the hazard function approximation
summarize_events <- function(data, trunc.time = NULL){
  n.events.i <- sum(data$right != Inf & data$treat==1) 
  n.events.c <- sum(data$right != Inf & data$treat==0)
  f <- est_frail(data, "lognormal")
  curves <- est_curve(data, "lognormal", f$theta)
  l.t <- est_hazard(curves$treat, end.point=trunc.time)*52
  l.c <- est_hazard(curves$control, end.point=trunc.time)*52
  hr <- l.t/l.c
  out <- matrix(c(n.events.i, n.events.c, round(l.t, 5), round(l.c, 5), round(hr, 3), f[[2]]), nrow=1)
  colnames(out) <- c("events.in.treated", "events.in.control", " hazard.t", " hazard.c", "HR", "f.theta")
  print(out)
}

### Summarizes the conditional power and data generation results for a particular multiplicative change in the 
### future conditional hazard functions
  # cp.out         list         est_cp() object
  # l0.delta       numeric      multiplicative change in the hazard within the control arm
  # l1.delta       numeric      multiplicative change in the hazard within the intervention arm
  # digits         numeric      number of decimal places to be printed in the summary table
  # mean           character    type of mean to report ("arithmetic"/"geometric")
summarize_results <- function(cp.out, l0.delta, l1.delta, digits=3, mean="arithmetic"){
  out <- matrix(NA, nrow=2, ncol=10)
  out[, 1] <- c(l0.delta, ".")
  out[, 2] <- c(round(cp.out$int.hazard$control*52*l0.delta, digits), ".")
  out[, 3] <- c(l1.delta/l0.delta, ".")
  out[, 4] <- c(round((cp.out$int.hazard$treat*l1.delta)/(cp.out$int.hazard$control*l0.delta), digits), ".")
  out[, 5] <- c(round(mean(cp.out$full.events$treat), 1), paste0("(", min(cp.out$full.events$treat), ", ", max(cp.out$full.events$treat), ")"))
  out[, 6] <- c(round(mean(cp.out$full.events$control), 1), paste0("(", min(cp.out$full.events$control), ", ", max(cp.out$full.events$control), ")"))
  out[, 7] <- c(switch(mean, "arithmetic"=round(mean(cp.out$full.hr), digits), "geometric"=round(exp(mean(log(cp.out$full.hr)))), digits), round(sd(cp.out$full.hr), digits))
  out[, 8] <- c(round(mean(cp.out$full.person.time$treat), 1), round(sd(cp.out$full.person.time$treat), 1))
  out[, 9] <- c(round(mean(cp.out$full.person.time$control), 1), round(sd(cp.out$full.person.time$control), 1))
  out[, 10] <- c(cp.out$cp, ".")
  colnames(out) <- c("l.delta", "l.2", "hr.delta", "HR.2", "events.in.t", "events.in.c", "full.HR", "time.in.t", "time.in.c", "cp")
  out <- as.data.frame(out)
  print(noquote(out), row.names=FALSE)
}


#-----------------------------------------------------------------------------------------------------
# BCPP Conditional Power Analysis
#-----------------------------------------------------------------------------------------------------

# Reading in the interim data
int.data <- read.csv("data/simulated-bcpp.csv")
int.data$proj <- ifelse(int.data$right == Inf & int.data$left > 78-4, 1, 0)
int.data$start.date <- rep(0, nrow(int.data))
summarize_events(int.data, 78)



#-----------------------------------------------------------------
# Analysis: \Delta_lambda = 1; \Delta_HR = 0.8
#-----------------------------------------------------------------

set.seed(101)
start.1.0.8 <- Sys.time()
cp.1.0.8 <- est_cp(int.data, 78, c(52, 104, 156), 4, 1, 0.8, matched_analysis, 0.05, 500, cens=0.002, trunc.time=78)
end.1.0.8 <- Sys.time()

end.1.0.8 - start.1.0.8
summarize_results(cp.1.0.8, 1, 0.8)


#-----------------------------------------------------------------
# Analysis: \Delta_lambda = 1; \Delta_HR = 0.9
#-----------------------------------------------------------------

set.seed(102)
start.1.0.9 <- Sys.time()
cp.1.0.9 <- est_cp(int.data, 78, c(52, 104, 156), 4, 1, 0.9, matched_analysis, 0.05, 500, cens=0.002, trunc.time=78)
end.1.0.9 <- Sys.time()

end.1.0.9 - start.1.0.9
summarize_results(cp.1.0.9, 1, 0.9)


#-----------------------------------------------------------------
# Analysis: \Delta_lambda = 1; \Delta_HR = 1
#-----------------------------------------------------------------

set.seed(103)
start.1.1 <- Sys.time()
cp.1.1 <- est_cp(int.data, 78, c(52, 104, 156), 4, 1, 1, matched_analysis, 0.05, 500, cens=0.002, trunc.time=78)
end.1.1 <- Sys.time()

end.1.1-start.1.1
summarize_results(cp.1.1, 1, 1)


#-----------------------------------------------------------------
# Analysis: \Delta_lambda = 1; \Delta_HR = 1.1
#-----------------------------------------------------------------

set.seed(104)
start.1.1.1 <- Sys.time()
cp.1.1.1 <- est_cp(int.data, 78, c(52, 104, 156), 4, 1, 1.1, matched_analysis, 0.05, 500, cens=0.002, trunc.time=78)
end.1.1.1 <- Sys.time()

end.1.1.1-start.1.1.1
summarize_results(cp.1.1.1, 1, 1.1)


#-----------------------------------------------------------------
# Analysis: \Delta_lambda = 1; \Delta_HR = 1.2
#-----------------------------------------------------------------

set.seed(105)
start.1.1.2 <- Sys.time()
cp.1.1.2 <- est_cp(int.data, 78, c(52, 104, 156), 4, 1, 1.2, matched_analysis, 0.05, 500, cens=0.002, trunc.time=78)
end.1.1.2 <- Sys.time()

end.1.1.2 - start.1.1.2
summarize_results(cp.1.1.2, 1, 1.2)


#-----------------------------------------------------------------
# Analysis: \Delta_lambda = 0.9; \Delta_HR = 0.8
#-----------------------------------------------------------------

set.seed(106)
start.0.9.0.8 <- Sys.time()
cp.0.9.0.8 <- est_cp(int.data, 78, c(52, 104, 156), 4, 0.9, 0.9*0.8, matched_analysis, 0.05, 500, cens=0.002, trunc.time=78)
end.0.9.0.8 <- Sys.time()

end.0.9.0.8 - start.0.9.0.8
summarize_results(cp.0.9.0.8, 0.9, 0.9*0.8)


#-----------------------------------------------------------------
# Analysis: \Delta_lambda = 0.9; \Delta_HR = 0.9
#-----------------------------------------------------------------

set.seed(107)
start.0.9.0.9 <- Sys.time()
cp.0.9.0.9 <- est_cp(int.data, 78, c(52, 104, 156), 4, 0.9, 0.9*0.9, matched_analysis, 0.05, 500, cens=0.002, trunc.time=78)
end.0.9.0.9 <- Sys.time()

end.0.9.0.9 - start.0.9.0.9
summarize_results(cp.0.9.0.9, 0.9, 0.9*0.9)


#-----------------------------------------------------------------
# Analysis: \Delta_lambda = 0.9; \Delta_HR = 1
#-----------------------------------------------------------------

set.seed(108)
start.0.9.1 <- Sys.time()
cp.0.9.1 <- est_cp(int.data, 78, c(52, 104, 156), 4, 0.9, 0.9*1, matched_analysis, 0.05, 500, cens=0.002, trunc.time=78)
end.0.9.1 <- Sys.time()

end.0.9.1 - start.0.9.1
summarize_results(cp.0.9.1, 0.9, 0.9*1)


#-----------------------------------------------------------------
# Analysis: \Delta_lambda = 0.9; \Delta_HR = 1.1
#-----------------------------------------------------------------

set.seed(109)
start.0.9.1.1 <- Sys.time()
cp.0.9.1.1 <- est_cp(int.data, 78, c(52, 104, 156), 4, 0.9, 0.9*1.1, matched_analysis, 0.05, 500, cens=0.002, trunc.time=78)
end.0.9.1.1 <- Sys.time()

end.0.9.1.1 - start.0.9.1.1
summarize_results(cp.0.9.1.1, 0.9, 0.9*1.1)


#-----------------------------------------------------------------
# Analysis: \Delta_lambda = 0.9; \Delta_HR = 1.2
#-----------------------------------------------------------------

set.seed(110)
start.0.9.1.2 <- Sys.time()
cp.0.9.1.2 <- est_cp(int.data, 78, c(52, 104, 156), 4, 0.9, 0.9*1.2, matched_analysis, 0.05, 500, cens=0.002, trunc.time=78)
end.0.9.1.2 <- Sys.time()

end.0.9.1.2 - start.0.9.1.2
summarize_results(cp.0.9.1.2, 0.9, 0.9*1.2)

