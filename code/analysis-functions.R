#-----------------------------------------------------------------------------------------------------#
#   Simulation Studies: Analysis Functions                                                            #
#-----------------------------------------------------------------------------------------------------#

library(ggplot2)
library(grid)
library(gridExtra)
library(ggmosaic)
library(rootSolve)
library(viridis)


#-----------------------------------------------------------------------------------------------------
# General Utility Functions
#-----------------------------------------------------------------------------------------------------

### Combines plots with a shared legend
### Taken from https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
  # plots           list          list of ggplot objects to be combined
  # ncol            numeric       number of columns in the grid arrangement
  # nrow            numeric       number of rows in the grid arrangement
  # position        character     location of the shared legend in the final arrangement
grid_arrange_shared_legend <- function(plots, ncol = length(plots), nrow = 1, position = c("bottom", "right"), title) {
    
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[2]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }


#-----------------------------------------------------------------------------------------------------
# Data Generation Functions
#-----------------------------------------------------------------------------------------------------

### Prints the mean number of events, accumulated person-time, and hazard (stratified by trial arm), as well as the 
### estimated full-trial hazard ratio and log-frailty variance, in both the original and projected complete-trial datasets
  # data            dataframe     conditional power results from any of the simulation scripts 
  # digits          numeric       number of decimal places to be printed in the summary table
  # beta            numeric       log hazard ratios for which to print the data generation performance 
summarize_data_gen <- function(data, digits=3, beta=NULL){
  
  if(!is.null(beta)){
    data <- data[round(log(data$hr), 2) %in% beta, ]
  }
  
  for (l in unique(data$hazard)){
    h.data <- data[data$hazard==l, ]
    
    cat("Hazard: "); cat(l); cat("\n\n")
    print.out <- matrix(NA, nrow=2*length(unique(h.data$f.theta))*length(unique(h.data$hr)),
                        ncol=11)
    i <- 1
    for (v in sort(unique(h.data$f.theta))){
      for (h in sort(unique(h.data$hr))){
        dat <- h.data[h.data$f.theta==v & h.data$hr==h, ]
        results.actual <- apply(cbind(dat$actual.events.0, dat$actual.pt.0, dat$actual.hazard.0,
                                      dat$actual.events.1, dat$actual.pt.1, dat$actual.hazard.1,
                                      dat$actual.hazard.1/dat$actual.hazard.0, dat$actual.f.theta),
                                2, mean)
        results.proj <- apply(cbind(dat$proj.mean.events.0, dat$proj.mean.pt.0, 
                                    dat$proj.mean.hazard.0, dat$proj.mean.events.1,
                                    dat$proj.mean.pt.1, dat$proj.mean.hazard.1,
                                    dat$proj.mean.hazard.1/dat$proj.mean.hazard.0, dat$proj.mean.f.theta),
                              2, mean)
        print.out[i,1] <- ifelse(i %% (2*length(unique(h.data$hr)))==1, v, ".")
        print.out[i, 2] <- round(log(h), 2)
        print.out[i, 3] <- "Original Trial"
        print.out[i, 4:ncol(print.out)] <- unlist(lapply(results.actual, FUN=function(x)round(x, digits=digits)))
        print.out[i+1, c(1,2)] <- "."
        print.out[i+1, 3] <- "Projected Trial"
        print.out[i+1, 4:ncol(print.out)] <- unlist(lapply(results.proj, FUN=function(x)round(x, digits=digits)))
        i <- i+2
      }
    }
    colnames(print.out) <- c("S^2", "Beta", "", "Events (X=0)", "PT (X=0)", "Hazard (X=0)", 
                             "Events (X=1)", "PT (X=1)", "Hazard (X=1)", "HR", "RE Var")
    print.out <- as.data.frame(print.out)
    print(noquote(print.out), row.names=FALSE)
  }
}

### Formats data for the bias plots
  # data            dataframe     conditional power results from any of the simulation scripts
format_bias_data <- function(data){
  
  comp.dat <- data.frame()
  i <- 1
  for (l in sort(unique(data$hazard))){
    for (v in sort(unique(data$f.theta))){
      for (b in sort(unique(data$hr))){
        dat <- data[which(data$hazard==l & data$f.theta==v & data$hr==b),]
        if (nrow(dat) > 0){
          n <- nrow(dat)
          proj.bias.h0 <- dat$proj.mean.hazard.0 - dat$actual.hazard.0
          known.bias.h0 <- dat$known.mean.hazard.0 - dat$actual.hazard.0
          proj.bias.h1 <- dat$proj.mean.hazard.1 - dat$actual.hazard.1
          known.bias.h1 <- dat$known.mean.hazard.1 - dat$actual.hazard.1
          proj.bias.var <- dat$proj.mean.f.theta - dat$actual.f.theta
          known.bias.var <- dat$known.mean.f.theta - dat$actual.f.theta
          
          temp <- data.frame(est=c(proj.bias.h0, known.bias.h0, proj.bias.h1, known.bias.h1,
                                   proj.bias.var, known.bias.var),
                             param=c(rep("h0",2*n), rep("h1",2*n), rep("var",2*n)),
                             method=rep(c(rep("proj", n), rep("known", n)), 3),
                             beta=rep(log(b), 6*n), var=rep(v, 6*n), lambda=rep(l, 6*n),
                             x.pos=rep(i, 6*n))
          
          comp.dat <- rbind(comp.dat, temp)
          i <- i+1
        }
      }
    }
  }
  
  return(comp.dat)
}

### Formats data for the MSE plots
  # data            dataframe     conditional power results from any of the simulation scripts
format_mse_data <- function(data){
  
  comp.dat <- data.frame()
  i <- 1
  for (l in sort(unique(data$hazard))){
    for (v in sort(unique(data$f.theta))){
      for (b in sort(unique(data$hr))){
        dat <- data[which(data$hazard==l & data$f.theta==v & data$hr==b),]
        if (nrow(dat) > 0){
          n <- nrow(dat)
          proj.mse.h0 <- dat$proj.mse.hazard.0.actual
          known.mse.h0 <- dat$known.mse.hazard.0.actual
          proj.mse.h1 <- dat$proj.mse.hazard.1.actual
          known.mse.h1 <- dat$known.mse.hazard.1.actual
          proj.mse.var <- dat$proj.mse.f.theta.actual
          known.mse.var <- dat$known.mse.f.theta.actual
          
          temp <- data.frame(est=c(proj.mse.h0, known.mse.h0, proj.mse.h1, known.mse.h1,
                                   proj.mse.var, known.mse.var),
                             param=c(rep("h0",2*n), rep("h1",2*n), rep("var",2*n)),
                             method=rep(c(rep("proj", n), rep("known", n)), 3),
                             beta=rep(log(b), 6*n), var=rep(v, 6*n), lambda=rep(l, 6*n),
                             x.pos=rep(i, 6*n))
          
          comp.dat <- rbind(comp.dat, temp)
          i <- i+1
        }
      }
    }
  }
  
  return(comp.dat)
}

### Visualizes bias and MSE of the projected-trial quantities as estimators of their complete-trial counterparts
  # data            dataframe     conditional power results from any of the simulation scripts
  # beta            numeric       log hazard ratios for which to calculate data generation performance 
visualize_data_gen <- function(data, beta){
  
  b.dat <- format_bias_data(data[round(log(data$hr), 2) %in% beta, ])
  levels(b.dat$param) <- c(expression(bar(lambda)[0]), expression(bar(lambda)[1]), expression(hat(sigma)^2))
  b.dat$hline <- 0
  m.dat <- format_mse_data(data[round(log(data$hr), 2) %in% beta, ])
  levels(m.dat$param) <- c(expression(bar(lambda)[0]), expression(bar(lambda)[1]), expression(hat(sigma)^2))
  m.dat$hline <- 0
  
  n.breaks <- length(unique(b.dat$x.pos))
  x.labels <- character(n.breaks)
  for (i in 1:n.breaks){
    x.name <- unique(b.dat$x.pos)[i]
    x.setting <- ifelse(unique(b.dat$lambda[b.dat$x.pos==x.name])==0.001, "i", "ii")
    x.beta <- format(round(unique(b.dat$beta[b.dat$x.pos==x.name]), 1), nsmall=1)
    x.labels[i] <- ifelse(nchar(x.beta)==4, paste0(x.setting, ": ", x.beta), paste0(x.setting, ":   ", x.beta))
    names(x.labels)[i] <- x.name
  }
  
  l.b <- ggplot(data=b.dat, aes(x=as.factor(x.pos), y=est, fill=method))
  l.b <- l.b + geom_boxplot() + scale_fill_viridis(name="Projected Values   ", labels=c(expression(paste("Data-Generating (", lambda[0], ", ", lambda[1], ", ", sigma^2, ", ", eta, ")   ")), expression(paste("Estimated (", tilde(lambda)[0], ", ", tilde(lambda)[1], ", ", hat(sigma)^2, ", ", hat(eta), ")"))),begin=0.5, end=0.85, discrete=TRUE, direction=-1) + 
    scale_x_discrete(name=" ", labels=x.labels) + scale_y_continuous(name="Bias") + theme_classic() +
    theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1, size=20), axis.text.y = element_text(size=20), axis.title.x=element_text(size=28), axis.title.y=element_text(margin=margin(t = 0, r = 20, b = 0, l = 0), size=28), legend.text=element_text(size=28), strip.background=element_blank(), strip.text.x = element_text(size = 28), panel.border = element_rect(colour = "black", fill=NA, size=1))
  l.b <- l.b + facet_wrap(~ as.factor(param), ncol=3, scales="free", labeller = label_parsed) + geom_hline(data=b.dat, aes(yintercept=hline), color="red", linetype="dashed")
  
  l.m <- ggplot(data=m.dat, aes(x=as.factor(x.pos), y=est, fill=method))
  l.m <- l.m + geom_boxplot() + scale_fill_viridis(name="Projected Values   ", labels=c(expression(paste(" Data-Generating (", lambda[0], ", ", lambda[1], ", ", sigma^2, ", ", eta, ")   ")), expression(paste(" Estimated (", tilde(lambda)[0], ", ", tilde(lambda)[1], ", ", hat(sigma)^2, ", ", hat(eta), ")"))),begin=0.5, end=0.85, discrete=TRUE, direction=-1) + 
    scale_x_discrete(name=expression(paste("CRT Setting: ", beta)), labels=x.labels) + scale_y_continuous(name="Mean Squared Error") + theme_classic() +
    theme(legend.position="bottom", axis.text.x = element_text(angle = 45, hjust=1, size=20), axis.text.y = element_text(size=20), axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 20, l = 0), size=28), axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=28), legend.text=element_text(size=26), legend.title=element_text(size=26), strip.background=element_blank(), strip.text.x = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1))
  l.m <- l.m + facet_wrap(~ as.factor(param), ncol=3, scales="free", labeller = label_parsed) + geom_hline(data=m.dat, aes(yintercept=hline), color="red", linetype="dashed") 
  
  p <- grid_arrange_shared_legend(list(l.b, l.m), nrow=2, ncol=1) 
  
  p
}


#-----------------------------------------------------------------------------------------------------
# Conditional Power Functions
#-----------------------------------------------------------------------------------------------------

### Formats data for the conditional power plots (and optionally prints a table summarizing the results)
  # data            dataframe     conditional power results from any of the simulation scripts
  # print           boolean       should a summary table also be printed?
  # digits          numeric       number of decimal places to be printed in the summary table
format_plot_data <- function(data, print=TRUE, digits=3){
  out <- data.frame()
  
  # For each unique baseline hazard simulation setting...
  for (l in unique(data$hazard)){
    if (print) cat("Hazard: "); cat(l); cat("\n\n")
    temp <- data[data$hazard==l, ]
    if (print) print.out <- matrix(NA, ncol=5, nrow=length(unique(temp$hr))*length(unique(temp$f.theta))); i <- 1
    
    # And for each unique clustering effect...
    for (v in sort(unique(temp$f.theta))){
      
      # And hazard ratio...
      for (h in sort(unique(temp$hr))){
        
        # Calculate the mean (conditional) power across all simulation replicates
        dat <- temp[temp$f.theta==v & temp$hr==h, ]
        results.mean <- apply(cbind(dat$power, dat$cp, dat$known.cp), 2, mean)
        
        # And optionally format the summary table for printing
        if (print){
          results.sd <- apply(cbind(dat$power, dat$cp, dat$known.cp), 2, sd)
          print.out[i, ] <- c(ifelse(length(unique(temp$hr))==1 | i %% length(unique(temp$hr))==1, v, "."), 
                              round(log(h), digits), 
                              paste0(round(results.mean[1], digits), " (", round(results.sd[1], digits), ")"),
                              paste0(round(results.mean[2], digits), " (", round(results.sd[2], digits), ")"), 
                              paste0(round(results.mean[3], digits), " (", round(results.sd[3], digits), ")"))
          i <- i + 1
        }
        
        res.out <- data.frame(hazard=rep(l, 3), f.theta=rep(v, 3), hr=rep(h, 3), est=results.mean, method=c("Power", "Conditional Power", "Known Hazards & Frailty"))
        out <- rbind(out, res.out)
      }
    }
    
    # Printing out the summary table
    if (print){
      colnames(print.out) <- c("S^2", "Beta", "Power", "CP (Original)", "CP (All Known)")
      print.out <- as.data.frame(print.out)
      print(noquote(print.out), row.names=FALSE)
    }
  }
  
  return(out)
}

### Generates plot with power and conditional power curves
  # plot.data       dataframe     formatted data from format_plot_data()
visualize_cp_curves <- function(plot.data){
  plot.data$beta <- round(log(plot.data$hr), 2)
  plot.data$method <- droplevels(plot.data$method)
  plot.data$method <- factor(plot.data$method, levels=levels(plot.data$method)[c(3, 1, 2)])
  title.lab <- substitute(paste(lambda[0], " = ", l, ", ", sigma^2, " = ", v), list(l=unique(plot.data$hazard), v=format(round(unique(plot.data$f.theta), 2), nsmall=2)))
  
  p <- ggplot(data=plot.data, aes(x=beta, y=est, color=method))
  p <- p + geom_point() + geom_line(size=0.7) + theme_classic() + scale_color_viridis(labels=c("Power", "Conditional Power", expression(paste("True (", lambda[0], ", ", lambda[1], ", ", sigma^2, ", ", eta, ")"))), begin=0, end=0.8, discrete=TRUE) + guides(color = guide_legend(reverse = TRUE))
  p <- p + geom_hline(yintercept=0.05, linetype=2, color="red") + scale_y_continuous(breaks = c(seq(0.25,1,by = 0.25), 0.05)) + scale_x_continuous(breaks=sapply(seq(min(unique(plot.data$beta)), max(unique(plot.data$beta)), by=0.1), FUN=function(x)round(x, 2))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p <- p + labs(x=expression(paste("Data-Generating ", beta)), y="Estimated (Conditional) Power", color="") 
  p <- p + theme(legend.position = c(0.185, 0.18), legend.text.align=0, legend.key.size = unit(2,"line"), plot.title = element_text(hjust = 0.5, size=23), axis.text=element_text(size=16), axis.title.x=element_text(margin=margin(t = 15, r = 0, b = 15, l = 0), size=21), axis.title.y=element_text(margin=margin(t = 0, r = 5, b = 0, l = 0), size=21), legend.text=element_text(size=20), legend.background=element_rect(fill="transparent"), panel.border = element_rect(colour = "black", fill=NA, size=1)) + ggtitle(title.lab)
  p
}


#-----------------------------------------------------------------------------------------------------
# Classification Accuracy Functions
#-----------------------------------------------------------------------------------------------------

### Prints false positive rate, false negative rate, accuracy rate, sensitivity, and specificity of the interim futility 
### classifications
  # data            dataframe     conditional power results from any of the simulation scripts
  # alpha           numeric       significance level (final analysis)
  # threshold       numeric       futility threshold (interim analysis)
  # digits          numeric       number of decimal places to be printed in the summary table
  # beta            numeric       (optional) log hazard ratios over which to stratify the results
summarize_accuracy <- function(data, alpha=0.05, threshold=0.20, digits=3, beta=NULL){
  
  cat("Baseline Hazard: "); cat(unique(data$hazard)); cat("\n")
  cat("Coefficient of Variation: "); cat(round(sqrt(exp(unique(data$f.theta))-1), 2)); cat("\n\n")
  
  if (!is.character(beta) & !is.null(beta)){
    data <- data[round(log(data$hr), 2) %in% beta, ]
  }
  
  # Computing measures of classification performance
  n.futile <- sum(data$cp < threshold)
  n.type.II <- sum(data$cp < threshold & data$actual.p.value <= alpha)
  n.continue <- sum(data$cp >= threshold)
  n.type.I <- sum(data$cp > threshold & data$actual.p.value > alpha)
  correct <- sum((data$cp > threshold) == (data$actual.p.value < alpha))
  sens.n <- sum(data$cp > threshold & data$actual.p.value < alpha)
  sens.d <- (sum(data$cp > threshold & data$actual.p.value < alpha) + sum(data$cp < threshold & data$actual.p.value < alpha))
  sens <- sum(data$cp > threshold & data$actual.p.value < alpha)/(sum(data$cp > threshold & data$actual.p.value < alpha) + sum(data$cp < threshold & data$actual.p.value < alpha))
  spec.n <- sum(data$cp < threshold & data$actual.p.value > alpha)
  spec.d <- (sum(data$cp < threshold & data$actual.p.value > alpha) + sum(data$cp > threshold & data$actual.p.value > alpha))
  spec <- sum(data$cp < threshold & data$actual.p.value > alpha)/(sum(data$cp < threshold & data$actual.p.value > alpha) + sum(data$cp > threshold & data$actual.p.value > alpha))
  overall.results <- data.frame(n.futile = n.futile, n.type.ii.errors = n.type.II, type.ii.rate = n.type.II/n.futile,
                                n.not.futile=n.continue, n.type.i.errors = n.type.I, type.i.rate = n.type.I/n.continue,
                                n.correct.decisions=correct, accuracy=correct/nrow(data),
                                sensitivity = sens, sens.n = sens.n, sens.d = sens.d,
                                specificity = spec, spec.n = spec.n, spec.d = spec.d,
                                youdon = sens + spec - 1)
  
  # Printing classification performance results, stratifying on the given beta values
  if(!is.null(beta)){
    out.results <- data.frame()
    if (is.character(beta)){
      beta <- sort(unique(round(log(data$hr), 2)))
    }
    for (b in beta){
      dat <- data[round(log(data$hr), 2) == b, ]
      n.futile <- sum(dat$cp < threshold)
      n.type.II <- sum(dat$cp < threshold & dat$actual.p.value <= alpha)
      n.continue <- sum(dat$cp >= threshold)
      n.type.I <- sum(dat$cp > threshold & dat$actual.p.value > alpha)
      correct <- sum((dat$cp > threshold) == (dat$actual.p.value < alpha))
      sens.n <- sum(dat$cp > threshold & dat$actual.p.value < alpha)
      sens.d <- (sum(dat$cp > threshold & dat$actual.p.value < alpha) + sum(dat$cp < threshold & dat$actual.p.value < alpha))
      sens <- sum(dat$cp > threshold & dat$actual.p.value < alpha)/(sum(dat$cp > threshold & dat$actual.p.value < alpha) + sum(dat$cp < threshold & dat$actual.p.value < alpha))
      spec.n <- sum(dat$cp < threshold & dat$actual.p.value > alpha)
      spec.d <- (sum(dat$cp < threshold & dat$actual.p.value > alpha) + sum(dat$cp > threshold & dat$actual.p.value > alpha))
      spec <- sum(dat$cp < threshold & dat$actual.p.value > alpha)/(sum(dat$cp < threshold & dat$actual.p.value > alpha) + sum(dat$cp > threshold & dat$actual.p.value > alpha))
      out <- data.frame(beta = b, n.futile = n.futile, n.type.ii.errors = n.type.II, type.ii.rate = n.type.II/n.futile,
                        n.not.futile=n.continue, n.type.i.errors = n.type.I, type.i.rate = n.type.I/n.continue,
                        n.correct.decisions=correct, accuracy=correct/nrow(dat),
                        sensitivity = sens, sens.n = sens.n, sens.d = sens.d,
                        specificity = spec, spec.n = spec.n, spec.d = spec.d,
                        youdon = sens + spec - 1)
      out.results <- rbind(out.results, out)
    }
    colnames(out.results) <- c("Beta", "Futile Trials (N)", "Type II Errors (N)", "Error Rate (%)", "Continued Trials (N)", "Type I Errors (N)", "Error Rate (%)", "Correct Decisions", "Accuracy (%)", "Sensitivity", "Sensitivity (n)", "Sensitivity (N)", "Specificity", "Specificity (n)", "Specificity (N)", "Youdon Index")
    out.results <- apply(out.results, 2, FUN=function(x)round(x, digits=digits))
    cat("Stratified Results"); cat("\n\n")
    print(noquote(out.results), row.names=FALSE)
    cat("\n\n")
  }
  
  # Printing overall classification performance results
  cat("Overall Results"); cat("\n\n")
  colnames(overall.results) <- c("Futile Trials (N)", "Type II Errors (N)", "Error Rate (%)", "Continued Trials (N)", "Type I Errors (N)", "Error Rate (%)", "Correct Decisions", "Accuracy (%)", "Sensitivity", "Sensitivity (n)", "Sensitivity (N)", "Specificity", "Specificity (n)", "Specificity (N)", "Youdon Index")
  overall.results <- apply(overall.results, 2, FUN=function(x)round(x, digits=digits))
  print(noquote(overall.results), row.names=FALSE)
}

### Visualizes accuracy of the interim futility classifications, startifying on matched/unmatched design
  # matched.data    dataframe     conditional power results from any of the simulation scripts under a matched design
  # unmatched.data  dataframe     conditional power results from any of the simulation scripts under an unmatched design
  # alpha           numeric       significance level (final analysis)
  # threshold       numeric       futility threshold (interim analysis)
visualize_accuracy_matched <- function(matched.data, unmatched.data, alpha=0.05, threshold=0.2){
  
  matched.data$matched <- 0
  unmatched.data$matched <- 1
  data <- rbind(matched.data, unmatched.data)
  n.beta <- length(unique(data$hr))
  data$actual.sig <- ifelse(data$actual.p.value < alpha, 1, 0)
  data$cp.above.t <- ifelse(data$cp > threshold, 1, 0)
  accuracy.dat <- data.frame(matched=rep(c("Matched CRT Design", "Unmatched CRT Design"), n.beta*2), setting=rep(rep(c(2, 1), each=2), n.beta),
                             beta=rep(sort(unique(log(data$hr))), each=(2*2)),
                             accuracy=as.vector(by(data, list(data$matched, data$f.theta, data$hr), FUN=function(x)mean(x$actual.sig==x$cp.above.t))))
  accuracy.dat$beta <- round(accuracy.dat$beta, 2)
  index <- which(accuracy.dat$beta %in% c(-0.15, -0.05, 0.05, 0.15))
  if (length(index) > 0) accuracy.dat <- accuracy.dat[-index, ]

  p <- ggplot(accuracy.dat, aes(x=as.factor(beta), y=as.factor(setting), fill=accuracy)) + geom_tile(colour = "white", size=0.1) + geom_text(aes(label = round(accuracy, 2)), col="white", size=8)
  p <- p + scale_fill_viridis(name="Accuracy", breaks=c(0.8, 0.9, 1), direction=-1, begin=0, end=0.9) + facet_wrap(~matched, nrow=2) + scale_x_discrete(labels=sapply(seq(min(unique(accuracy.dat$beta)), max(unique(accuracy.dat$beta)), by=0.1), FUN=function(x)format(round(x, 1), nsmall=1)))
  p <- p + theme_classic() + theme(axis.ticks.y = element_blank(), axis.text.y = element_text(size=26), axis.ticks.x = element_blank(), axis.text.x = element_text(size=24), axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 30, l = 0), size=30), legend.text=element_text(size=26), legend.title=element_text(size=30), panel.border=element_blank(), line = element_blank(), strip.background = element_blank(), strip.text.x = element_text(size = 34)) + scale_y_discrete(labels=c(bquote(atop(lambda[0]==0.01 * ",", sigma^2==0.22)), bquote(atop(lambda[0]==0.001 * ",", sigma^2==0.06)))) + labs(x=expression(paste("Data-Generating ", beta)), y=NULL)
  p
}

### Visualizes accuracy of the interim futility classifications, startifying on baseline incidence
  # data            dataframe     conditional power results from any of the simulation scripts
  # alpha           numeric       significance level (final analysis)
  # threshold       numeric       futility threshold (interim analysis)
visualize_accuracy_incidence <- function(data, alpha=0.05, threshold=0.2){
  
  n.lambda <- length(unique(data$hazard))
  n.beta <- length(unique(data$hr))
  n.var <- length(unique(data$f.theta))
  data$actual.sig <- ifelse(data$actual.p.value < alpha, 1, 0)
  data$cp.above.t <- ifelse(data$cp > threshold, 1, 0)
  accuracy.dat <- data.frame(k=rep(sort(unique(round(sqrt(exp(data$f.theta)-1), 2))), n.beta*n.lambda),
                             lambda=rep(rep(sort(unique(data$hazard)), each=n.var), n.beta),
                             beta=rep(sort(unique(log(data$hr))), each=(n.lambda*n.var)),
                             accuracy=as.vector(by(data, list(data$f.theta, data$hazard, data$hr), FUN=function(x)mean(x$actual.sig==x$cp.above.t))))
  accuracy.dat$beta <- round(accuracy.dat$beta, 2)
  index <- which(accuracy.dat$beta %in% c(-0.15, -0.05, 0.05, 0.15))
  if (length(index) > 0) accuracy.dat <- accuracy.dat[-index, ]
  accuracy.dat$lambda <- as.factor(accuracy.dat$lambda)
  levels(accuracy.dat$lambda) <- c(expression(lambda[0]==0.001), expression(lambda[0]==0.01))
  accuracy.dat$k <- factor(accuracy.dat$k, levels=c("0.5", "0.25", "0"))
  
  p <- ggplot(accuracy.dat, aes(x=as.factor(beta), y=k, fill=accuracy)) + geom_tile(colour = "white", size=0.1) + geom_text(aes(label = round(accuracy, 2)), col="white", size=8)
  p <- p + scale_fill_viridis(name="Accuracy", direction=-1, begin=0, end=0.9) + facet_wrap(~lambda, nrow=2, labeller = label_parsed) + scale_x_discrete(labels=sapply(seq(min(unique(accuracy.dat$beta)), max(unique(accuracy.dat$beta)), by=0.1), FUN=function(x)format(round(x, 1), nsmall=1)))
  p <- p + theme_classic() + theme(axis.ticks.y = element_blank(), axis.text.y = element_text(size=24), axis.ticks.x = element_blank(), axis.text.x = element_text(size=20), axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 30, l = 0), size=28), legend.text=element_text(size=24), legend.title=element_text(size=28), panel.border=element_blank(), line = element_blank(), strip.background = element_blank(), strip.text.x = element_text(size = 28)) + scale_y_discrete(labels=c(bquote(sigma^2 * " = 0.22"), bquote(sigma^2 * " = 0.06"), bquote(sigma^2 * " = 0.00"))) + labs(x=expression(paste("Data-Generating ", beta)), y=NULL)
  
  p
}

### Visualizes the distribution of conditional power estimates and the correspondence between the interim futility
### classifications and the final study significance
  # data            dataframe     conditional power results from any of the simulation scripts
  # beta            list          log hazard ratios over which to stratify the results
  # alpha           numeric       significance level (final analysis)
  # threshold       numeric       futility threshold (interim analysis)
  # titlecombined   boolean       should each column of the ploting grid have a single combined title?
  # position        character     position of the shared legend 
combined_distribution_mosaic_plots <- function(data, beta, alpha=0.05, threshold=0.2, titlecombined=TRUE, position="bottom"){
  n.plots <- length(beta)*2
  plot.list <- list()
  title.labels <- list()
  titleindividual <- !titlecombined
  
  for (b in beta){
    dat <- data[round(log(data$hr), 2) == b, ]
    dat$actual.sig <- ifelse(dat$actual.p.value > alpha, 0, 1)
    title.labels[[as.character(b)]] <- substitute(paste(beta, " = ", bvalue), list(bvalue=b))
    
    # Plotting the distribution of the conditional power estimates (stratified by the final study significance)
    p <- ggplot(dat, aes(x=cp, group=as.factor(actual.sig), fill=as.factor(actual.sig)))
    p <- p + geom_density(alpha=0.5) + geom_vline(xintercept=threshold, color="red", linetype="dashed") + scale_fill_viridis(name="Result of Completed Trial  ", labels=c("Not Significant  ", "Significant"), begin=0.5, end=0.85, discrete=TRUE)
    p <- p + theme_classic() + theme(legend.position="bottom", axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 20, l = 0), size=26), axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=26), legend.text=element_text(size=28), legend.title=element_text(size=28), plot.title = element_text(hjust = 0.5, size=30), panel.border = element_rect(colour = "black", fill=NA, size=1))
    p <- p + scale_x_continuous(name="Conditional Power") + scale_y_continuous(name="Density", limits=c(0, NA)) 
    if (titleindividual) p <- p + ggtitle(title.labels[[as.character(b)]])
    plot.list[[paste0("density", as.character(b))]] <- p
    
    # Plotting the correspondence between the interim futility classifications and final study significance
    plot.dat <- data.frame(final=c("Significant", "Not Significant", "Significant", "Not Significant"),
                           interim=c("Futile", "Futile", "Not Futile", "Not Futile"),
                           count=rep(NA, 4))
    plot.dat$count[1] <- sum(dat$actual.p.value < alpha & dat$cp < threshold)
    plot.dat$count[2] <- sum(dat$actual.p.value >= alpha & dat$cp < threshold)
    plot.dat$count[3] <- sum(dat$actual.p.value < alpha & dat$cp >= threshold)
    plot.dat$count[4] <- sum(dat$actual.p.value >= alpha & dat$cp >= threshold)
    h.just.val <- 1
    if (plot.dat$count[2]/(plot.dat$count[1] + plot.dat$count[2]) < 0.5){
      h.just.val <- 0
    }
    p <- ggplot(data=plot.dat) + geom_mosaic(aes(weight=count, x=product(interim), fill=final)) + scale_fill_viridis(name="Result of Completed Trial  ", labels=c("Not Significant  ", "Significant"), begin=0.5, end=0.85, discrete=TRUE)
    p <- p + theme(legend.position="bottom", axis.text.x = element_text(size=20), axis.ticks = element_blank(), axis.text.y = element_text(size=20, angle=90, hjust=h.just.val), axis.title.x=element_text(margin = margin(t = 20, r = 0, b = 30, l = 0), size=26), axis.title.y=element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=26), legend.text=element_text(size=28), legend.title=element_text(size=28), plot.title = element_text(hjust = 0.5, size=30), panel.background = element_blank())
    p <- p + scale_x_productlist("Interim Result") + scale_y_productlist("Final Result")
    p <- p + geom_text(data = ggplot_build(p)$data[[1]], aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label=.wt), size=6)
    if (titleindividual) p <- p + ggtitle(title.labels[[as.character(b)]])
    plot.list[[paste0("mosaic", as.character(b))]] <- p
  }
  
  # Combining the two sets of plots
  g <- ggplotGrob(plot.list[[2]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plot.list, function(x) x + theme(legend.position = "none"))
  grob.list=list()
  for (i in 1:length(beta)){
    if (titlecombined){
      grob.list[[as.character(beta[i])]] <- arrangeGrob(grobs=list(gl[[2*i-1]], gl[[2*i]]), top=textGrob(title.labels[[i]], gp=gpar(fontsize=30)), ncol=1)
    } else{
      grob.list[[as.character(beta[i])]] <- arrangeGrob(grobs=list(gl[[2*i-1]], gl[[2*i]]), ncol=1)
    }
  }
  grob.list <- c(grob.list, ncol = length(beta), nrow = 1)
  combined <- arrangeGrob(do.call(arrangeGrob, grob.list), legend, ncol = 1, heights = unit.c(unit(1, "npc") - lheight, lheight))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}


#-----------------------------------------------------------------------------------------------------
# Auxiliary Functions: Web Figures 1 & 2
#-----------------------------------------------------------------------------------------------------

### Computes the absolute difference between two survival functions at time t
  # t               numeric       time at which to compute the difference in survival 
  # curve1          function      first survival function to be compared  
  # curve2          function      second survival function to be compared
  # params1         list          list of inputs for the first survival function
  # params2         list          list of inputs for the second survival function
diff_curve <- function(t, curve1, curve2, params1, params2){
  
  params1[["t"]] <- t
  params2[["t"]] <- t
  c1 <- do.call(Vectorize(curve1), params1)
  c2 <- do.call(Vectorize(curve2), params2)
  
  return(abs(c1-c2))
}

### Returns the true cluster-conditional survival at time t assuming an exponential survival distribution
  # t               numeric       time 
  # lambda          numeric       rate parameter
cc_curve_exp <- function(t, lambda){
  exp(-lambda*t)
}

### Returns the marginal survival at time t assuming an exponential survival distribution and lognormal frailty distribution
  # t               numeric       time 
  # lambda          numeric       rate parameter
  # theta           numeric       frailty variance parameter
m_curve_exp <- function(t, lambda, theta){
  
  integrand <- function(x, t, lambda, theta){
  if (t==0){
      return(dnorm(x, 0, sd=sqrt(theta)))
    } else{
      return(exp(-lambda*t*exp(x))*dnorm(x, 0, sd=sqrt(theta)))
    }
  }
  
  integrate(integrand, lower=-Inf, upper=Inf, t=t, lambda=lambda, theta=theta, subdivisions = 2000)$value
}


### Returns the recovered cluster-conditional survival at time t assuming an exponential survival distribution and lognormal frailty distribution
  # t               numeric       time 
  # lambda          numeric       rate parameter
  # theta           numeric       frailty variance parameter
approx_cc_curve_exp <- function(t, lambda, theta){
  
  obj <- function(x){
    x*(1 + theta/2*log(x)*(log(x)+1)) - m_curve_exp(t, lambda, theta)
  }
  y <- uniroot.all(obj, interval=c(0,1))
  return(ifelse(length(y)==0, 1, y))
    
}

### Visually compares the true cluster-conditional, approximated cluster-conditional, and marginal survival functions
### assuming an exponential (cluster-conditional) survival distribution and a lognormal frailty distribution
  # lambda          numeric       rate parameter
  # theta           numeric       frailty variance parameter
  # stop            numeric       upper limit of the time range for the curve comparison
  # title           boolean       should a title be produced for the plot?
curve_comparison_exp <- function(lambda, theta, stop, title=FALSE){
  
  cc.dat <- data.frame(surv = sapply(seq(0, stop, by=1), FUN=function(x)cc_curve_exp(x, lambda)),
                       time = seq(0, stop, by=1),
                       curve = "cluster conditional")
  m.dat <- data.frame(surv = sapply(seq(0, stop, by=1), FUN=function(x)m_curve_exp(x, lambda, theta)),
                      time = seq(0, stop, by=1),
                      curve = "marginal")
  approx.cc.dat <- data.frame(surv = sapply(seq(0, stop, by=1), FUN=function(x)approx_cc_curve_exp(x, lambda, theta)),
                              time = seq(0, stop, by=1),
                              curve="approximate cluster conditional")
  plot.dat <- rbind(cc.dat, m.dat, approx.cc.dat)
  
  m.diff <- integrate(diff_curve, lower=0, upper=stop, curve1=cc_curve_exp, curve2=m_curve_exp, params1=list("lambda"=lambda), params2=list("lambda"=lambda, "theta"=theta), subdivisions=2000)$value
  approx.diff <- integrate(diff_curve, lower=0, upper=stop, curve1=cc_curve_exp, curve2=approx_cc_curve_exp, params1=list("lambda"=lambda), params2=list("lambda"=lambda, "theta"=theta), subdivisions=2000)$value
  m.diff <- format(round(m.diff, 3), nsmall = 3)
  approx.diff <- format(round(approx.diff, 3), nsmall = 3)
  
  colpal <- viridis(2, begin=0.5, end=0.85)
  y.stop <- ifelse(lambda==0.01, 0, 0.77)
  y.2 <- ifelse(lambda==0.01, 0.96, 0.9905)
  y.3 <- ifelse(lambda==0.01, 0.9053, 0.9775)
  y.min <- ifelse(lambda==0.01, 0.867, 0.9685)
  
  annotation.label.1 <- substitute(paste("d(S,",bar(S), ") :    ", diff), list(diff=m.diff))
  annotation.label.2 <- substitute(paste("d(S,",hat(S), ") :    ", diff), list(diff=approx.diff))
  
  p <- ggplot(data=plot.dat, aes(x=time, y=surv, color=curve))
  p <- p + geom_line(aes(linetype=curve, color=curve), size=1) + scale_linetype_manual(values=c(1, 1, 2)) + guides(linetype=FALSE)
  p <- p + scale_color_manual(name=" ", labels=c("S(t|X=0)", expression(paste(bar(S), "(t|X=0)")),  expression(paste(hat(S), "(t|X=0)"))), values=c(colpal[2], colpal[1], "black"))
  p <- p + theme_classic() + theme(legend.position = c(0.15, 0.12), legend.text.align = 0, legend.text=element_text(size=18), legend.key.size = unit(1.8, 'lines'), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.x=element_text(margin=margin(t = 15, r = 0, b = 15, l = 0), size=22), axis.title.y=element_text(margin=margin(t = 0, r = 10, b = 0, l = 0), size=22), plot.title = element_text(margin=margin(t = 0, r = 0, b = 15, l = 0), hjust = 0.5, size=28), panel.border = element_rect(colour = "black", fill=NA, size=1))
  p <- p + scale_x_continuous(name="Time (in Weeks)") + scale_y_continuous(name="Survival", limits=c(y.stop, 1))
  p <- p + 
    annotate('text', x = 133, y = y.2, 
             label = deparse(annotation.label.1), parse = TRUE, hjust = 0, size = 6) +
    annotate('text', x = 133, y = y.3,
             label = deparse(annotation.label.2), parse = TRUE, hjust = 0, size = 6) +
    annotate('rect', xmin=127, xmax=205, ymin=y.min, ymax=1, col="dark gray", alpha=0)
  if (title){
    k <- format(round(sqrt(exp(theta)-1), 2), nsmall=2)
    p <- p + ggtitle(substitute(paste("k=", val), list(val=k)))
  }
  
  return(p)
}

### Selects Weibull parameters so that the survival at t weeks matches that for a given exponential distribution
  # kappa           numeric       shape parameter
  # exp.lambda      numeric       rate parameter of the corresponding exponential distribution
  # t               numeric       time point at which to match the survival
get_param_lambda <- function(kappa, exp.lambda, t){
  return(t*(-log(1-pexp(t, exp.lambda)))^(-1/kappa))
}

### Returns the true cluster-conditional survival at time t assuming a Weibull survival distribution
  # t               numeric       time 
  # lambda          numeric       scale parameter
  # kappa           numeric       shape parameter
cc_curve_weib <- function(t, lambda, kappa){
  exp(-(t/lambda)^(kappa))
}

### Returns the marginal survival at time t assuming a Weibull survival distribution and lognormal frailty distribution
  # t               numeric       time 
  # lambda          numeric       scale parameter
  # kappa           numeric       shape parameter
  # theta           numeric       frailty variance parameter
m_curve_weib <- function(t, lambda, kappa, theta){
  
  integrand <- function(x, t, lambda, kappa, theta){
    if (t==0){
      return(dnorm(x, 0, sd=sqrt(theta)))
    } else{
      return(exp(-(t/lambda)^(kappa)*exp(x))*dnorm(x, 0, sd=sqrt(theta)))
    }
  }
  integrate(integrand, lower=-Inf, upper=Inf, t=t, lambda=lambda, kappa=kappa, theta=theta, subdivisions=2000)$value
}


### Returns the recovered cluster-conditional survival at time t assuming a Weibull survival distribution and lognormal frailty distribution
  # t               numeric       time 
  # lambda          numeric       scale parameter
  # kappa           numeric       shape parameter
  # theta           numeric       frailty variance parameter
approx_cc_curve_weib <- function(t, lambda, kappa, theta){
  
  obj <- function(x){
    x*(1 + theta/2*log(x)*(log(x)+1)) - m_curve_weib(t, lambda, kappa, theta)
  }
  y <- uniroot.all(obj, interval=c(0,1))
  return(ifelse(length(y)==0, 1, y))
}

### Visually compares the true cluster-conditional, approximated cluster-conditional, and marginal survival functions
### assuming a Weibull (cluster-conditional) survival distribution and a lognormal frailty distribution
  # lambda          numeric       scale parameter
  # kappa           numeric       shape parameter
  # theta           numeric       frailty variance parameter
  # stop            numeric       upper limit of the time range for the curve comparison
  # title           boolean       should a title be produced for the plot?
curve_comparison_weib <- function(lambda, kappa, theta, stop, title=FALSE){
  
  cc.dat <- data.frame(surv = sapply(seq(0, stop, by=1), FUN=function(x)cc_curve_weib(x, lambda, kappa)),
                       time = seq(0, stop, by=1),
                       curve = "cluster conditional")
  m.dat <- data.frame(surv = sapply(seq(0, stop, by=1), FUN=function(x)m_curve_weib(x, lambda, kappa, theta)),
                      time = seq(0, stop, by=1),
                      curve = "marginal")
  approx.cc.dat <- data.frame(surv = sapply(seq(0, stop, by=1), FUN=function(x)approx_cc_curve_weib(x, lambda, kappa, theta)),
                              time = seq(0, stop, by=1),
                              curve="approximate cluster conditional")
  plot.dat <- rbind(cc.dat, m.dat, approx.cc.dat)
  
  m.diff <- integrate(diff_curve, lower=0, upper=stop, curve1=cc_curve_weib, curve2=m_curve_weib, params1=list("lambda"=lambda, "kappa"=kappa), params2=list("lambda"=lambda, "kappa"= kappa, "theta"=theta), subdivisions=2000)$value
  approx.diff <- integrate(diff_curve, lower=0, upper=stop, curve1=cc_curve_weib, curve2=approx_cc_curve_weib, params1=list("lambda"=lambda, "kappa"=kappa), params2=list("lambda"=lambda, "kappa"=kappa, "theta"=theta), subdivisions=2000)$value
  m.diff <- format(round(m.diff, 3), nsmall = 3)
  approx.diff <- format(round(approx.diff, 3), nsmall = 3)
  
  colpal <- viridis(2, begin=0.5, end=0.85)
  y.stop <- ifelse(lambda < 200, 0, 0.77)
  y.2 <- ifelse(lambda < 200, 0.96, 0.9905)
  y.3 <- ifelse(lambda < 200, 0.9053, 0.9775)
  y.min <- ifelse(lambda < 200, 0.867, 0.9685)
  
  annotation.label.1 <- substitute(paste("d(S,",bar(S), ") :    ", diff), list(diff=m.diff))
  annotation.label.2 <- substitute(paste("d(S,",hat(S), ") :    ", diff), list(diff=approx.diff))
  
  p <- ggplot(data=plot.dat, aes(x=time, y=surv, color=curve))
  p <- p + geom_line(aes(linetype=curve, color=curve), size=1) + scale_linetype_manual(values=c(1, 1, 2)) + guides(linetype=FALSE)
  p <- p + scale_color_manual(name=" ", labels=c("S(t|X=0)", expression(paste(bar(S), "(t|X=0)")),  expression(paste(hat(S), "(t|X=0)"))), values=c(colpal[2], colpal[1], "black"))
  p <- p + theme_classic() + theme(legend.position = c(0.15, 0.12), legend.text.align = 0, legend.text=element_text(size=18), legend.key.size = unit(1.8, 'lines'), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.x=element_text(margin=margin(t = 15, r = 0, b = 15, l = 0), size=22), axis.title.y=element_text(margin=margin(t = 0, r = 10, b = 0, l = 0), size=22), plot.title = element_text(margin=margin(t = 0, r = 0, b = 15, l = 0), hjust = 0.5, size=28), panel.border = element_rect(colour = "black", fill=NA, size=1))
  p <- p + scale_x_continuous(name="Time (in Weeks)") + scale_y_continuous(name="Survival", limits=c(y.stop, 1))
  p <- p + 
    annotate('text', x = 133, y = y.2, 
             label = deparse(annotation.label.1), parse = TRUE, hjust = 0, size = 6) +
    annotate('text', x = 133, y = y.3,
             label = deparse(annotation.label.2), parse = TRUE, hjust = 0, size = 6) +
    annotate('rect', xmin=127, xmax=205, ymin=y.min, ymax=1, col="dark gray", alpha=0)
  if (title){
    k <- format(round(sqrt(exp(theta)-1), 2), nsmall=2)
    p <- p + ggtitle(substitute(paste("k=", val), list(val=k)))
  }
  
  return(p)
}
