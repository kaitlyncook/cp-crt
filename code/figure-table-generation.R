#-----------------------------------------------------------------------------------------------------#
#   Simulation Studies: Producing the Figures and Tables in Cook and Wang (2020)                      #
#-----------------------------------------------------------------------------------------------------#

source("code/analysis-functions.R")


#-----------------------------------------------------------------------------------------------------
# Reading in the Simulation Results
#-----------------------------------------------------------------------------------------------------

## Unmatched Design Results

unmatched.low.ind <- read.csv("data/simulation-results/unmatched_0.001_0.csv")
unmatched.low.mild <- read.csv("data/simulation-results/unmatched_0.001_0.06.csv")
unmatched.low.mod <- read.csv("data/simulation-results/unmatched_0.001_0.22.csv")
unmatched.low <- rbind(unmatched.low.ind, unmatched.low.mild, unmatched.low.mod)

unmatched.high.ind <- read.csv("data/simulation-results/unmatched_0.01_0.csv")
unmatched.high.mild <- read.csv("data/simulation-results/unmatched_0.01_0.06.csv")
unmatched.high.mod <- read.csv("data/simulation-results/unmatched_0.01_0.22.csv")
unmatched.high <- rbind(unmatched.high.ind, unmatched.high.mild, unmatched.high.mod)


## Matched Design Results

matched.low <- read.csv("data/simulation-results/matched_0.001_0.06.csv")
matched.high <- read.csv("data/simulation-results/matched_0.01_0.22.csv")


## Sensitivity Analysis Results

moderate.m <- read.csv("data/simulation-results/moderate_m_0.001_0.06.csv")
large.m <- read.csv("data/simulation-results/large_m_0.001_0.06.csv")
frequency <- read.csv("data/simulation-results/freq_0.001_0.06.csv")
re <- read.csv("data/simulation-results/re_0.001_0.25.csv")


#-----------------------------------------------------------------------------------------------------
# Main Manuscript: Figures & Tables
#-----------------------------------------------------------------------------------------------------

## Figure 3: Conditional Power Curves

# Matched CRT design
matched.l <- format_plot_data(matched.low[abs(log(matched.low$hr)) < 0.9, ])
matched.h <- format_plot_data(matched.high)
p1 <- visualize_cp_curves(matched.l)
p2 <- visualize_cp_curves(matched.h)
matched.p <- grid.arrange(p1, p2, nrow=1, ncol=2, top = textGrob("Matched CRT Design", gp=gpar(fontsize=26)))

# Unmatched CRT design
unmatched.l <- format_plot_data(unmatched.low[unmatched.low$f.theta==0.06 & abs(log(unmatched.low$hr)) < 0.9, ])
unmatched.h <- format_plot_data(unmatched.high[unmatched.high$f.theta==0.22, ])
p3 <- visualize_cp_curves(unmatched.l)
p4 <- visualize_cp_curves(unmatched.h)
unmatched.p <- grid.arrange(p3, p4, nrow=1, ncol=2, top = textGrob("Unmatched CRT Design", gp=gpar(fontsize=26)))

figure.3 <- grid.arrange(matched.p, unmatched.p, nrow=2, ncol=1)
ggsave("figure-3-cp.pdf", figure.3, width=18, height=17, units="in")


## Figure 4: Classification Accuracy

figure.4 <- visualize_accuracy_matched(rbind(matched.low, matched.high), 
                                       rbind(unmatched.low[unmatched.low$f.theta==0.06, ], 
                                             unmatched.high[unmatched.high$f.theta==0.22,]))
ggsave("figure-4-accuracy.pdf", plot=figure.4, width=24, height=16, unit="in")


## Summary of Classification Performance

# Matched CRT design
summarize_accuracy(rbind(matched.low[!(round(log(matched.low$hr), 2) %in% c(-0.15, -0.05, 0.05, 0.15)),],
                         matched.high[!(round(log(matched.high$hr), 2) %in% c(-0.15, -0.05, 0.05, 0.15)), ]))

# Unmatched CRT design
summarize_accuracy(rbind(unmatched.low[unmatched.low$f.theta==0.06 & !(round(log(unmatched.low$hr), 2) %in% c(-0.15, -0.05, 0.05, 0.15)),],
                         unmatched.high[unmatched.high$f.theta==0.22 & !(round(log(unmatched.high$hr), 2) %in% c(-0.15, -0.05, 0.05, 0.15)), ]))

# Overall
summarize_accuracy(rbind(matched.low[!(round(log(matched.low$hr), 2) %in% c(-0.15, -0.05, 0.05, 0.15)),],
                         matched.high[!(round(log(matched.high$hr), 2) %in% c(-0.15, -0.05, 0.05, 0.15)), ],
                         unmatched.low[unmatched.low$f.theta==0.06 & !(round(log(unmatched.low$hr), 2) %in% c(-0.15, -0.05, 0.05, 0.15)),],
                         unmatched.high[unmatched.high$f.theta==0.22 & !(round(log(unmatched.high$hr), 2) %in% c(-0.15, -0.05, 0.05, 0.15)), ]))


#-----------------------------------------------------------------------------------------------------
# Supporting Information: Web Figures & Tables
#-----------------------------------------------------------------------------------------------------

## Web Figures 1 & 2: Performance of the Approximated Cluster-Conditional Survival Function

# Exponential failure times
el1 <- curve_comparison_exp(0.001, 0.22, 208)
el2 <- curve_comparison_exp(0.001, 0.45, 208)
el3 <- curve_comparison_exp(0.001, 0.69, 208)
eh1 <- curve_comparison_exp(0.01, 0.22, 208)
eh2 <- curve_comparison_exp(0.01, 0.45, 208)
eh3 <- curve_comparison_exp(0.01, 0.69, 208)
web.figure.1 <- grid.arrange(arrangeGrob(grobs=list(el1, eh1), top=textGrob("k=0.50", gp=gpar(fontsize=25)), ncol=1),
                             arrangeGrob(grobs=list(el2, eh2), top=textGrob("k=0.75", gp=gpar(fontsize=25)), ncol=1),
                             arrangeGrob(grobs=list(el3, eh3), top=textGrob("k=1.00", gp=gpar(fontsize=25)), ncol=1),
                             ncol=3)

# Weibull failure times
wl1 <- curve_comparison_weib(get_param_lambda(0.5, 0.001, 208), 0.5, 0.22, 208)
wl2 <- curve_comparison_weib(get_param_lambda(0.5, 0.001, 208), 0.5, 0.45, 208)
wl3 <- curve_comparison_weib(get_param_lambda(0.5, 0.001, 208), 0.5, 0.69, 208)
wh1 <- curve_comparison_weib(get_param_lambda(2, 0.01, 208), 2, 0.22, 208)
wh2 <- curve_comparison_weib(get_param_lambda(2, 0.01, 208), 2, 0.45, 208)
wh3 <- curve_comparison_weib(get_param_lambda(2, 0.01, 208), 2, 0.69, 208)
web.figure.2 <- grid.arrange(arrangeGrob(grobs=list(wl1, wh1), top=textGrob("k=0.50", gp=gpar(fontsize=25)), ncol=1),
                             arrangeGrob(grobs=list(wl2, wh2), top=textGrob("k=0.75", gp=gpar(fontsize=25)), ncol=1),
                             arrangeGrob(grobs=list(wl3, wh3), top=textGrob("k=1.00", gp=gpar(fontsize=25)), ncol=1),
                             ncol=3)

ggsave("figure-1-exp-curve.png", plot=web.figure.1, width=20, height=16, units="in")
ggsave("figure-2-weib-curve.png", plot=web.figure.2, width=20, height=16, units="in")


## Web Figures 3 & 4: Individual Study-Specific Concordance

web.figure.3 <- combined_distribution_mosaic_plots(matched.low, beta=list(-0.1, -0.2, -0.3))
web.figure.4 <- combined_distribution_mosaic_plots(unmatched.low[unmatched.low$f.theta==0.06,], beta=list(-0.1, -0.2, -0.3))

ggsave("figure-3-matched-concordance.pdf", plot=web.figure.3, width=22, height=16, unit="in")
ggsave("figure-4-unmatched-concordance.pdf", plot=web.figure.4, width=22, height=16, unit="in")


## Web Figures 5 & 6: Data Generation (Bias + MSE)

web.figure.5 <- visualize_data_gen(rbind(matched.low, matched.high), beta=c(-0.2, 0, 0.2))
web.figure.6 <- visualize_data_gen(rbind(unmatched.low[unmatched.low$f.theta==0.06,],
                                         unmatched.high[unmatched.high$f.theta==0.22, ]), beta=c(-0.2, 0, 0.2))

ggsave("figure-5-matched-data-generation.pdf", plot=web.figure.5, width=20, height=16, units="in")
ggsave("figure-6-unmatched-data-generation.pdf", plot=web.figure.6, width=20, height=16, units="in")


## Web Table 2: Data Generation

summarize_data_gen(rbind(matched.low, matched.high), beta=c(-0.2, 0, 0.2))
summarize_data_gen(rbind(unmatched.low[unmatched.low$f.theta==0.06, ], 
                         unmatched.high[unmatched.high$f.theta==0.22, ]), beta=c(-0.2, 0, 0.2))


## Web Figure 7: Full Conditional Power Curves (Unmatched CRT Design)

l.i  <- format_plot_data(unmatched.low[unmatched.low$f.theta==0 & abs(log(unmatched.low$hr)) < 0.5, ], TRUE)
l.mi <- format_plot_data(unmatched.low[unmatched.low$f.theta==0.06 & abs(log(unmatched.low$hr)) < 0.9, ], TRUE)
l.mo <- format_plot_data(unmatched.low[unmatched.low$f.theta==0.22, ], TRUE)
h.i  <- format_plot_data(unmatched.high[unmatched.high$f.theta==0 & abs(log(unmatched.high$hr)) < 0.3, ], TRUE)
h.mi <- format_plot_data(unmatched.high[unmatched.high$f.theta==0.06 & abs(log(unmatched.high$hr)) < 0.7, ], TRUE)
h.mo <- format_plot_data(unmatched.high[unmatched.high$f.theta==0.22, ], TRUE)

cp1 <- visualize_cp_curves(l.i)
cp2 <- visualize_cp_curves(l.mi)
cp3 <- visualize_cp_curves(l.mo)
cp4 <- visualize_cp_curves(h.i)
cp5 <- visualize_cp_curves(h.mi)
cp6 <- visualize_cp_curves(h.mo)

web.figure.7 <- grid.arrange(cp1, cp4, cp2, cp5, cp3, cp6, nrow=3, ncol=2)
ggsave("figure-7-cp.pdf", web.figure.7, width=18.5, height=24, units="in")


## Web Figure 8: Classification Accuracy (Unmatched CRT Design)

web.figure.8 <- visualize_accuracy_incidence(rbind(unmatched.low, unmatched.high))
ggsave("figure-8-accuracy.pdf", plot=web.figure.8, width=22, height=16, unit="in")


## Web Table 3: Cluster Size Sensitivity Analysis

setting.1 <- format_plot_data(moderate.m, TRUE) 
setting.2 <- format_plot_data(large.m, TRUE) 


## Web Table 4: Frequency of Monitoring Visits Sensitivity Analysis

temp <- format_plot_data(frequency, TRUE)


## Web Table 5: Misspecification of the Frailty Distribution Sensitivity Analysis

# Table creation
temp <- matrix(NA, ncol=5, nrow=3); i <- 1
for (h in sort(unique(re$hr))){
  dat <- re[re$hr==h, ]
  results.mean <- apply(cbind(dat$power, dat$ln.cp, dat$gamma.cp, dat$known.cp), 2, mean)
  results.sd <- apply(cbind(dat$power, dat$ln.cp, dat$gamma.cp, dat$known.cp), 2, sd)
  temp[i, 1] <- round(log(h), 3)
  temp[i, 2:5] <- c(paste0(round(results.mean[1], 3), " (", round(results.sd[1], 3), ")"),
                    paste0(round(results.mean[2], 3), " (", round(results.sd[2], 3), ")"), 
                    paste0(round(results.mean[3], 3), " (", round(results.sd[3], 3), ")"),
                    paste0(round(results.mean[4], 3), " (", round(results.sd[4], 3), ")"))
  i <- i + 1
}
colnames(temp) <- c("Beta", "Power", "CP (Lognormal)", "CP (Gamma)", "CP (All Known)")
print(noquote(temp), row.names=FALSE)

# Correspondence between interim futility classifications under gamma and lognormal frailties
gamma.futile <- ifelse(re$gamma.cp < 0.2, 1, 0)
ln.futile <- ifelse(re$ln.cp < 0.2, 1, 0)
sum(gamma.futile==ln.futile)
mean(gamma.futile==ln.futile)

