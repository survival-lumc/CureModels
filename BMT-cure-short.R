#' This file contains the code supporting Section 2 of the paper
#' "Cure models in survival analysis: Issues with identifiability and
#' competing risks" by Hein Putter and Per Kragh Andersen.
#' 
#' It uses the BMT data as used in the book by Andersen & Ravn (2023),
#' available at "https://multi-state-book.github.io/companion/data/bmt.csv".

# Packages used in this analysis
library(survival)
library(prodlim)
library(smcure)
library(intsurv)

# Read in the data
bmt <- read.csv("https://multi-state-book.github.io/companion/data/bmt.csv")
# Have a look at the events of interest, here relapse and death
head(bmt)
table(bmt$death)
table(bmt$rel)
table(bmt$rel, bmt$death)
# The only covariates that we will use are all (selecting AML, all=0) and age
table(bmt$all)
summary(bmt$age)
hist(bmt$age)

#' A cure model would focus on relapse and fit a cure model where subjects that
#' died are treated as censored observations. We favour a competing risks model
#' where death without prior relapse is considered as a competing risk.
#' 
#' First prepare the data for a competing risks analysis by defining time and
#' status variables.
bmt$relcr <- bmt$rel
bmt$relcr[bmt$death == 1 & bmt$rel==0] <- 2
table(bmt$relcr)
bmt$timecr <- bmt$timedeath
bmt$timecr[bmt$relcr == 1] <- bmt$timerel[bmt$relcr == 1]
#' We will use age very crudely as above or below 30 years. Define this variable,
#' select the AML patients, retain only the relevant variables and look at the
#' first six lines.
bmt$Age <- cut(bmt$age, c(0, 30, 100))
aml <- subset(bmt, all==0)
table(aml$Age, aml$relcr)

#' Start with the incorrect naive Kaplan-Meiers for relapse.
kmfit <- prodlim(Hist(timecr, rel) ~ Age, data = aml)
plot(kmfit,
     xlim = c(0, 144), xlab = "Months since transplantation", ylim = c(0.6, 1),
     axis1.at = seq(0, 144, 24), axis2.at = seq(0.6, 1, 0.1),
     atrisk.at = seq(0, 144, 24), atrisk.title = "",
     legend.x = "topright", legend.cex = 1.25)
#' There are the plateaus of the KM's.
summary(kmfit, times=120)

#' Fit a cure model with `smcure`.
curefit <- smcure(formula = Surv(timecr, rel) ~ Age, 
                  cureform = ~ Age, data = aml, 
                  model = "ph", Var = TRUE)

#' I tried the same thing now using another package that can fit
#' cure models, called `intsurv`. Adapting some code from the help
#' from `intsurv`, I was able to obtain much more plausible
#' estimates, at least of the cure probability model.
mm <- model.matrix(~ Age, data = aml)[, -1, drop=FALSE] # without intercept
set.seed(2024) # To get reproducible results for the variance (using bootstrap)
fit1 <- cox_cure.fit(mm, mm, aml$timecr, aml$rel, bootstrap = 100)
summary(fit1)
#' Checking whether the KM plateaus agree with the cure probabilities.
b <- coef(fit1)$cure
1 - plogis(b[1]) # very close to the KM plateau of AML
1 - plogis(b[1] + b[2]) # very close to the KM plateau of ALL

#' Let's stick with the competing risks analysis and plot the non-parametric
#' estimates of the cumulative incidence curves, using prodlim.
ajfit <- prodlim(Hist(timecr, relcr) ~ Age, data = aml)
plot(ajfit, # relapse
     xlim = c(0, 144), xlab = "Months since transplantation", ylim = c(0, 0.4),
     axis1.at = seq(0, 144, 24), axis2.at = seq(0, 0.4, 0.1),
     atrisk.at = seq(0, 144, 24), atrisk.title = "",
     legend.x = "topright", legend.cex = 1.25)
title(main = "Relapse")
plot(ajfit, cause=2, # death in remission
     xlim = c(0, 144), xlab = "Months since transplantation", ylim = c(0, 0.4),
     axis1.at = seq(0, 144, 24), axis2.at = seq(0, 0.4, 0.1),
     atrisk.at = seq(0, 144, 24), atrisk.title = "",
     legend.x = "bottomright", legend.cex = 1.25)
title(main = "Death in remission")

#' The following fits proportional cause-specific hazards models in the AML
#' data with age (categorized older or younger than 30) as sole covariate.
# Cause 1, relapse
coxph(Surv(timecr, relcr==1) ~ Age, data = aml)
# Cause 2, non-relapse mortality
coxph(Surv(timecr, relcr==2) ~ Age, data = aml)

