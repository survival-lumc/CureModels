library(simsurv)
library(intsurv)
library(tictoc)
library(dynpred)
library(progress)

logit <- function(p) log(p / (1-p))
expit <- function(x) exp(x) / (1 + exp(x))

gd <- function(n, tau, lambda, xdist=c("binary", "uniform")) {

  b <- function(t) log(2) - lambda * log(1 + exp(2 * t))
  hazard <- function(t, x, betas) {
    exp(x[["x"]] * b(t))
  }

  gam0 <- 0
  gam1 <- -log(2)
  if (xdist == "binary")
    x <- rbinom(n, size = 1, prob = 0.5)
  else
    x <- runif(n)
  tt <- cured <- rep(NA, n)
  pcure <- expit(gam0 + gam1*x)
  cured <- rbinom(n, 1, prob=pcure)

  dfr <- data.frame(id = 1:n, x = x)
  sdat <- simsurv(hazard = hazard, x = dfr, maxt = 100, interval = c(1e-8, 101))
  dfr <- merge(dfr, sdat, sort = FALSE)
  names(dfr)[names(dfr) == "eventtime"] <- "time" # same names as dfr earlier on
  # Add cure to the data
  dfr$time[cured==1] <- tau
  dfr$status[cured==1] <- 0
  return(dfr)
}

###
# Running simulation for binary x
###

M <- 1000
n <- 1000

set.seed(12345)

lambdaseq <- seq(0, 1, by=0.2)
nlambdas <- length(lambdaseq)

tic("All lambdas")
res <- NULL
for (i in 1:nlambdas) {

  lambda <- lambdaseq[i]
  deb(lambda, method = "cat")

  resi <- matrix(NA, M, 3)
  
  pb <- progress_bar$new(total = M)
  
  for (m in 1:M) {
    pb$tick()
    dd <- gd(n=1000, tau=100, lambda = lambda, xdist="binary")
    # Fit with smcure (gives same results as cox_cure below, but much slower)
    # curefit <- smcure(formula = Surv(time, status) ~ x,
    #                   cureform = ~ x, data = dd,
    #                   model = "ph", Var = FALSE)
    # resi[m, 1:2] <- -curefit$b
    # resi[m, 3] <- curefit$beta
    # Fit with cox_cure from intsurv package
    fit1 <- cox_cure(surv_formula = ~ x, cure_formula = ~ x,
                     time = time, event = status, data = dd)
    resi[m, 1:2] <- -coef(fit1)$cure
    resi[m, 3] <- coef(fit1)$surv
  }
  resi <- cbind(resi, rep(lambda, M))
  res <- rbind(res, resi)
}
toc()

save(res, file="resbinary.Rdata")

###
# Running simulation for uniform x
###

M <- 1000
n <- 1000

set.seed(12345)

lambdaseq <- seq(0, 1, by=0.2)
nlambdas <- length(lambdaseq)

tic("All lambdas")
res <- NULL
for (i in 1:nlambdas) {

  lambda <- lambdaseq[i]
  deb(lambda, method = "cat")
  
  resi <- matrix(NA, M, 3)
  pb <- progress_bar$new(total = M)
  
  for (m in 1:M) {
    pb$tick()
    dd <- gd(n=1000, tau=100, lambda = lambda, xdist="uniform")
    # Fit with smcure (gives same results as cox_cure below, but much slower)
    # curefit <- smcure(formula = Surv(time, status) ~ x,
    #                   cureform = ~ x, data = dd,
    #                   model = "ph", Var = FALSE)
    # resi[m, 1:2] <- -curefit$b
    # resi[m, 3] <- curefit$beta
    # Fit with cox_cure from intsurv package
    fit1 <- cox_cure(surv_formula = ~ x, cure_formula = ~ x,
                     time = time, event = status, data = dd)
    resi[m, 1:2] <- -coef(fit1)$cure
    resi[m, 3] <- coef(fit1)$surv
  }
  resi <- cbind(resi, rep(lambda, M))
  res <- rbind(res, resi)
}
toc()

save(res, file="resuniform.Rdata")

##### Reading and plotting results #####

library(tidyverse)
library(ggplot2)

### Binary x
load("resbinary.Rdata")

res <- as.data.frame(res)
names(res) <- c("gamma0", "gamma1", "logHR", "lambda")
res$lambda <- as.factor(res$lambda)
res <- as_tibble(res)

# gamma_0
ggp <- res %>% group_by(lambda) %>% 
  ggplot(aes(x = lambda, y = gamma0)) +
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=0.25, dotsize=0.25) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=3, color="red") +
  labs(title="gamma_0", x="Lambda", y = "Coefficient") + 
  theme_bw()
ggp

# gamma_1
ggp <- res %>% group_by(lambda) %>% 
  ggplot(aes(x = lambda, y = gamma1)) +
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=0.1, dotsize=0.25) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=3, color="red") +
  labs(title="gamma_1", x="Lambda", y = "Coefficient") + 
  theme_bw()
ggp
ggsave("gamma1_binary.pdf", device = "pdf")

ggp <- res %>% group_by(lambda) %>% 
  ggplot(aes(x = lambda, y = logHR)) +
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=0.15, dotsize=0.25) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=3, color="red") +
  labs(title="Log hazard ratio", x="Lambda", y = "Coefficient") + 
  theme_bw()
ggp

### Uniform x
load("resuniform.Rdata")

res <- as.data.frame(res)
names(res) <- c("gamma0", "gamma1", "logHR", "lambda")
res$lambda <- as.factor(res$lambda)
res <- as_tibble(res)

ggp <- res %>% group_by(lambda) %>% 
  ggplot(aes(x = lambda, y = gamma0)) +
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=0.25, dotsize=0.25) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=3, color="red") +
  labs(title="gamma_0", x="Lambda", y = "Coefficient") + 
  theme_bw()
ggp

ggp <- res %>% group_by(lambda) %>% 
  ggplot(aes(x = lambda, y = gamma1)) +
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=0.1, dotsize=0.25) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=3, color="red") +
  labs(title="gamma_1", x="Lambda", y = "Coefficient") + 
  theme_bw()
ggp
ggsave("gamma1_uniform.pdf", device = "pdf")

ggp <- res %>% group_by(lambda) %>% 
  ggplot(aes(x = lambda, y = logHR)) +
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=0.15, dotsize=0.25) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=3, color="red") +
  labs(title="Log hazard ratio", x="Lambda", y = "Coefficient") + 
  theme_bw()
ggp
