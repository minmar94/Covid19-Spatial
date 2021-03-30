
# Packages ----------------------------------------------------------------

require(rstan)
require(tidyverse)
require(magrittr)
require(lubridate)
require(zeallot)
require(parallel)
require(mvtnorm)
require(ggforce)
source("OtherFuns/DrichFuns.R")
source("OtherFuns/DataFuns.R")
source("OtherFuns/StanFuns.R")
set.seed(130494)

# Reading Data ------------------------------------------------------------

load("Data/DatiStan.RData")

# Stan model compilation --------------------------------------------------

mc.cores = parallel::detectCores()
# Stan options
rstan_options(auto_write = TRUE)

stan_Multip1 <- stan_model("IndTRnd_MultipAll.stan")
stan_Multip2 <- stan_model("IndTRnd2_MultipAll.stan")

# Data preparation for New Positives--------------------------------------------------------

Y <- as.integer(dati$NP)
N <- as.integer(length(Y))
Nreg <- as.integer(length(unique(dati$denominazione_regione)))
timeIdx <- as.integer(dati$WW)+1
regIdx <- as.integer(droplevels(dati$denominazione_regione))
t <- unique(dati$WW)
Ntimes <- as.integer(length(t))

X1 <- model.matrix(~.-1, data=dati %>% dplyr::select(Swabs) %>% mutate(Swabs=scale(Swabs)))
k1 <- ncol(X1)
lOff1 <- rep(logE, each=Ntimes)

lOff2 <- log(dati %$% Swabs + 1)

# Useful data summaries
ySums <-  dati %>% group_by(denominazione_regione) %>% summarise(y=sum(NP)) %$% y
yMins <-  dati %>% group_by(denominazione_regione) %>% summarise(y=min(NP)) %$% y
tflex <- dati %>% group_by(denominazione_regione) %>% summarise(WW=WW[which.max(NP)]) %$% WW

# Prepare data
dat1 <- list(
  "N" = N, 
  "Y" = Y,
  "Ntimes" = Ntimes,
  "t" = t,
  "tId" = timeIdx,
  "Nreg" = Nreg, 
  "rId" = regIdx,
  "lOff" = lOff1,
  "k" = k1,
  "X" = X1)
dat2 <- list(
  "N" = N, 
  "Y" = Y,
  "Ntimes" = Ntimes,
  "t" = t,
  "tId" = timeIdx,
  "Nreg" = Nreg, 
  "rId" = regIdx,
  "lOff" = lOff2)

# Stan Settings -----------------------------------------------------------

# Chains
n_chains <- 2
M <- 10000
n_cores <- (mc.cores)

# Define a function to generate initial values

init <- function(chain_id = 1)
{
  list(# Richards
    logr = rnorm(Nreg, log(ySums)-logE, 1), 
    logh = rnorm(Nreg, log(.5), .1),
    p = rnorm(Nreg, tflex, 5), 
    logs = rnorm(Nreg, log(1), .1),
    # Baseline
    logbase = rnorm(Nreg, log(yMins+0.0000001)-logE, .1))
} 


fit_Stan1 <- sampling(stan_Multip1, data = dat1, chains = n_chains, iter = M, 
                      cores = n_cores, init=init,
                      control = list(adapt_delta = 0.9, max_treedepth = 12)
)
fit_Stan2 <- sampling(stan_Multip2, data = dat2, chains = n_chains, iter = M, 
                      cores = n_cores, init=init,
                      control = list(adapt_delta = 0.9, max_treedepth = 12)
)
c(postsamples_Stan1, ypreds_Stan1, ypredsQ_Stan1) %<-% extract_postY(fit_Stan1)
c(postsamples_Stan2, ypreds_Stan2, ypredsQ_Stan2) %<-% extract_postY(fit_Stan2)
# Diagnostic check: yellow or red points are bad!

save.image(file="WS/IndTRndMultipAll_TrueData.RData")
