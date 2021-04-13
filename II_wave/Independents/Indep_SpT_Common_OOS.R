
# Packages ----------------------------------------------------------------

require(rstan)
require(tidyverse)
require(magrittr)
require(lubridate)
require(zeallot)
require(parallel)

# Reading Data ------------------------------------------------------------

load("Data/data_I_wave.RData")

# Stan model compilation --------------------------------------------------

mc.cores <- parallel::detectCores()
# Stan options
rstan_options(auto_write = TRUE)

stan_Multip1 <- stan_model("Indep_SpT_Common_OOS.stan")

# Split -------------------------------------------------------------------
set.seed(130494)
dati %<>% group_by(denominazione_regione) %>% 
  mutate(
    Train=sample(c(rep(1, round(0.85*n())), rep(0, round(0.15*n()))))
  ) %>% 
  ungroup %>%
  mutate(denominazione_regione = droplevels(factor(denominazione_regione)))
datiTr <- dati %>% filter(Train==1)
datiTe <- dati %>% filter(Train==0)

# Data preparation for New Positives--------------------------------------------------------

Y <- as.integer(datiTr$NP)
N <- as.integer(length(Y))
NTe <- nrow(datiTe)
Nreg <- as.integer(length(unique(datiTr$denominazione_regione)))
timeIdx <- as.integer(datiTr$WW) + 1
regIdx <- as.integer(datiTr$denominazione_regione)
timeIdxTe <- as.integer(datiTe$WW) + 1
regIdxTe <- as.integer(datiTe$denominazione_regione)
t <- sort(unique(datiTr$WW)) + 1
Ntimes <- as.integer(length(t))

trScaled <- datiTr %>% dplyr::select(Swabs) %>% scale()
teScaled <- datiTe %>% dplyr::select(Swabs) %>% scale(center = attr(trScaled, "scaled:center"),
                                                      scale = attr(trScaled, "scaled:scale"))

X1 <- model.matrix(~.-1, data=as_tibble(trScaled))
X1Te <- model.matrix(~.-1, data=as_tibble(teScaled))
k1 <- ncol(X1)
lOff1 <- rep(logE, each = round(Ntimes*.85))
lOff1Te <- rep(logE, each = round(Ntimes*.15))

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
  "X" = X1,
  "NTest" = NTe, 
  "tIdTest" = timeIdxTe,
  "rIdTest" = regIdxTe,
  "lOffTest" = lOff1Te,
  "XTest" = X1Te)


# Stan Settings -----------------------------------------------------------

# Chains
n_chains <- 2
M <- 5000
n_cores <- mc.cores - 2

# Define a function to generate initial values

init <- function(chain_id = 1)
{
  list(# Richards
    logr = rnorm(1, log(ySums)-unique(lOff1), 1), 
    logh = rnorm(1, log(.5), .1),
    p = rnorm(1, tflex, 1), 
    logs = rnorm(1, log(1), .1),
    # Baseline
    logbase = rnorm(1, log(yMins+0.0000001)-unique(lOff1), 1))
} 


# Fit 
fit_Stan1 <- sampling(stan_Multip1, data = dat1, chains = n_chains, iter = M, 
                      cores = n_cores, init=init,
                      control = list(adapt_delta = 0.9, max_treedepth = 15)
)
