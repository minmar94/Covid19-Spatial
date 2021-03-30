
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

load("Data/DatiStan19.RData")
dati <- dati19
rm(dati19)
adj <- adjGeo
rm(adj19)
logE <- logE19
rm(logE19)

# Stan model compilation --------------------------------------------------

mc.cores = parallel::detectCores()
# Stan options
rstan_options(auto_write = TRUE)

stan_Multip1 <- stan_model("CovStCARSp_MultipAll.stan")

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

# Useful data summaries
ySums <-  dati %>% group_by(denominazione_regione) %>% summarise(y=sum(NP)) %$% y
yMins <-  dati %>% group_by(denominazione_regione) %>% summarise(y=min(NP)) %$% y
tflex <- dati %>% group_by(denominazione_regione) %>% summarise(WW=WW[which.max(NP)]) %$% WW

# Prepare data
dat1 <- list(
  "N" = N, 
  "W" = adj,
  "W_n" = sum(adj[upper.tri(adj)]>0),
  "Y" = Y,
  "Ntimes" = Ntimes,
  "t" = t,
  "tId" = timeIdx,
  "Nreg" = Nreg, 
  "rId" = regIdx,
  "lOff" = lOff1,
  "k" = k1,
  "X" = X1)

# Stan Settings -----------------------------------------------------------

# Chains
n_chains <- 2
M <- 10000
n_cores <- (mc.cores)

# Define a function to generate initial values

init <- function(chain_id = 1)
{
  list(# Richards
    logr = rnorm(Nreg, log(ySums)-logE, .5), 
    logh = rnorm(Nreg, log(.5), .1),
    p = rnorm(Nreg, tflex, 1), 
    logs = rnorm(Nreg, log(1), .1),
    # Baseline
    logbase = rnorm(Nreg, log(yMins+0.0000001)-logE, .1))
} 


fit_Stan1 <- sampling(stan_Multip1, data = dat1, chains = n_chains, iter = M, 
                     cores = n_cores, init=init,
                     control = list(adapt_delta = 0.9, max_treedepth = 12)
)

c(postsamples_Stan1, ypreds_Stan1, ypredsQ_Stan1) %<-% extract_postY(fit_Stan1)

save.image(file="WS/StCARMultipAll_TrueData.RData")
