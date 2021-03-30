
# Packages ----------------------------------------------------------------

require(rstan)
require(tidyverse)
require(magrittr)
require(lubridate)
require(zeallot)
require(parallel)
require(mvtnorm)
source("OtherFuns/DrichFuns.R")
source("OtherFuns/DataFuns.R")
source("OtherFuns/StanFuns.R")
set.seed(130494)

# Reading Data ------------------------------------------------------------

load("Data/DatiStan18.RData")
dati <- dati18
rm(dati18)
adj <- adjTr
rm(adjTr)


# Stan model compilation --------------------------------------------------

mc.cores = parallel::detectCores()
# Stan options
rstan_options(auto_write = TRUE)

stan_CARSpMultip <- stan_model("CovCARSp_MultipAll.stan")


# Data preparation for New Positives--------------------------------------------------------

Y <- as.integer(dati$NP)
N <- as.integer(length(Y))
Nreg <- as.integer(length(unique(dati$denominazione_regione)))
timeIdx <- as.integer(dati$WW)+1
regIdx <- as.integer(droplevels(dati$denominazione_regione))
t <- unique(dati$WW)
Ntimes <- as.integer(length(t))
X <- model.matrix(~., data=dati %>% dplyr::select(Swabs))
k <- ncol(X)
lOff <- rep(0, N)

# Useful data summaries
#ySums <-  dati %>% group_by(denominazione_regione) %>% summarise(y=sum(NP)) %$% y
tflex <- dati %>% group_by(denominazione_regione) %>% summarise(WW=WW[which.max(NP)]) %$% WW

# Prepare data
datSp <- list(
  "N" = N, 
  "W" = adj,
  "W_n" = sum(adj[upper.tri(adj)]>0),
  "Y" = Y,
  "Ntimes" = Ntimes,
  "t" = t,
  "tId" = timeIdx,
  "Nreg" = Nreg, 
  "rId" = regIdx,
  "lOff" = lOff,
  "k" = k,
  "X" = X)

# Stan Settings -----------------------------------------------------------

# Chains
n_chains <- 2
M <- 10000
n_cores <- (mc.cores)

# Define a function to generate initial values
load("Data/Inits.RData")
ppp %<>% filter(Regione%in%dati$denominazione_regione)

init <- function(chain_id = 1)
{
  list(# Richards
    logr = ppp$logr, 
    logh = ppp$logh,
    p = ppp$p, 
    logs = log(ppp$s),
    # Baseline
    base = exp(ppp$beta0)
)
}
fit_StanSp <- sampling(stan_CARSpMultip, data = datSp, chains = n_chains, iter = M, 
                     cores = n_cores, init=init,
                     control = list(adapt_delta = 0.9, max_treedepth = 15)
)
c(postsamples_StanSp, ypreds_StanSp, ypredsQ_StanSp) %<-% extract_postY(fit_StanSp)
# Diagnostic check: yellow or red points are bad!

# Saving the WS ----------------------------------------------------------

save.image(file="WS/SwabsCovCARSp.RData")
