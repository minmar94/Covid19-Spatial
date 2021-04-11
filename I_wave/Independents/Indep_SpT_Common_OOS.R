
# Packages ----------------------------------------------------------------

require(rstan)
require(tidyverse)
require(magrittr)
require(lubridate)
require(zeallot)
require(parallel)
source("OtherFuns/DrichFuns.R")
source("OtherFuns/DataFuns.R")
source("OtherFuns/StanFuns.R")
set.seed(130494)

# Reading Data ------------------------------------------------------------

load("Data/DatiStanIIwave.RData")
dati <- datasave
dati$denominazione_regione <- factor(dati$denominazione_regione)
rm(datasave)

# Stan model compilation --------------------------------------------------

mc.cores = parallel::detectCores()
# Stan options
rstan_options(auto_write = TRUE)

stan_Multip1 <- stan_model("IndTRnd_MultipAll_Common.stan")

# Data preparation for New Positives--------------------------------------------------------

Y <- as.integer(dati$NP)
N <- as.integer(length(Y))
Nreg <- as.integer(length(unique(dati$denominazione_regione)))
timeIdx <- as.integer(dati$WW-min(dati$WW))+1
regIdx <- as.integer(droplevels(dati$denominazione_regione))
t <- unique(dati$WW-min(dati$WW)) + 1
Ntimes <- as.integer(length(t))

X1 <- model.matrix(~.-1, data=scale(dati %>% dplyr::select(NewSwabsSett)) %>% as.data.frame)
k1 <- ncol(X1)
lOff1 <- log(dati$totale/10000)

# Useful data summaries
ySums <-  dati %>% group_by(denominazione_regione) %>% summarise(y=sum(NP)) %$% y
yMins <-  dati %>% group_by(denominazione_regione) %>% summarise(y=min(NP)) %$% y
tflex <- dati %>% group_by(denominazione_regione) %>% summarise(WW=(WW-min(WW))[which.max(NP)]) %$% WW

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


fit_Stan1 <- sampling(stan_Multip1, data = dat1, chains = n_chains, iter = M, 
                      cores = n_cores, init=init,
                      control = list(adapt_delta = 0.9, max_treedepth = 15)
)

c(postsamples_Stan1, ypreds_Stan1, ypredsQ_Stan1, ypredsT_Stan1, ypredsTQ_Stan1) %<-% extract_postY(fit_Stan1)

save.image()
