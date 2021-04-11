
# Packages ----------------------------------------------------------------

require(tidyverse)
require(magrittr)
require(lubridate)
require(zeallot)
require(mvtnorm)
require(ggforce)
source("OtherFuns/DataFuns.R")


# Reading Data ------------------------------------------------------------

# Reading data
dati_reg <- read_regional("https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv") %>% spread(Key, Value)

# First wave
dati_regI <- dati_reg %>% filter(data <= "2020-06-28") 
dati <- dati_regI %>% 
  mutate(denominazione_regione = factor(denominazione_regione, levels = unique(denominazione_regione))
         #,tnum = as.numeric(data) - min(as.numeric(data))
  ) %>% 
  group_by(denominazione_regione, WW = week(data), Y = year(data)) %>% 
  summarise(data = first(data), NP = sum(`New positives`)
            #, tnum = first(tnum)
  ) %>% 
  #ggplot(aes(x = data, y = NP, color = denominazione_regione)) + geom_line()
  ungroup() %>% 
  mutate(WW = WW - min(WW)) %>% 
  arrange(denominazione_regione, data) %>% filter(!denominazione_regione%in%c("Molise", "Valle d'Aosta"))

Resid <- read_residents("Data/residenti2019.csv") %>% filter(Territorio %in% unique(dati$denominazione_regione)) %>% 
  arrange(Territorio)
logE <- log(Resid$totale/10000)

# Graph
GraphData <- readxl::read_xlsx("Data/mobility_data.xlsx", sheet = "total",
                               skip=1) %>% dplyr::select(-1) %>% 
  mutate_at(vars(-Regions), .funs = function(x) as.numeric(x>0)) %>% 
  gather(region2, value, -Regions) %>% arrange(Regions, region2) %>% 
  filter(!Regions%in%c("Molise", "Valle d'Aosta")) %>%  filter(!region2%in%c("Molise", "Valle d'Aosta")) %>% 
  spread(region2, value)

adj <- GraphData %>% dplyr::select(-1) %>% as.matrix()
colnames(adj) <- NULL
isSymmetric(adj)

save(dati, logE, adj, file = "Data/DatiStan.RData")