require(tidyverse)
require(magrittr)
require(lubridate)


# Functions ---------------------------------------------------------------
# Read data fun
read_italian <- function(path){
  
  dati_Ita <- data.table::fread(path) %>% as_tibble() %>% 
    mutate(data = date(data)) %>% 
    dplyr::select(data:casi_testati, -starts_with("note"), -starts_with("casi_da"), -starts_with("variazione")) %>% 
    mutate_if(is.numeric, abs)
  
  colnames(dati_Ita)[-c(1,2)] <- c("Hospitalized with symptoms", "Intensive care", "Current hospitalized",
                                   "Home isolation", "Current positives", "New positives", 
                                   "Discharged recovered", "Deceased",  "Cumulative positives", "Swabs", "Tested cases")
  dati_Ita$`Tested cases`[dati_Ita$data == "2020-12-03"] <- 13264937
  dati_Ita$`Tested cases`[dati_Ita$data == "2020-12-05"] <- 13435624
  return(dati_Ita)
  
}

# Read data fun
read_regional <- function(path){
  
  dati_reg <- data.table::fread(path) %>% as_tibble() %>% 
    dplyr::select(-note, -variazione_totale_positivi, -casi_da_sospetto_diagnostico, -casi_da_screening) %>% 
    mutate_if(is.numeric, abs) %>% 
    # 
    mutate(denominazione_regione = ifelse(denominazione_regione %in% c("P.A. Bolzano", "P.A. Trento"), "Trentino-Alto Adige",
                                          denominazione_regione),
           data = date(data)) %>% 
    gather(Key, Value, ricoverati_con_sintomi:casi_testati) %>% 
    group_by(Key, data, denominazione_regione) %>% 
    dplyr::summarise(Value = sum(Value)) %>% 
    ungroup() 
  
  dati_reg$Key <- factor(dati_reg$Key, levels = unique(dati_reg$Key), 
                         labels = c("Tested cases","Deceased", "Discharged recovered", "Home isolation", "New positives", "Hospitalized with symptoms", 
                                    "Swabs", "Intensive care", "Cumulative positives", "Current hospitalized", "Current positives"))
  
  dati_reg$Value[dati_reg$Key == "Casi testati" & dati_reg$denominazione_regione == "Campania" & dati_reg$data == "2020-12-03"] <- 1118787
  dati_reg$Value[dati_reg$Key == "Casi testati" & dati_reg$denominazione_regione == "Molise" & dati_reg$data == "2020-12-05"] <- 89725 
  return(dati_reg)
}

read_residents <- function(path){
  
  residents <- data.table::fread(input = path, sep = ",") %>% as_tibble()
  
  # Fix Labels
  # Regions
  residents %<>% 
    mutate(
      Territorio = ifelse(Territorio == "Valle Aosta", "Valle d'Aosta", Territorio),
      Territorio = ifelse(Territorio == "Emilia Romagna", "Emilia-Romagna", Territorio),
      Territorio = ifelse(Territorio == "Friuli V. G.", "Friuli Venezia Giulia", Territorio),
      Territorio = ifelse(Territorio == "TrentinoAltoAdige", "Trentino-Alto Adige", Territorio),
    ) 
  # Provinces
  to_change <- c("AscoliPiceno", "ProvinciaAutonomaTrento", "MonzaedellaBrianza", "ReggiodiCalabria",
                 "ProvinciaAutonomaBolzano/Bozen", "ViboValentia", "LaSpezia", "ReggionellEmilia", "Massa-Carrara")
  
  residents$Territorio[residents$Territorio %in% to_change] <- c("La Spezia", "Monza e della Brianza", "Bolzano", "Trento", "Reggio nell'Emilia",
                                                                 "Massa Carrara", "Ascoli Piceno", "Reggio di Calabria", "Vibo Valentia")
  
  return(residents)
}

