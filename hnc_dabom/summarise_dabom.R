# Author: Kevin See
# Purpose: sumarise DABOM results for hatchery no-clip steelhead from several years at Lower Granite
# Created: 1/7/2020
# Last Modified: 1/7/2020
# Notes: 

#-------------------------------
# load needed packages
library(tidyverse)
library(magrittr)
library(lubridate)
library(STADEM)
library(DABOM)
library(jagsUI)

#-------------------------------
spp = 'Steelhead'

# where are STADEM results stored?
stademFolder = '/Users/seek/Documents/GitProjects/MyProjects/SnakeBasinFishStatus/STADEM_results'

bootstrap_samp = 2000

sthd_hnc_summ = 2017:2019 %>%
  as.list() %>%
  purrr::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           
           # load STADEM JAGS model
           load(paste0(stademFolder,'/LGR_STADEM_', spp, '_', x, '.rda'))
           
           # load DABOM JAGS model
           load(paste0('hnc_dabom/ModelFits/LGR_DABOM_HNC_', spp, '_', x,'.rda'))
           
           escape_summ = calcTribEscape_LGD(dabom_mod, 
                                            stadem_mod,
                                            stadem_param_nm = 'X.tot.new.hnc',
                                            node_order = proc_list$NodeOrder,
                                            time_varying = F)
           
           
           return(escape_summ)
           
         }) %>%
  mutate(Species = spp,
         Type = 'HNC') %>%
  select(Year, Species, Type, everything())

sthd_hnc_summ %>%
  filter(!is.na(cv)) %>%
  write_csv('hnc_dabom/ModelFits/Sthd_HNC_abund.csv')
