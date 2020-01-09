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
library(readxl)


# source some needed functions
source('R/definePopulations.R')

#-------------------------------
spp = 'Steelhead'

# where are STADEM results stored?
stademFolder = '/Users/seek/Documents/GitProjects/MyProjects/SnakeBasinFishStatus/STADEM_results'

bootstrap_samp = 2000

set.seed(5)
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

# group things by population
pop_df = definePopulations(spp)

sthd_hnc_trt = sthd_hnc_summ %>%
  left_join(pop_df) %>%
  filter(!is.na(TRT)) %>%
  mutate_at(vars(TRT),
            list(fct_explicit_na)) %>%
  mutate_at(vars(Year),
            list(as.numeric)) %>%
  group_by(Year, Species, TRT) %>%
  summarise(est_hnc = sum(mean),
            sd_hnc = sqrt(sum(sd^2)),
            cv_hnc = sd_hnc / est_hnc) %>%
  filter(est_hnc > 0)

sthd_wild_trt = read_excel("../SnakeBasinFishStatus/Abundance_results/LGR_AllSummaries_Steelhead.xlsx",
                           sheet = "Pop Total Esc") %>%
  select(Year = spawn_yr,
         Species = species,
         TRT,
         est_wild = mean,
         sd_wild = sd, 
         cv_wild = cv)

sthd_hnc_trt %>%
  inner_join(sthd_wild_trt) %>%
  rowwise() %>%
  mutate(pHOS = est_hnc / (est_wild + est_hnc),
         pHOS_se = msm::deltamethod(~ x1 / (x1 + x2),
                                    mean = c(est_hnc, est_wild),
                                    cov = diag(c(sd_hnc, sd_wild)^2)),
         pHOS_cv = pHOS_se / pHOS) %>%
  ungroup() %>%
  select(Year:TRT, 
         HNC = est_hnc,
         Wild = est_wild,
         starts_with("pHOS")) %>%
  write_csv('hnc_dabom/ModelFits/Sthd_HNC_pHOS.csv')
  