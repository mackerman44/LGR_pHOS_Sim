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
           
           # #------------------------------------------------------------------------------
           # ## Gather Detection Probabilities
           # #------------------------------------------------------------------------------
           # detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
           #                                    capHist_proc = proc_list$ProcCapHist) %>%
           #   mutate(spawn_yr = yr,
           #          species = spp,
           #          cv = sd/median) %>%
           #   select(spawn_yr, species, Node, n_tags, estimate = median, sd, cv, lowerCI, upperCI)
           
           #------------------------------------------------------------------------------
           ## Gather All Transition Probabilities for STADEM
           #------------------------------------------------------------------------------
           # abundance at main branches, upstream sites and black-boxes
           escape_post = as.matrix(stadem_mod$samples,
                                   iters = T,
                                   chains = T) %>%
             as_tibble() %>%
             select(CHAIN, ITER, matches("X.tot.new.hnc")) %>%
             tidyr::gather(param, value, -CHAIN, -ITER) %>%
             arrange(CHAIN, ITER) %>%
             sample_n(size = bootstrap_samp,
                      replace = T) %>%
             mutate(iter = 1:n()) %>%
             select(-CHAIN, -ITER) %>%
             select(iter, tot_escape = value) %>%
             # get transition probabilities
             left_join(compileTransProbs_LGD(dabom_mod,
                                             time_varying = F) %>%
                         select(-chain) %>%
                         group_by(param) %>%
                         mutate(iter = 1:n()) %>%
                         ungroup() %>%
                         arrange(param, iter) %>%
                         tidyr::spread(param, value) %>%
                         sample_n(size = bootstrap_samp,
                                  replace = T) %>%
                         mutate(iter = 1:n()) %>%
                         ungroup()) %>%
             mutate_at(vars(-iter, -tot_escape),
                       list(~ . * tot_escape)) %>%
             select(-tot_escape) %>%
             tidyr::gather(area, escape, -iter)
           
           # estimate the credible interval for each parameter
           credInt = escape_post %>%
             spread(area, escape) %>%
             select(-iter) %>%
             coda::as.mcmc() %>%
             coda::HPDinterval(prob = 0.95) %>%
             as_tibble(rownames = 'area') %>%
             rename(lowerCI = lower,
                    upperCI = upper)
           
           escape_summ = escape_post %>%
             group_by(area) %>%
             summarise(mean = mean(escape),
                       median = median(escape),
                       mode = estMode(escape),
                       sd = sd(escape),
                       cv = sd / mean) %>%
             mutate_at(vars(mean, median, mode, sd),
                       list(~ ifelse(. < 0, 0, .))) %>%
             left_join(credInt,
                       by = 'area') %>%
             mutate_at(vars(mean, median, mode),
                       list(round)) %>%
             mutate_at(vars(sd, lowerCI, upperCI),
                       list(round),
                       digits = 1) %>%
             mutate_at(vars(cv),
                       list(round),
                       digits = 3)
           
           escape_summ %<>%
             filter(area != "deviance")
           
           return(escape_summ)
           
         }) %>%
  mutate(Species = spp,
         Type = 'HNC') %>%
  select(Year, Species, Type, everything())

sthd_hnc_summ %>%
  filter(!is.na(cv)) %>%
  write_csv('hnc_dabom/ModelFits/Sthd_HNC_abund.csv')
