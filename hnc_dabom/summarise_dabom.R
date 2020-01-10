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
library(PITcleanr)
library(jagsUI)
library(readxl)
library(sf)

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
  
#-------------------------------
# produce summaries of all HNC tags
#-------------------------------
# Ryan Kinzer put this data together in the SnakeBasinFishStatus repo
load('data/Snake_POP_metadata.rda')
pop_sf <- SR_st_pop

ignore <- c('USE', 'USI', 'SFG')
sfclw <-  c('SC1', 'SC2')


sthd_hnc_tag_summ = 2017:2019 %>%
  as.list() %>%
  purrr::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           
           # load some configuration data
           load(paste0("data/DABOM/DABOM_preppd_LGR_", spp, '_', x, '.rda'))
           
           config <- configuration %>%
             filter(SiteID %in% site_df$SiteID) %>%
             filter(SiteID != 'GRA') %>%
             select(SiteID, SiteType, SiteName, RKM, Latitude, Longitude, Node) %>%
             distinct() %>%
             sf::st_as_sf(coords = c('Longitude', 'Latitude'),
                          crs = 4326) %>%
             st_transform(st_crs(pop_sf))
           
           site_pop = config %>%
             sf::st_join(pop_sf %>%
                           select(ESU_DPS, MPG, POP_NAME, TRT = TRT_POPID, GSI_Group)) %>%
             mutate(POP_NAME = ifelse(SiteID %in% ignore, NA, POP_NAME),
                    TRT = ifelse(SiteID %in% ignore, NA, TRT)) %>%
             mutate(POP_NAME = ifelse(SiteID %in% sfclw, 'South Fork Clearwater River', POP_NAME),
                    TRT = ifelse(SiteID %in% sfclw, 'CRSFC-s', TRT))
           
           rm(proc_list)
           
           # load DABOM JAGS model
           load(paste0('hnc_dabom/ModelFits/LGR_DABOM_HNC_', spp, '_', x,'.rda'))
           
           tag_summ = summariseTagData(capHist_proc = proc_list$ProcCapHist %>%
                                         mutate(UserProcStatus = AutoProcStatus),
                                       trap_data = proc_list$ValidTrapData) %>%
             left_join(site_pop %>%
                         as_tibble() %>%
                         select(Node, MPG, POP_NAME, TRT),
                       by = c('AssignSpawnNode' = 'Node'))
         })

# some tags appear more than once in the LGR trap database. Try to only keep one record per tag & year
sthd_hnc_tag_summ %<>%
  group_by(TagID, Year) %>%
  filter(CollectionDate == max(CollectionDate)) %>%
  slice(1) %>%
  ungroup()
  

# save to share with IDFG
sthd_hnc_tag_summ %>%
  rename(DABOM_branch = Group) %>%
  select(-BranchNum) %>%
  write_csv("outgoing/Sthd_HNC_SpawnSites.csv")


# for fish detected somewhere, which hatchery did they come from, and what model branch were they detected in?
sthd_hnc_tag_summ %>%
  filter(!is.na(Group)) %>%
  mutate_at(vars(Group, GenParentHatchery, TRT),
            list(fct_explicit_na)) %>%
  rename(Branch = Group) %>%
  janitor::tabyl(Branch, GenParentHatchery,
                 show_missing_levels = F) %>%
  janitor::adorn_totals(where = c("row", "col")) %>%
  janitor::adorn_percentages(denominator = "col") %>%
  janitor::adorn_pct_formatting() %>%
  janitor::adorn_ns()
  

# these are all the fish that fell into the main black box, by hatchery.
# Maybe they were recovered by the hatchery, so no 
sthd_hnc_tag_summ %>%
  mutate(Main_bb = if_else(!is.na(Group), T, F)) %>%
  mutate_at(vars(Group, GenParentHatchery, TRT),
            list(fct_explicit_na)) %>%
  janitor::tabyl(GenParentHatchery, Main_bb) %>%
  # janitor::adorn_totals(where = c("row", "col")) %>%
  janitor::adorn_percentages(denominator = "row") %>%
  janitor::adorn_pct_formatting() %>%
  janitor::adorn_ns()


sthd_hnc_tag_summ %>%
  filter(!is.na(Group)) %>%
  # select(AssignSpawnSite, PtagisEventLastSpawnSite) %>%
  filter(AssignSpawnSite != PtagisEventLastSpawnSite) %>%
  select(Year,
         DABOM = AssignSpawnSite, Ptagis = PtagisEventLastSpawnSite, 
         TagPath, PtagisEventSites)
