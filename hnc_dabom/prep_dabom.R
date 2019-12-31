# Author: Kevin See
# Purpose: Process capture histories for hatchery no-clip steelhead from several years at Lower Granite, and run DABOM on the observations
# Created: 12/31/19
# Last Modified: 12/31/19
# Notes: 

#-------------------------------
# load needed packages
library(tidyverse)
library(magrittr)
library(lubridate)
library(STADEM)
library(WriteXLS)
library(PITcleanr)
library(janitor)

#-------------------------------
# use STADEM results from SnakeBasinFishStatus Git repo
# set up folder structure
stademFolder = '/Users/seek/Documents/GitProjects/MyProjects/SnakeBasinFishStatus/STADEM_results'

# pull out total HNC abundance for each year
sthd_hnc = 2010:2019 %>%
  as.list() %>%
  purrr::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           
           load(paste0(stademFolder,'/LGR_STADEM_Steelhead_', x, '.rda'))
           
           res = stadem_mod$summary %>%
             as_tibble(rownames = 'param') %>%
             filter(param == 'X.tot.new.hnc')
           rm(stadem_mod, stadem_list)
           return(res)
         })

# save results
write_rds(sthd_hnc,
          path = 'data/DABOM/Sthd_HNC_STADEM.rds')

#-------------------------------
# build detection site configuration 
# build configuration table (requires internet connection)
org_config = buildConfig()

# customize some nodes based on DABOM framework
configuration = org_config %>%
  mutate(Node = ifelse(SiteID %in% c('VC2', 'VC1', 'LTR', 'MTR', 'UTR'),
                       SiteID,
                       Node),
         Node = ifelse(SiteID == 'SC2',
                       'SC2B0',
                       Node),
         Node = ifelse(SiteID %in% c('CROTRP',
                                     'CRT',
                                     'REDTRP',
                                     'REDR',
                                     'RRT'),
                       'SC2A0',
                       Node),
         Node = ifelse(Node == 'ACB',
                       'ACBB0',
                       Node),
         Node = ifelse(Node == 'CCA',
                       'CCAB0',
                       Node),
         Node = ifelse(SiteID == 'AFC',
                       ifelse(grepl('MAINSTEM', AntennaGroup),
                              'AFCB0',
                              'AFCA0'),
                       Node),
         Node = ifelse(SiteID == 'HBC',
                       'HYCA0',
                       Node),
         Node = ifelse(SiteID %in% c('TUCH', 'TFH'),
                       'TUCH',
                       Node),
         Node = ifelse(SiteID == 'MCCA',
                       'STR',
                       Node),
         Node = ifelse(SiteID == 'CARMEC',
                       'CRCA0',
                       Node),
         Node = ifelse(SiteID == 'BIG2C',
                       'TAYA0',
                       Node),
         Node = ifelse(SiteID == 'WIMPYC',
                       'WPCA0',
                       Node),
         Node = ifelse(SiteID == 'IML' & ConfigID == 130 & AntennaID == '09',
                       'IMLA0',
                       Node),
         Node = str_replace(Node, '^BTC', 'BTL'),
         Node = ifelse(SiteID %in% c('YANKFK', 'CEY'),
                       'YFKA0',
                       Node),
         Node = ifelse(SiteID == 'SAWT',
                       'STL',
                       Node),
         Node = ifelse(SiteID == 'LOOH',
                       'LOOKGC',
                       Node),
         Node = ifelse(SiteID == 'RPDTRP',
                       'RAPH',
                       Node),
         Node = ifelse(SiteID == 'CHARLC',
                       'CCAB0',
                       Node),
         Node = ifelse(Node == 'KEN',
                       'KENB0',
                       Node),
         Node = ifelse(Node == 'HYC',
                       'HYCB0',
                       Node),
         Node = ifelse(Node == 'YFK',
                       'YFKB0',
                       Node),
         Node = ifelse(Node == 'LLR',
                       'LLRB0',
                       Node),
         Node = ifelse(Node == 'LRW',
                       'LRWB0',
                       Node),
         Node = ifelse(SiteID == '18M',
                       str_replace(Node, '18M', 'HEC'),
                       Node)) %>%
  distinct()


# Node network for DABOM
site_df = writeLGRNodeNetwork()

# remove some sites that have been combined with others (see the modifications to the configuration file)
site_df = site_df %>%
  filter(!SiteID %in% c('TFH',
                        'MCCA',
                        'WIMPYC',
                        'YANKFK', 'CEY',
                        'SAWT',
                        'LOOH',
                        'CARMEC',
                        'BIG2C',
                        'RPDTRP'))

# Save file.
save(configuration, site_df, file = 'data/DABOM/site_config.rda')

#-------------------------------
# file paths to some data (from SnakeBasinFishStatus Git repo)
trap_path = "/Users/seek/Documents/GitProjects/MyProjects/SnakeBasinFishStatus/data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv"
# trap_path = "/Users/seek/OneDrive - Merck Sharp & Dohme, Corp/Data/LGR_trap_database/20191231/tblLGDMasterCombineExportJodyW.csv"
ptagis_obs = "/Users/seek/Documents/GitProjects/MyProjects/SnakeBasinFishStatus/data/CompleteTagHistories/"

# read in trap database
trap_df = read_csv(trap_path)

# all the detection histories for HNC fish were already downloaded from PTAGIS, because they were initially marked "valid" in the trap database. So we just need to re-process them with PITcleanr
all_tags = 2010:2019 %>%
  as.list() %>%
  purrr::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           read_csv(paste0(ptagis_obs, "LGR_Steelhead_", x, ".csv")) %>%
                      select(LGDNumPIT = `Tag Code`) %>%
                      distinct() %>%
                      left_join(trap_df %>%
                                  filter(SpawnYear == paste0('SY', x)) %>%
                                  select(LGDNumPIT, SRR, GenRear, LGDMarkAD, LGDValid)) %>%
                      distinct()
         }) %>%
  filter(grepl('^3', SRR)) %>%
  mutate_at(vars(SRR, GenRear, LGDMarkAD),
            list(as.factor))

# pull out the hatchery no-clip tags
hnc_tags = all_tags %>%
  filter(SRR == '32H',
         LGDMarkAD == "AI")


spp = 'Steelhead'

for(yr in 2017:2019) {
  
  cat(paste('Starting year', yr, '\n'))
  
  # start date is July 1 of the previous year for steelhead, or March 1 of current year for Chinook
  startDate = if_else(spp == 'Steelhead',
                      paste0(yr-1, '0701'),
                      if_else(spp == 'Chinook',
                              paste0(yr, '0301'),
                              NULL))
  
  # build parent-child table
  parent_child = createParentChildDf(site_df,
                                     configuration,
                                     startDate = startDate)
  
  # get raw observations from PTAGIS
  # These come from running a saved query on the list of tags to be used
  observations = read_csv(paste0(ptagis_obs, '/LGR_', spp, '_', yr, '.csv')) %>%
    inner_join(hnc_tags %>%
                 filter(Year == yr) %>%
                 select(`Tag Code` = LGDNumPIT))
  
  proc_list = processCapHist_LGD(species = spp,
                                 spawnYear = yr,
                                 configuration = configuration,
                                 trap_path = trap_path,
                                 filter_by_PBT = F,
                                 parent_child = parent_child,
                                 observations = observations,
                                 site_df = site_df,
                                 truncate = T,
                                 step_num = 3,
                                 save_file = T,
                                 file_name = paste0("data/DABOM/Process_CH_LGR_", spp, '_', yr, '.xlsx'))
  
  save(startDate, site_df, configuration, parent_child, proc_list,
       file = paste0("data/DABOM/DABOM_preppd_LGR_", spp, '_', yr, '.rda'))
  
}

