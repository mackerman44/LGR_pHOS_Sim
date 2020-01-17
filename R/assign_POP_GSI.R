# Purpose: Helper function to assign MPG, TRTs and GSI groups
# to Sites and Nodes.
# Author: Ryan Kinzer
# Created: 7/12/19
# Modified:

assign_POP_GSI <- function(species, configuration, site_df){

  load('data/Snake_POP_metadata.rda')

  config <- configuration %>%
    filter(SiteID %in% site_df$SiteID) %>%
    filter(SiteID != 'GRA') %>%
    select(SiteID, SiteType, SiteName, RKM, Latitude, Longitude, Node) %>%
    distinct() %>%
    sf::st_as_sf(coords = c('Longitude', 'Latitude'),
                 crs = 4326)

  if(species == 'Steelhead'){
    pop_sf <- SR_st_pop
    
    ignore <- c('USE', 'USI', 'SFG')
    sfclw <-  c('SC1', 'SC2')
    
    df <- config %>%
      sf::st_join(pop_sf %>%
                    select(ESU_DPS, MPG, POP_NAME, TRT = TRT_POPID, GSI_Group)) %>%
      mutate(POP_NAME = ifelse(SiteID %in% ignore, NA, POP_NAME),
             TRT = ifelse(SiteID %in% ignore, NA, TRT)) %>%
      mutate(POP_NAME = ifelse(SiteID %in% sfclw, 'South Fork Clearwater River', POP_NAME),
             TRT = ifelse(SiteID %in% sfclw, 'CRSFC-s', TRT))

      }
  
  if(species == 'Chinook'){
    pop_sf <- SR_ch_pop
    
    ignore <- c('USE', 'USI', 'SFG', 'IR1', 'IR2','WR1', 'UGR')
    sfclw <-  c('SC1', 'SC2')
    
    df <- config %>%
      sf::st_join(pop_sf %>%
                    select(ESU_DPS, MPG, POP_NAME, TRT = TRT_POPID, GSI_Group)) %>%
      mutate(POP_NAME = ifelse(SiteID %in% ignore, NA, POP_NAME),
             TRT = ifelse(SiteID %in% ignore, NA, TRT)) %>%
      mutate(POP_NAME = ifelse(SiteID %in% sfclw, 'Upper South Fork Clearwater', POP_NAME),
             TRT = ifelse(SiteID %in% sfclw, 'SCUMA', TRT))
    
  }
  

    #as_tibble() %>%
    #select(-geometry)
    #   mutate_at(vars(obsGSI, TRT, modBranch, MPG, NWR_POPID),
    #             list(~ if_else(AssignSpawnNode == 'GRA',
    #                            as.character(NA),
    #                            .)))
  
  return(list(df, pop_sf))

}
