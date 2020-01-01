# Author: Kevin See
# Purpose: Run DABOM for hatchery no-clip steelhead from several years at Lower Granite
# Created: 12/31/19
# Last Modified: 12/31/19
# Notes: 

#-------------------------------
# load needed packages
library(tidyverse)
library(magrittr)
library(lubridate)
# library(STADEM)
# library(WriteXLS)
# library(PITcleanr)
library(DABOM)
library(jagsUI)
# library(janitor)

#-------------------------------
# should initial movement probabilities be time-varying?
time_varying = TRUE

# file path to the default and initial model 
basic_modNm = 'hnc_dabom/ModelFiles/LGR_DABOM.txt'

writeDABOM_LGD(file_name = basic_modNm,
               time_varying = time_varying)

#-------------------------------
spp = 'Steelhead'
yr = 2019

for(yr in 2017:2019) {
  
  load(paste0("data/DABOM/DABOM_preppd_LGR_", spp, '_', yr, '.rda'))
  
  proc_ch <- proc_list$ProcCapHist %>%
    filter(AutoProcStatus)
  
  #------------------------------------------------------------------------------
  # Because HLM was giving us problems...
  # Switch Potlatch detections - move POTREF and POTRWF to HLMA0 with det = 1.0
  #------------------------------------------------------------------------------
  proc_ch %<>%
    mutate(Node = ifelse(Node %in% c('POTREF', 'POTRWF'), 'HLMA0', Node))
  
  #------------------------------------------------------------------------------
  # In 2010 and 2011, the Asotin Creek weir was upstream of ACB. For those years, 
  # take all detections at ASOTIC and move them to ACBA0
  #------------------------------------------------------------------------------
  if(spp == 'Steelhead' & yr %in% c(2010, 2011)) {
    proc_ch <- proc_ch %>%
      mutate(Node = if_else(Node == 'ASOTIC',
                            'ACBA0',
                            Node))
  }
  
  # filepath for specific JAGS model code for species and year
  mod_path = paste0('hnc_dabom/ModelFiles/LGR_DABOM_HNC_', spp, '_', yr, '.txt')
  
  #writes species and year specific jags code
  fixNoFishNodes(basic_modNm,
                 mod_path,
                 proc_ch,
                 proc_list$NodeOrder)
  
  #------------------------------------------------------------------------------
  # Create capture history matrices for each main branch to be used in 
  # the JAGS data list
  #------------------------------------------------------------------------------
  dabom_list = createDABOMcapHist(proc_ch,
                                  proc_list$NodeOrder,
                                  split_matrices = T)
  #------------------------------------------------------------------------------
  # Only Used to Debug
  #------------------------------------------------------------------------------
  # full_dabom = createDABOMcapHist(proc_ch,
  #                                proc_list$NodeOrder,
  #                                split_matrices = F)
  #------------------------------------------------------------------------------
  
  # Creates a function to spit out initial values for MCMC chains
  init_fnc = setInitialValues_LGD(dabom_list)
  
  #Create all the input data for the JAGS model
  jags_data = createJAGSinputs_LGD(dabom_list)
  
  
  # CHANGE TIME VARY DATE FOR SPECIES!!!!!!
  
  if(time_varying) {
    if(spp == 'Steelhead'){
      jags_data = c(jags_data,
                    addTimeVaryData(proc_ch,
                                    node_order = proc_list$NodeOrder,
                                    start_date = paste0(yr-1,'0701'), 
                                    end_date = paste0(yr,'0630')))
    } else {
      jags_data = c(jags_data,
                    addTimeVaryData(proc_ch,
                                    node_order = proc_list$NodeOrder,
                                    start_date = paste0(yr-1,'0301'), 
                                    end_date = paste0(yr,'0817')))
    }
  }
  
  #------------------------------------------------------------------------------
  # Tell JAGS which parameters in the model that it should save.
  # the fnc is hard coded and needs to be updated if there are changes!
  #------------------------------------------------------------------------------
  
  jags_params = setSavedParams_LGD(time_varying = time_varying)
  
  #------------------------------------------------------------------------------
  # Run the model
  
  # Recommended MCMC parameters are:
  #   
  #   * `n.chains`: 4
  # * `n.iter`: 5,000
  # * `n.burnin`: 2,500
  # * `n.thin`: 10
  # 4*(5000+2500) = 30000
  # 1 iteration takes about .18 minutes
  # about 3.75 days!!!!!
  
  set.seed(12)
  dabom_mod <- jags.basic(data = jags_data,
                          inits = init_fnc,
                          parameters.to.save = jags_params,
                          model.file = mod_path,
                          n.chains = 2,
                          n.iter = 4,
                          n.burnin = 2,
                          DIC = F)
                          # n.chains = 4,
                          # n.iter = 5000,
                          # n.burnin = 2500,
                          # n.thin = 10,
                          # DIC = T)
  
  
  
  #--------------------------------------------------------------------------------
  # Save the results
  #--------------------------------------------------------------------------------
  save(dabom_mod, dabom_list, proc_list,
       file = paste0('hnc_dabom/ModelFits/LGR_DABOM_HNC_', spp, '_', yr,'.rda'))
  
}
