#' @title Define Popuation Groups - LGD
#'
#' @description Define which main branches or detection sites fall into TRT populations.
#'
#' @author Ryan Kinzer
#'
#' @param spp Species, either "Chinook" or "Steelhead"
#' @param node_order output of function \code{createNodeOrder}.
#'
#' @import dplyr stringr
#' @export
#' @return NULL
#' @examples definePopulations()

definePopulations = function(spp = c('Chinook', 'Steelhead')) {
  
  spp = match.arg(spp)
  
  if(spp == 'Chinook') {
    # NEEDS FIXED
    report_df = list('CRLAP' = 'LAP',
                     'CRLOC' = 'LRL',
                     'CRLOL' = 'LC1',
                     'CRPOT' = 'JUL',
                     'GRCAT' = 'CATHEW',
                     'GRLOO' = 'LOOKGC',
                     'GRLOS' = 'WR2',
                     'GRLOS/GRMIN' = c('WR1_bb', 'WR2'),
                     'GRMIN' = 'MR1',
                     'GRUMA' = 'GRANDW',
                     'GRWEN' = 'WEN',
                     'IRBSH' = c('CMP', 'BSC', 'LSHEEF'),
                     'IRMAI' = c('IR1_bb', 'IR3', 'CowCreek'),
                     'MFBEA' = 'BRC',
                     'MFBIG' = 'TAY',                 
                     'SCUMA' = 'SC1',
                     'SEUMA/SEMEA/SEMOO' = 'SW1',
                     'SFEFS' = 'ESS',
                     'SFSEC' = 'ZEN',
                     'SFSMA' = c('KRS', 'SFG_bb'),
                     'SNASO' = c('ACM', 'ALMOTC', 'ALPOWC', 'TENMC2'),
                     'SNTUC' = c('LTR', 'PENAWC'),
                     'SREFS' = 'SALEFT',
                     'SRLEM' = c('LLR','CRC'),
                     'SRLSR' = 'RAPH',
                     'SRNFS' = 'NFS',
                     'SRPAH' = 'PAHH',
                     'SRPAN' = 'PCA',
                     'SRLMA' = 'USI_bb',
                     'SRUMA' = 'STL',
                     'SRVAL' = 'VC2',
                     'SRYFS' = 'YFK') %>%
      stack() %>%
      tbl_df() %>%
      select(TRT = ind,
             site = values)
  }
  
  if(spp == 'Steelhead') {
    
    report_df = list('CRLOC-s' = 'LRL',
                     'CRLOL-s' = 'LC1',
                     'CRSEL-s' = 'SW1',
                     'CRSFC-s' = 'SC1',
                     'GRUMA-s' = c('UGR', 'LOOKGC'),
                     'GRJOS-s' = 'JOC',
                     'GRLMT-s' = 'WEN',
                     'GRWAL-s' = 'WR1',
                     'IRMAI-s' = c('IR1', 'COC'),
                     'SNASO-s' = c('ACM', 'PENAWC', 'ALMOTC', 'ALPOWC', 'TENMC2'),
                     'SREFS-s' = 'SALEFT', # need to think about USI_bb
                     'SRLEM-s' = c('LLR','CRC'),
                     'SRNFS-s' = 'NFS',
                     'SRPAH-s' = 'PAHH',
                     'SRPAN-s' = 'PCA',
                     'SRUMA-s' = c('YFK', 'VC2', 'STL'),
                     'SFSEC-s' = 'ZEN',
                     'SFMAI-s' = c('KRS', 'ESS', 'SFG_bb'),
                     'CRLMA-s' = c('JUL', 'LAP', 'CLC'),
                     'SNTUC-s' = 'LTR',
                     'MFBIG-s' = 'TAY',
                     'SRLSR-s' = 'RAPH') %>%
      stack() %>%
      tbl_df() %>%
      select(TRT = ind,
             site = values)
    
  }
  
  return(report_df)
  
}