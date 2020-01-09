# this function defines which estimates from DABOM correspond to which populations

definePopulations = function(spp = c('Chinook', 'Steelhead')) {
  
  spp = match.arg(spp)
  
  if(spp == 'Chinook') {
    # NEEDS FIXED
    report_df = list('CRLAP' = 'Lapwai',
                     'CRLOC' = 'Lochsa',
                     'CRLOL' = 'Lolo',
                     'CRPOT' = 'Potlatch',
                     'GRCAT' = 'CATHEW',
                     'GRLOO' = 'LookingGlass',
                     'GRLOS' = c('Wallowa_bb', 'past_WR2'),
                     'GRMIN' = 'past_MR1',
                     'GRUMA' = 'GRANDW',
                     'GRWEN' = 'Wenaha',
                     'IRBSH' = c('past_CMP', 'past_BSC', 'LSHEEF'),
                     'IRMAI' = c('ImnahaRiver_bb', 'past_IR3', 'CowCreek'),
                     'MFBEA' = 'BearValley',
                     'MFBIG' = 'BigCreek',                 
                     'SCUMA' = 'SFClearwater',
                     'SEMEA' = 'Selway',
                     'SFEFS' = 'past_ESS',
                     'SFSEC' = 'past_ZEN',
                     'SFMAI' = c('past_KRS', 'SFSalmon_bb'),
                     'SNASO' = c('Asotin', 'Almota', 'Alpowa', 'TenMileCreek'),
                     'SNTUC' = c('Tucannon', 'Penawawa'),
                     'SREFS' = 'SALEFT',
                     'SRLEM' = c('Lemhi','CarmenCreek'),
                     'SRLSR' = 'RapidRiver',
                     'SRNFS' = 'NFSalmon',
                     'SRPAH' = 'PAHH',
                     'SRPAN' = 'Panther',
                     'SRLMA' = 'USI_bb',
                     'SRUMA' = 'past_STL',
                     'SRVAL' = 'past_VC2',
                     'SRYFS' = 'past_YFK') %>%
      stack() %>%
      tbl_df() %>%
      select(TRT = ind,
             area = values)
  }
  
  if(spp == 'Steelhead') {
    
    report_df = list('CRLOC-s' = 'Lochsa',
                     'CRLOL-s' = 'Lolo',
                     'CRSEL-s' = 'Selway',
                     'CRSFC-s' = 'SFClearwater',
                     'GRUMA-s' = c('GrandeRonde', 'LookingGlass'),
                     'GRJOS-s' = 'JosephCreek',
                     'GRLMT-s' = 'Wenaha',
                     'GRWAL-s' = 'Wallowa',
                     'IRMAI-s' = c('ImnahaRiver', 'CowCreek'),
                     'SNASO-s' = c('Asotin', 'Penawawa', 'Almota', 'Alpowa', 'TenMileCreek'),
                     'SREFS-s' = 'SALEFT', # need to think about USI_bb
                     'SRLEM-s' = c('Lemhi','CarmenCreek'),
                     'SRNFS-s' = 'NFSalmon',
                     'SRPAH-s' = 'PAHH',
                     'SRUMA-s' = c('past_YFK', 'past_VC2', 'past_STL'),
                     'SFSEC-s' = 'past_ZEN',
                     'SFMAI-s' = c('past_KRS', 'past_ESS', 'SFSalmon_bb'),
                     'CRLMA-s' = c('Potlatch', 'Lapwai', 'ClearCreek'),
                     'SNTUC-s' = 'Tucannon',
                     'MFBIG-s' = 'BigCreek',
                     'SRLSR-s' = 'RapidRiver') %>%
      stack() %>%
      tbl_df() %>%
      select(TRT = ind,
             area = values)
    
  }
  
  return(report_df)
  
}

