query_inat <- function(x) {
  sp_bi=vector("list", nrow(x))
  
  for(i in 1:nrow(x)) {
    site.coord = x[i, ]
    loc.name = x[i, ]
    
    sp_bi[[i]] = NA
    boundary = c(site.coord[1] - 0.01,
                 site.coord[2] - 0.01,
                 site.coord[1] + 0.01,
                 site.coord[2] + 0.01) ## Approximately a 1000x1000m area. Distorted at higher latitudes
    #print(boundary)
    
    out = tryCatch({
      #require(rinat)
      #require(rgbif)
      
      inatq = get_inat_obs(taxon_name = "Viridiplantae",
                           bounds = boundary,
                           maxresults  = 5000)
      inatq = inatq[which(inatq$captive_cultivated == 'false'),]
      inatq = inatq[which(inatq$quality_grade == 'research'),]
      inatq = inatq[which(inatq$positional_accuracy < 200),]
      
      key <- name_backbone(name = '', kingdom = 'plantae')$kingdomKey
      s.gbif = occ_data(kingdomKey = key,
                        geometry = boundary[c(2, 1, 4, 3)],
                        limit = 2000)
      t_names = unique(c(unique(inatq$scientific_name),
                         unique(s.gbif$data$name)))
      print(i); print(t_names)
      sp_bi[[i]] = t_names[grep(' ', t_names)]
      print(sp_bi[[i]])
      
    },
    error = function(cond) {
      # Choose a return value in case of error
      sp_bi[[i]] = NA
    }, warning = function(cond) {
      sp_bi[[i]] = NA
    })
  }
  print(nrow(x))
  print(i)
  print(length(sp_bi))
  return(sp_bi)
}
