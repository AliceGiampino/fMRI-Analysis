threshold <- function(stat, min_t, max_t, by_t, size, coord, k){
  
  thresh <- seq(min_t, max_t, by_t)
  
  tdp_mean <- NULL
  
  for(t in 1:length(thresh)){
    
    clusters_t = cluster.threshold(stat, nmat=nmat, 
                                   level.thr = thresh[t], 
                                   size.thr = size)
    clusters_t = clusters_t[mask3D==1]==1
    
    dd_t = dist(coord[clusters_t,])
    hc = hclust(dd_t, "single")
    
    ct_t = cutree(hc, k=k)
    
    clusters2 = clusters_t
    clusters2[clusters_t] = ct_t
    
    map <- mask3D
    for (i in 1:nrow(coord)){
      map[coord[i,1],coord[i,2],coord[i,3]] <- clusters2[i]
    }
    
    tdp <- rep(0, length(4))
    
    for (j in 1:k){
      id_selected=which( map[mask3D==1]==j)
      cat(length(id_selected),"\n")
      tdp[j] <- round(tdp(res,alpha=.05, ix=id_selected)*100,1)
    }
    
    tdp_mean[[t]] <- list(mean_tdp = mean(tdp))
    
  }
  
  return(tdp_mean)
  
}