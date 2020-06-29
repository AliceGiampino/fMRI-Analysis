threshold_k <- function(stat, min_t, max_t, by_t, size, coord, k){
  
  thresh <- seq(min_t, max_t, by_t)
  
  tdp_mean <- NULL
  
  for(t in 1:length(thresh)){
    
    clusters_t = cluster.threshold(stat, nmat=nmat, 
                                   level.thr = thresh[t], 
                                   size.thr = size)
    clusters_t = clusters_t[mask3D==1]==1
    
    dd_t = dist(coord[clusters_t,])
    hc = hclust(dd_t, "single")
    
    k = k
    
    tdp_mean_k <- NULL
    
    for(num in 1:length(k)){
      
      ct_t = cutree(hc, k=k[num])
      
      clusters2 = clusters_t
      clusters2[clusters_t] = ct_t
      
      map <- mask3D
      for (i in 1:nrow(coord)){
        map[coord[i,1],coord[i,2],coord[i,3]] <- clusters2[i]
      }
      
      tdp <- rep(0, length(4))
      
      for (j in 1:length(k[num])){
        id_selected=which( map[mask3D==1]==j)
        cat(length(id_selected),"\n")
        tdp[j] <- round(tdp(res,alpha=.05, ix=id_selected)*100,1)
      }
      
      tdp_mean_k[[num]] <- list(mean_tdp = mean(tdp), k = k[num])
      
    }
    
    tdp_mean[[t]] <- list(mean_tdp = tdp_mean_k)
    
  }
  
  return(tdp_mean)
  
}