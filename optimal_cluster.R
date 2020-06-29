conc_pvalue <- function(pval_cluster, m_cluster){
  
  sum(pval_cluster)/m_cluster
  
}

optimal_cluster <- function(map, coord, hc, min_h, max_h, by_h, clusters, pvalues){
  
  h <- seq(min_h, max_h, by_h)
  
  K <- rep(0,length=length(h))
  
  G <- list()
  
  for(i in 1:length(h)){
    
    ct_h <- cutree(hc, h=h[i])
    
    K[i] <- length(unique(ct_h))
    
    conc <- rep(0, K[i])
    
    attivation <- list()
    
    clusters_h = clusters
    
    clusters_h[clusters] = ct_h
    
    for (j in 1:nrow(coord)){
      
      map[coord[j,1],coord[j,2],coord[j,3]] <- clusters_h[j]
      
    }
    
    for(k in 1:K[i]){
      
      id_selected = which( map[mask3D==1]==k)
      
      discov_cluster <- discoveries(pvalues,alpha=0.05, ix=id_selected)
      
      tdp_cluster <- round(tdp(pvalues,alpha=.05, ix=id_selected)*100,1)
      
      pval_cluster <- pvalues@adjusted[id_selected] <= 0.05
      
      m_cluster <- length(pval_cluster)
      
      conc[k] <- conc_pvalue(pval_cluster, m_cluster)
      
      attivation[[k]] <- c(discov_cluster, tdp_cluster)
      
    }
    G[[i]] <-  c(gini=list(Gini(conc)), attivation=list(attivation), h=h[i])
  }
  
  return(G)
}