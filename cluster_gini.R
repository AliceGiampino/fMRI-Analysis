
Cluster_Gini <- function(hc, min_h, max_h, by_h, clusters){
  
  cluster_h <- clusters
  
  h <- seq(min_h, max_h, by_h)
  
  K <- rep(0,length=length(h))
  
  G <- list()
  
  for(i in 1:length(h)){
    ct_h <- cutree(hc, h=h[i])
    
    cluster_h[clusters] <- ct_h
    
    K[i] <- length(unique(ct_h))
    
    print(c(K[i], h[i]))
    
    coef_Gini_in <- NULL
    coef_Gini_out <- NULL
    
    mean_ck <- rep(0, K[i])
    
    for(k in 1:K[i]){
      ck <- cluster_h==k
      coef_Gini_in[k] <- Gini(stat[which(stat!=-1)[ck]])
      
      if(coef_Gini_in[k]=="NaN"){
        coef_Gini_in[k] <- 0
      }
      
      mean_ck[k] <- mean(stat[which(stat!=-1)[ck]])
    }
    
    coef_Gini_out <- Gini(mean_ck)
    
    G[[i]] <- c(gini_in = list(coef_Gini_in), gini_out = list(coef_Gini_out))
  }
  return(G)
}
