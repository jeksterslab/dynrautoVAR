# Subgrouping Algorithm Based on Dynr Summary Output
  # 3 Major Tasks:
    # 1. Read in dynr summary object
    # 2. Significance Test
    # 3. Cluster based on either hard or fuzzy approaches
# library(fclust)
# inputs:
  # list of dynr.cook objects or output from dynr.var()
  # p.val and params are arguments required if supplying dynr.cook objects
    # p.val is for setting non-sig paths to 0.00
    # params is for selecting specific parameters by which you want to cluster on
  # var.opt
    # Argument when supplying dynr.var output.
      # Clustering can be done on the betas, residual covariance matrix, 
        # PDC's, PCC's, or both.
  # method
    # method for clustering can be hard (WalkTrap) or fuzzy (FKM)
    # If FKM is used, user must specify k [groups] and m [fuzzifier]
      # Defaults are 2 subgroups and 2 for fuzzifier
sub4dynr = function(inputs, type = c('dynr.cook', 'dynr.var'), 
                    p.val = 0.05, params = NULL,
                    var.opt = c('Betas', 'ResidCov', 
                                'PDC', 'PCC', 'Both'),
                    method = c('hard', 'fuzz'),
                    k = 2, m = 2){
  if(type == 'dynr.var'){
    comp.list = NULL
    N = length(inputs)
      if(var.opt == 'Both'){
        var.opt = c('PDC', 'PCC')
      }
    for(i in 1:N){
      comp.list = cbind(comp.list, unlist(inputs[[i]][var.opt]))
    }
    adj.mat = matrix(NA, N, N)
    for(i in 1:N){
      for(j in 1:N){
        adj.mat[i,j] = sum(abs(comp.list[,i]) > 0 &
                             abs(comp.list[,j]) > 0 &
          sign(comp.list[,i]) == sign(comp.list[,j]))
      }
    }
    adj.mat = adj.mat - min(adj.mat)
    diag(adj.mat) = 0
    if(method == 'hard'){
      grapp = graph_from_adjacency_matrix(adj.mat, weighted = TRUE, mode = 'undirected')
      sub.sol = cluster_walktrap(grapp)  
    }else{if(method == 'fuzz'){
      sub.sol = FKM(adj.mat,k=k,m=m)
    }}
    
  }
  if(type == 'dynr.cook'){
    N = length(inputs)
    comp.list = NULL
    for(i in 1:N){
      temp = summary(inputs[[i]])$Coefficients[,c(1,6)]
      noms = rownames(temp)
      comp.list = cbind(comp.list, temp[,1] * as.numeric(temp[,2] < p.val))
    }
    row.names(comp.list) = noms
    if(!is.null(params)){
      comp.list[which(row.names(comp.list) %in% params),]
    }else{
      message(paste0('Params is NULL so we are using all parameters provided to cluster.'))
      comp.list = comp.list      
    }

    adj.mat = matrix(NA, N, N)
    for(i in 1:N){
      for(j in 1:N){
        adj.mat[i,j] = sum(abs(comp.list[,i]) > 0 &
                             abs(comp.list[,j]) > 0 &
          sign(comp.list[,i]) == sign(comp.list[,j]))
      }
    }
    adj.mat = adj.mat - min(adj.mat)
    diag(adj.mat) = 0
    if(method == 'hard'){
      grapp = graph_from_adjacency_matrix(adj.mat, weighted = TRUE, mode = 'undirected')
      sub.sol = cluster_walktrap(grapp)  
    }else{if(method == 'fuzz'){
      sub.sol = FKM(adj.mat,k=k,m=m)
    }}
  }
  return(sub.sol)
}