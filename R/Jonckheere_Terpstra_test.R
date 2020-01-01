#' @title Jonckheere-Terpstra Test
#' @description Two tests can be implemented based on permutation or based on asymptotic assumptions
#'
#' @param rank The rank output from DOOR.rank function, need to include two columns: "DOOR.rank" and "label", "label" should be the ordinal group for each individual, need to be numeric from 1,2,...n
#' @param rsam Number of permutations to perform, default = 1000
#' @return The results of Jonckheere-Terpstra Test

JT_test=function(rank, #The rank output from DOOR.rank function, need to include two columns: "DOOR.rank" and "label", "label" should be the ordinal group for each individual, need to be numeric from 1,2,...n
                 rsam = 1000 #Number of permutations to perform, default = 1000
){


  n_sub = as.numeric(table(rank$label)) #Number of subjects for each group, ordered by group label level from small to large
  n_group = length(n_sub) #Number of groups
  dose_level = sort(unique(rank$label))


  #Construct the u statistics
  ustat=function(x,y){
    u=matrix(nrow=dim(x)[1],ncol=dim(y)[1])
    for (i in 1:dim(x)[1]){
      for (j in 1:dim(y)[1]){

        if(x$DOOR.rank[i]<y$DOOR.rank[j]) {
          u[i,j]=1
        } else if(x$DOOR.rank[i]==y$DOOR.rank[j]) {
          u[i,j]=0.5
        } else {u[i,j]=0}

      }
    }
    return(u)
  }


  #By aymptotic formula

  u_stat = 0

  for(j in 1:(n_group-1)){
    for(z in 1:(n_group-j)){
      u_stat = c(u_stat,sum(ustat(rank[which(rank$label == dose_level[j]),],rank[which(rank$label == dose_level[j+z]),])))
    }
  }

  u_stat = u_stat[-1]
  u_all = sum(u_stat)
  N = sum(n_sub)

  expectu = (N^2-sum(n_sub^2))/4

  tie = as.numeric(table(rank$DOOR.rank))

  varu = ((N*(N-1)*(2*N+5)-sum(n_sub*(n_sub-1)*(2*n_sub+5)) - sum(tie*(tie-1)*(2*tie+5)))/72)+
    (sum(n_sub*(n_sub-1)*(n_sub-2))+sum(tie*(tie-1)*(tie-2)))/(36*N*(N-1)*(N-2))+
    (sum(n_sub*(n_sub-1))+sum(tie*(tie-1)))/(8*N*(N-1))

  z_u = (u_all-expectu) / sqrt(varu)
  p_jon_asymp = pnorm(z_u,lower.tail=F)


  #Permutation Jonckheere test
  library(clinfun)
  joncktest_perm = jonckheere.test(x=rank$DOOR.rank,g=rank$label,alternative = "increasing",nperm = rsam)
  test_stat_jon_perm = joncktest_perm$statistic
  p_jon_perm = joncktest_perm$p.value

  results = list(JT_score_asymp = u_all, test_stat_jon_asymp = z_u, variance_jon_asymp = varu, p_jon_asymp = p_jon_asymp, JT_score_perm = test_stat_jon_perm, p_jon_perm = p_jon_perm, jon_perm_all = joncktest_perm)
  return(results)
}


