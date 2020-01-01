#' @title Tukey's Trend Test
#' @description Two tests can be implemented based on permutation or based on asymptotic assumptions
#'
#' @param rank The rank output from DOOR.rank function, need to include two columns: "DOOR.rank" and "Dose". Make sure to have the actual dose levels for each subject, column name is "Dose", make sure it is numeric and controls have the lowest dose level
#' @param rsam Number of permutations to perform, default = 1000
#' @param alpha Desired type I error level, defalt = 0.05
#' @return The results of Tukey's Trend Test
#' @export

Tukey_trend_test=function(rank, #The rank output from DOOR.rank function, need to include two columns: "DOOR.rank" and "Dose". Make sure to have the actual dose levels for each subject, column name is "Dose", make sure it is numeric and controls have the lowest dose level
                          rsam = 1000, #Number of permutations to perform, default = 1000
                          alpha = 0.05 #Desired type I error level, defalt = 0.05
){
  library(mvtnorm)
  n_sub = as.numeric(table(rank$Dose)) #Number of subjects for each group, ordered by group dose level from small to large
  n_group = length(n_sub) #Number of groups


  dose_level = sort(unique(rank$Dose))
  m=rep(0,n_group)
  for(i in 1:n_group){
    m[i]=mean(rank$DOOR.rank[which(rank$Dose == dose_level[i])])
  }


  #Ordinal
  dose_level_ordinal = 0:(n_group-1)
  x = sum(dose_level_ordinal * n_sub) / sum(n_sub)
  c1 = n_sub * ( dose_level_ordinal - x )
  test_tukey_ordinal = sum(c1 * m)


  #Arithmetic
  dose_level_arith = dose_level
  dose_level_arith[1] = 0
  x = sum(dose_level_arith * n_sub) / sum(n_sub)
  c2 = n_sub * (dose_level_arith - x)
  test_tukey_arith = sum(c2 * m)


  #Logarithmic
  dose_level_arithlog=c(log10(dose_level[2])-(log10(dose_level[3])-log10(dose_level[2]))*(dose_level[2]-dose_level[1])/(dose_level[3]-dose_level[2]),log10(dose_level[-1]))
  x = sum( dose_level_arithlog * n_sub) / sum(n_sub)
  c3 = n_sub * (dose_level_arithlog - x)
  test_tukey_arithlog = sum(c3 * m)



  #Permutation Tukey's test
  perm=list()
  for (z in 1:rsam){
    rank_perm=rank
    rank_perm$Dose_sam=sample(rank$Dose, replace = F)

    m_perm=rep(0,n_group)
    for(i in 1:n_group){
      m_perm[i]=mean(rank_perm$DOOR.rank[which(rank_perm$Dose_sam == dose_level[i])])
    }


    test_tukey_ordinal_perm = sum(c1 * m_perm)
    test_tukey_arith_perm = sum(c2 * m_perm)
    test_tukey_arithlog_perm = sum(c3 * m_perm)
    perm[[z]]=list(test_tukey_ordinal_perm=test_tukey_ordinal_perm,test_tukey_arith_perm=test_tukey_arith_perm,test_tukey_arithlog_perm=test_tukey_arithlog_perm)
  }
  test_ordinal_perm=unlist(lapply(perm,"[[",1))
  test_arith_perm=unlist(lapply(perm,"[[",2))
  test_arithlog_perm=unlist(lapply(perm,"[[",3))


  p_perm_ordinal=1/(rsam+1)*(sum(test_ordinal_perm >= test_tukey_ordinal)+1)
  p_perm_arith=1/(rsam+1)*(length(which(test_arith_perm >= test_tukey_arith))+1)
  p_perm_arithlog=1/(rsam+1)*(length(which(test_arithlog_perm >= test_tukey_arithlog))+1)

  #Final p-value for Tukey's permutation test
  p_tukey_perm=min(p_perm_arith,p_perm_ordinal,p_perm_arithlog)


  #By asymptotic formula
  fit=lm(rank$DOOR.rank~rank$Dose)
  df = sum(n_sub)-n_group

  test_asymp_ordinal = sum(c1*m)^2 / (mean(fit$residuals^2) * sum(c1^2/n_sub))
  p_asymp_ordinal = 1-pf(test_asymp_ordinal,1,df)

  test_asymp_arithmetic = sum(c2*m)^2 / (mean(fit$residuals^2) * sum(c2^2/n_sub))
  p_asymp_arithmetic = 1-pf(test_asymp_arithmetic,1,df)

  test_asymp_arithlog = sum(c3*m)^2 / (mean(fit$residuals^2) * sum(c3^2/n_sub))
  p_asymp_arithlog = 1-pf(test_asymp_arithlog,1,df)

  #Final p-value for Tukey's test based on asymptotic formula
  tukey_p_asymp=min(p_asymp_ordinal,p_asymp_arithmetic,p_asymp_arithlog)


  #To adjust for tukey's trend test with the new alpha threshold at the desired type I error level, by default type I error level = 0.05
  cor12=cor(dose_level_ordinal,dose_level_arith)
  cor13=cor(dose_level_ordinal,dose_level_arithlog)
  cor23=cor(dose_level_arith,dose_level_arithlog)
  cormat1=matrix(c(1,cor12,cor13,cor12,1,cor23,cor13,cor23,1),3,3)
  #If desired type 1 error set to 0.05
  t=as.numeric(qmvt(1-alpha,tail="both.tails",maxiter=1000,corr=cormat1,interval=c(0,1),df=df)[1])
  alpha_adjusted=pt(t,df=df,lower.tail = FALSE)*2



  results = list(tukey_perm_statistics = perm, p_tukey_perm = p_tukey_perm, tukey_asymp_statistics = c(test_asymp_ordinal = test_asymp_ordinal,test_asymp_arithmetic = test_asymp_arithmetic,test_asymp_arithlog = test_asymp_arithlog ),tukey_p_asymp = tukey_p_asymp, alpha_adjusted = alpha_adjusted)
  return(results)


}
