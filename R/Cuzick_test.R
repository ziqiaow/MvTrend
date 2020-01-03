#' @title Cuzick's Test for Trend
#' @description Two tests can be implemented based on permutation or based on asymptotic assumptions
#'
#' @param rank The rank output from DOOR.rank function, need to include two columns: "DOOR.rank" and "label", "label" should be the ordinal group for each individual, need to be numeric from 1,2,...n
#' @param alternative Which type of test, default is "greater"
#' @param correct.ties Whether to correct ties when ties exist, default is TRUE
#' @param verbose Print whether ties exist or not
#' @param perm Whether to conduct permutation test, default is FALSE
#' @param rsam Number of permutations to perform, default = 1000
#' @return The results of Cuzick's Test for Trend
#' @export


cuzick.test=function(rank, #The rank output from DOOR.rank function, need to include two columns: "DOOR.rank" and "label", "label" should be the ordinal group for each individual, need to be numeric 1,2,...n
                     alternative = c("greater","less","two.sided"), #Which type of test, default is "greater"
                     correct.ties=TRUE, #Whether to correct ties when ties exist, default is TRUE
                     verbose = FALSE, #Print whether ties exist or not
                     perm = FALSE, #Whether to conduct permutation test, default is FALSE
                     rsam = 1000 #Number of permutations
){
  x=rank$DOOR.rank;g=rank$label
  alternative <- match.arg(alternative)
  t = sum(g*x) #Estimate
  n = length(x)  #number of observations
  ngroup = length(unique(g)) #number of groups
  expz = sum(as.numeric(table(g)) * 1:ngroup)/n
  expt = expz * n * (n + 1) / 2 #Expectation of the test statistics
  varz = sum(as.numeric(table(g)) * (1:ngroup)^2)/n - expz^2

  ## check for ties
  TIES <- FALSE
  TIES <- sum(table(x) - 1) != 0

  if (correct.ties == TRUE){

    if (verbose ==TRUE & TIES ) {print("Ties are present. Correct for ties.")}
    if (verbose ==TRUE & TIES == FALSE) {print("Ties are not present.")}
    n_ties=as.numeric(table(x))
    phi = (sum(n_ties^3-n_ties)) / (n * (n^2 - 1))
    vart = n^2*(n+1)/12*varz * (1 - phi)

  } else {

    if (verbose ==TRUE & TIES == FALSE) {print("Ties are not present.")}
    if (TIES == TRUE) {warning("Ties are present, but do not correct for ties.")}
    vart = n^2*(n+1)/12*varz

  }

  zscore=(t-expt)/sqrt(vart)

  if(alternative == "greater") {
    p = pnorm(zscore,lower.tail = FALSE)
  } else if (alternative =="less") {
    p = pnorm(zscore)
  } else {
    p = 2 * min(pnorm(abs(zscore),lower.tail = FALSE),0.5)
  }


  if (perm == TRUE){

    t.perm = 0
    for (z in 1:rsam){
      rank_perm=rank
      rank_perm$label_perm=sample(rank$label, replace = F)
      t.perm[z] = sum(x*rank_perm$label_perm) #Estimate
    }

    perm.p.value=1/(rsam+1)*(sum(t.perm >= t)+1)

  } else {t.perm=NULL;perm.p.value=NULL}


  results=list(estimates = t, variance = vart, z = zscore, p.value = p, test.type = alternative, estimates.perm = t.perm, perm.p.value = perm.p.value)
  return(results)
}




