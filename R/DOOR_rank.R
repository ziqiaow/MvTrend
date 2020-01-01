#' @title DOOR algorithm for ranking (Based on the algorithm of Evans et al (2015))
#' @description The higher DOOR.rank, the better clinical outcome the patient has (based on benefits and harms)
#'
#' @param df Dataset that contains all the variables that should be included for ranking and group label, class data.frame
#' @param var_order Desired order of variables in dataset df
#' @param ties.method By default, ties.method is average
#' @return The rank based on DOOR algorithm based on multivariate clinical variables
#' @export


DOOR.rank <- function( df, #Dataset that contains all the variables that should be included for ranking and group label, class data.frame
                       var_order, #Desired order of variables in dataset df
                       ties.method = c("average","random") #By default, ties.method is average
) {
  library(dplyr)
  if (is.null(ties.method) == 1) {ties.method <- "average"}
  average <- function(x) mean(x)
  random <- function(x) sample(x, length(x))
  out <- df %>% mutate(r = do.call(order, c(df[,var_order])), rr=order(r)) %>% group_by_(.dots = var_order,add=T) %>% mutate(DOOR.rank = get(ties.method)(rr)) %>% as.data.frame() %>% select(-r, -rr)
  rownames(out)=rownames(df)
  return(out)
}
