## ----setup,include=F----------------------------------------------------------
library(MvTrend)

## -----------------------------------------------------------------------------
# Load the example data
data("example_trend")
head(example_trend)

## -----------------------------------------------------------------------------
# Calculate the rank
rank.output <- DOOR.rank( example_trend, var_order = c("Severe","Benefit","Moderate","Mild") )
head(rank.output)

## -----------------------------------------------------------------------------
#Jonckheere Terpstra Test
JT.test.output <- JT_test(rank.output, rsam = 1000)
objects(JT.test.output)
JT.test.output$p_jon_perm #The p value for permutation test
JT.test.output$p_jon_asymp #The p value for the nonparametric test based on asymptotic assumptions

## -----------------------------------------------------------------------------
tukey.test.output <- Tukey_trend_test(rank.output, rsam = 1000, alpha = 0.05) #alpha is the desired type I error level, defalt is 0.05
objects(tukey.test.output)
tukey.test.output$p_tukey_perm #The p value for permutation test
tukey.test.output$p_tukey_asymp #The p value for the parametric test based on asymptotic assumptions
tukey.test.output$alpha_adjusted #The adjusted alpha threshold to reject the test

## -----------------------------------------------------------------------------
cuzick.test.output <- cuzick.test(rank.output, alternative = "greater", perm = TRUE)
objects(cuzick.test.output)
cuzick.test.output$perm.p.value #The p value for permutation test
cuzick.test.output$p.value #The p value for the nonparametric test based on asymptotic assumptions


