##########################################################################
#
#
# R-code for Ch.~8: Random Effects Model for Single-Factor Experiments
#
#
##########################################################################
ls()
rm(list = ls())
##########################################################################
# Loom Data 
##########################################################################
a <- 4 # number of factor levels
n <- 4 # number of replicates for each factor
loom <- factor(rep(c("1", "2", "3", "4"), each = n))
y <- c(98, 97, 99, 96,
       91, 90, 93, 92,
       96, 95, 97, 95,
       95, 96, 99, 98)
N <- n* a

# standard ANOVA:

aov_fixed <- aov(y ~ loom, contrasts = list(loom = "contr.sum"))
summary(aov_fixed)

# NOTE: we can re-use a good amount of quantities from 
# the anova output, e.g., the F-statistic and MSE (estimator 
# of sigma^2):

sigmasq_hat <- sigma(aov_fixed)^2

# the latter allows us to estimate sigmatau_sq:  
aov_fixed_sum <-  unlist(summary(aov_fixed)) 
sigmatausq_hat <- (aov_fixed_sum["Mean Sq1"] -  sigmasq_hat)/n

# ICC (intra-clas correlation coefficient):
ICC <- sigmatausq_hat / (sigmatausq_hat + sigmasq_hat)  

### confidence interval for intra-class correlation coefficient
F0 <- aov_fixed_sum["Mean Sq1"] / aov_fixed_sum["Mean Sq2"]
alpha <- 0.05
critF_lower <- qf(alpha/2, df1 = a-1, df2 = N - a)
critF_upper <- qf(1 - alpha/2, df1 = a-1, df2 = N - a)
Upper <- (F0 - critF_lower)/(F0 + (n - 1) * critF_lower)
Lower <- (F0 - critF_upper)/(F0 + (n - 1) * critF_upper)

print(paste("Lower CI limit: ", signif(Lower, 3), " Upper CI limit: ", signif(Upper, 3), sep = ""))

### confidence interval for sigma_sq
Upper_sigma_sq <- aov_fixed_sum["Sum Sq2"] / qchisq(alpha/2, N-a) 
Lower_sigma_sq <- aov_fixed_sum["Sum Sq2"] / qchisq(1-alpha/2, N-a)

print(paste("Lower CI limit: ", signif(Lower_sigma_sq, 3), " Upper CI limit: ", signif(Upper_sigma_sq, 3), sep = ""))

### grand mean and confidence interval 
ybar <- mean(y)
varhat_ybar <- aov_fixed_sum["Mean Sq1"]/(n*a)  
crit_t <- qt(1-alpha/2, df = a - 1)
Lower_mean <- ybar - crit_t * sqrt(varhat_ybar)
Upper_mean <- ybar + crit_t * sqrt(varhat_ybar)


### Analysis using the lme4 package in R ---
### this package implements a more modern and general approach
### to parameter estimation and inference.

library(lme4)
mem <- lmer(y ~ 1 + (1|loom)) 
sigmasq_reml <- summary(mem)$sigma^2 # exactly identical to sigmasq_hat
sigmatausq_reml <- as.numeric(summary(mem)$varcor) # almost exactly identical to momemt-based estimator

### we can obtain asymptotic confidence intervals using the profile function:
prof_mem <- profile(mem)
confint(prof_mem) # all parameters, including the mean


confint(prof_mem)[1:2,]^2

# we note that the exact Chi^2-interval is wider. Since the number
# of sample is small, the underlying asymptotics in the REML 
# estimation method are not justified. 
