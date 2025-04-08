##########################################################################
#
#
# R-code for Ch.~8: Random Effects and Mixed Effects Models
#                   for factor experiments with two factors
#
#
##########################################################################
ls()
rm(list = ls())
##########################################################################
# Gauge Data 
##########################################################################
a <- 20  # number of parts
b <- 3   # number of operators
n <- 2   # number of replicates for each factor 
N <- a * b * n # total number of measurements 

### create data set 
operator_no <- rep(c(1:b), each = a) 
part_no <- rep(1:a, times = b)

# replicate 1
y1 <- c(21, 24, 20, 27, 19, 23, 22, 19, 24, 25, 21, 18, 23, 24, 29, 26, 20, 19, 25, 19,
        20, 24, 19, 28, 19, 24, 22, 18, 25, 26, 20, 17, 25, 23, 30, 25, 19, 19, 25, 18,
        19, 23, 20, 27, 18, 23, 22, 19, 24, 24, 21, 18, 25, 24, 31, 25, 20, 21, 25, 19)

dat1 <- cbind(part_no, operator_no, y1)


# replicate 2    
y2 <- c(20, 23, 21, 27, 18, 21, 21, 17, 23, 23, 20, 19, 25, 24, 30, 26, 20, 21, 26, 19,
        20, 24, 21, 26, 18, 21, 24, 20, 23, 25, 20, 19, 25, 25, 28, 26, 20, 19, 24, 17,
        21, 24, 22, 28, 21, 22, 20, 18, 24, 25, 20, 19, 25, 25, 30, 27, 20, 23, 25, 17)


dat2 <- cbind(part_no, operator_no, y2)

dat <- rbind(dat1, dat2)
colnames(dat) <- c("part_no", "operator_no", "y")
data <- data.frame(part = as.factor(dat[,"part_no"]),
                   operator = as.factor(dat[,"operator_no"]),
                   y = dat[,"y"])

# standard ANOVA:

aov0 <- aov(y ~ part*operator, contrasts = list(part = "contr.sum", operator = "contr.sum"),
            data = data)
summary(aov0)

aov0_sum <- unlist(summary(aov0))

# estimation of variance components (assuming a two-way random effects model): 

sigmasq_hat <- sigma(aov0)^2

sigmataubetasq_hat <-  (aov0_sum["Mean Sq3"] - sigmasq_hat) / n # negative variance estimate
sigmabetasq_hat <- (aov0_sum["Mean Sq2"] - aov0_sum["Mean Sq3"])/(a * n)
sigmatausq_hat <- (aov0_sum["Mean Sq1"] - aov0_sum["Mean Sq3"])/(b * n)

# Note that the F-tests from the standard ANOVA table 
# are !*not*! the F-tests for testing variance components  

# F-tests and associated p-values: 
Ftaubeta <- aov0_sum["Mean Sq3"] / sigmasq_hat; p_taubeta <- 1-pf(Ftaubeta, df1 = aov0_sum["Df3"], df2 = aov0_sum["Df4"])
print(p_taubeta) # This F-test is the same as in the standard two-way ANOVA.

# Different from standard two-way ANOVA
F_beta <- aov0_sum["Mean Sq2"] / aov0_sum["Mean Sq3"]; p_beta <-  1 - pf(F_beta, df1 = aov0_sum["Df2"], df2= aov0_sum["Df3"])
print(p_beta)

# Different from standard two-way ANOVA
F_tau <- aov0_sum["Mean Sq1"] / aov0_sum["Mean Sq3"]; p_tau <- 1 - pf(F_tau, df1 = aov0_sum["Df1"], df2= aov0_sum["Df3"])   
print(p_tau)

################################################################################################
# simplifying to a main effect model: 
aov1 <- aov(y ~ part + operator, contrasts = list(part = "contr.sum", operator = "contr.sum"),
            data = data)

sigmasq_hat <- sigma(aov1)^2

summary(aov1)
aov1_sum <- unlist(summary(aov1))

# estimation of variance components (again, assuming random effects)

sigmabetasq_hat <- (aov1_sum["Mean Sq2"] - sigmasq_hat)/(a * n)
sigmatausq_hat <- (aov1_sum["Mean Sq1"] - sigmasq_hat)/(b * n)

# F-tests are the same as in the fixed factor w/ main effects model

################################################################################################

### Inference under the RESTRICTED mixed effects model, treating
### ``operator" as fixed factor and ``part" as random factor
### Variance components and F-tests

# (NOTE: "tau" and "beta" are flipped compared to the slides; here the tau's are the random effects) 

# Restoring the MS_E from the interaction model: 
sigmasq_hat <- sigma(aov0)^2

### Inference under the UNRESTRICTED mixed effects model:

# unchanged 
sigmataubetasq_hat <-  (aov0_sum["Mean Sq3"] - sigmasq_hat) / n
# changed compared to the two-way random effects case
sigmatausq_hat <- (aov0_sum["Mean Sq1"] - sigmasq_hat)/(b * n)
# F-test changes accordingly: 
Ftau <- aov0_sum["Mean Sq1"] / sigmasq_hat; p_tau <- 1 - pf(Ftau, df1 = aov0_sum["Df1"], df2 = aov0_sum["Df4"])
print(p_tau)
# F-test for fixed effects is the same as "Fbeta" above in the fully random effects case  

### Inference under the RESTRICTED mixed effects model:

### sigmatausq_hat reverts back to its estimate under the fully random effect case,
### and so does the corresponding F-statistic. 
