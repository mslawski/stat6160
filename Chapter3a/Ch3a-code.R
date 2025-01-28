##########################################################################
#
#
# R-code for Ch.~3a: One-way Analysis of Variance
#
#
##########################################################################
ls()
##########################################################################
# Power Setting vs. Etch Rate Experiment
##########################################################################
power_setting <- rep(c("160W", "180W", "200W", "220W"), each = 5)
etchrate  <- c(575, 542, 530, 539, 570, 565, 593, 590, 579, 610, 600, 651, 610, 637, 629, 725, 700, 715, 685, 710)

data <- data.frame(etchrate = etchrate, power_setting = as.factor(power_setting))
contrasts(data$power_setting) <- contr.sum(4)

# visual representation
boxplot(etchrate ~ power_setting, data = data)

# generate ANOVA table
anova_result <- aov(data = data, formula = etchrate ~ power_setting) 
print(anova_result)
summary(anova_result)

# diagnostics

# residuals vs. fitted
plot(fitted(anova_result), residuals(anova_result), xlab = "fitted", ylab = "residuals", pch = 16)

# standardized residuals 
plot(rstandard(anova_result), ylim = c(-4, 4), pch = 16)
abline(h = c(-3,3), col = "red")

# Normal Q-Q plot
qqnorm(rstandard(anova_result), pch = 16)
abline(0,1)

# Bartlett test for equality of variance
bartlett.test(etchrate ~ power_setting, data = data)

# (for reference: hand-calculation according to the expressions on the slides)
Sisq <- tapply(data$etchrate, INDEX = data$power_setting, FUN = var)
ni <-  tapply(data$etchrate, INDEX = data$power_setting, FUN = length)
k <- nlevels(data$power_setting)
N <- sum(ni)
Ssq <- sum((ni-1)* Sisq)/(N-k)
Q <- (N-k) * log10(Ssq) - sum((ni-1)*log10(Sisq))
C <- 1 + 1/(3 * (k-1)) * sum((1/(ni-1) - 1/(N-k)))
Chi0sq <- 2.3026 * Q / C
1 - pchisq(Chi0sq, df = k - 1) # p-value
# (results are not exactly identical but close) 

# The Levene test is not pre-implemented in R (but can be accessed via the package car)

# we here hand-calculate the Levene test statistic:

# calculate the absolute differences from the median in each group
Dij <- tapply(data$etchrate, INDEX = data$power_setting, FUN = function(x) abs(x - median(x)))
# run ANOVA on D_ij
Dij_data <- data.frame(D = unlist(Dij), power_setting = power_setting)
summary(aov(D ~ power_setting, data = Dij_data))
# again, p-value is large.


### (post-hoc) power calculation (for F-test) 

# we have 5 samples per group
power.anova.test(groups = k, n = 5, between.var = mean((fitted(anova_result) - mean(fitted(anova_result)))^2),
                 within.var = sigma(anova_result)^2)

# even with only 1/3 of between variance as observed, the power would have been nearly one: 
power.anova.test(groups = k, n = 5, between.var = 1000,
                 within.var = sigma(anova_result)^2)

# with about 1/10 of the between variance as observed, the power is still about .79
power.anova.test(groups = k, n = 5, between.var = 300,
                 within.var = sigma(anova_result)^2)

# now suppose that the within variance is four times higher (standard deviation 2 times higher)
power.anova.test(groups = k, n = 5, between.var = 300,
                 within.var = 4 * sigma(anova_result)^2)
# now the power is kind of low (.25)
