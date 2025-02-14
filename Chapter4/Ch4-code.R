##########################################################################
#
#
# R-code for Ch.~4: Block Design
#
#
##########################################################################
ls()
rm(list = ls())
##########################################################################
# Highway Paint Data 
##########################################################################

data <- data.frame(site = as.factor(rep(1:6, each = 4)), supplier = as.factor(rep(c("GS", "FD", "L", "ZK"), times = 6)),
                   response = c(69, 59, 55, 70,
                                83, 65, 65, 75,
                                74, 64, 59, 74,
                                61, 52, 59, 62,
                                78, 71, 67, 74,
                                69, 64, 58, 74))

# The above data frame was created in long format, which is more suitable for subsequent statistical analyses.
# Here, we re-create the corresponding table (as shown in the slides). This is an important double-check to
# verify that data were entered correctly. 

tapply(data, INDEX = list(data$site, data$supplier), FUN = function(z) z$response)[, c("GS", "FD", "L", "ZK")]


### enforce sum-to-zero constraints on the coefficients prior to running anova 

contrasts(data$site) <- contr.sum(6)
contrasts(data$supplier) <- contr.sum(4)


# generate ANOVA table
anova_result <- aov(data = data, formula = response ~ supplier + site) 
print(anova_result)
summary(anova_result)

# align the above anova table with the expressions shown on the slides
N <- nrow(data)  
SST <- var(data$response) * (N-1) # total sum of squares

ybar <- mean(data$response) 
sum((data$response - ybar)^2)

# 
SS_all <- summary(anova_result)[[1]][,"Sum Sq"]
all.equal(sum(SS_all), SST)

#
SSE <- sum(residuals(anova_result)^2)
all.equal(SSE, SS_all[3])

# SS_Bl, SS_Tr
ybar <- mean(data$response)

ybar_sites <- tapply(data, data$site, function(z) mean(z$response))
ybar_suppliers <- tapply(data, data$supplier, function(z) mean(z$response))

n <- nlevels(data$site)
k <- nlevels(data$supplier)

SS_Bl <- k * sum((ybar_sites - ybar)^2)
SS_Tr <- n * sum((ybar_suppliers - ybar)^2)

all.equal(SS_Bl, SS_all[2])
all.equal(SS_Tr, SS_all[1])

### pair-wise comparisons of suppliers 

gr <- expand.grid(1:4, 1:4)
gr <- gr[gr[,1] > gr[,2],] 

Contr <- t(contr.sum(4))

coef_suppliers <- coef(anova_result)[2:4]
# pairwise mean differences 
Deltas <- apply(gr, 1, function(z) sum(Contr[,z[1]] * coef_suppliers) - sum(Contr[,z[2]] * coef_suppliers))

V <- vcov(anova_result)[2:4, 2:4]

ses <- sqrt(apply(gr, 1, function(z) t(Contr[,z[1]] - Contr[,z[2]]) %*% V %*% (Contr[,z[1]] - Contr[,z[2]])))   
# standard errors are all identical, given balanced group sizes 

Ts <- Deltas / ses
names_suppliers <- levels(data$supplier)
names(Ts) <- apply(gr, 1, function(z) paste(names_suppliers[z[1]], names_suppliers[z[2]], sep = "vs"))

# Bonferroni cut-off
qt(1 - 0.025/choose(4,2), df = nrow(data) - n - k + 1)

# Tukey cut-off 
qtukey(1 - 0.05, k, nrow(data) - n - k + 1)/sqrt(2)

# Model adequacy checking
rstand <- rstandard(anova_result)
x11()
qqnorm(rstand, main = "Q-Q plot of residuals", pch = 16)
abline(0,1)

# they are two observations with larger residuals 
plot(rstand, ylim = c(-3, 3))

# what are these larger residuals:
# locator()
# Number 5, 15
data[c(5,15),]

# residuals vs. factor levels 

# by suppliers 
plot(as.numeric(data$supplier), residuals(anova_result),
     xlab = "supplier", ylab = "residuals", pch = 16)

# by sites 
plot(as.numeric(data$site), residuals(anova_result),
     xlab = "site", ylab = "residuals", pch = 16)

# Residuals vs. fitted
plot(fitted(anova_result), residuals(anova_result), xlab = "fitted", ylab = "residuals", pch = 16)


##########################################################################
# Gasoline Mixture Modification Example  
##########################################################################

data2 <- data.frame(drivers = as.factor(rep(1:4, each = 4)), cars = as.factor(rep(1:4, times = 4)),
                    treatment = c("A", "B", "D", "C", "D", "C", "A", "B", "B", "D", "C", "A", "C", "A", "B", "D"),  
                    response = c(19, 24, 23, 26, 23, 24, 19, 30, 15, 14, 15, 16, 19, 18, 19, 16))

# reproduce tables on slide 25
tapply(data2, INDEX = list(data2$drivers, data2$cars), function(z) z$treatment)
tapply(data2, INDEX = list(data2$drivers, data2$cars), function(z) z$response)

tapply(data2, INDEX = data2$drivers, FUN = function(z) mean(z$response))
tapply(data2, INDEX = data2$cars, FUN = function(z) mean(z$response))
tapply(data2, INDEX = data2$treatment, FUN = function(z) mean(z$response))

# data analysis 

# NOTE: this is not the preferred way of performing the analysis since we
# we have not specified proper contrasts. By default, R does not constrain 
# model coefficients for each term to sum to zero. Instead, the coefficient
# for the 'first' factor level (according to a lexicographical ordering) is
# set to zero. While this does not change the SS decomposition nor the F-test
# nor model fit, it changes interpretation of the coefficients significantly -- 
# they no longer represent differences from a grand mean but differences from
# the first factor level. This may make sense if the first factor level is a 
# designated reference such as 'placebo' (vs. different treatments)
# or 'established setup' of a manufactoring process (vs. innovated setup). In 
# general, however, there is no such designated reference.

aov0 <- aov(data = data2, formula = response ~ treatment + cars + drivers)
summary(aov0)
coef(aov0)

# sum-to-zero constraints
data3 <- data2
contrasts(data3$cars) <- contr.sum(4)
contrasts(data3$drivers) <- contr.sum(4)
data3$treatment <- as.factor(data3$treatment) # treatment is a character variable originally, so we convert it to a factor
contrasts(data3$treatment) <- contr.sum(4)

aov1 <- aov(data = data3, formula = response ~ treatment + cars + drivers)
summary(aov1)

coef(aov1)
# Interpretation of these coefficients:
# 
# intercept -- overall mean
mean(data2$response)
# treatment1 through treatment 3 -- difference from overall mean 
tapply(data2, INDEX = data2$treatment, FUN = function(z) mean(z$response)) - mean(data3$response)
# (similar for the blocking variables)

# Setting the contrast to 'contr.treatment' re-enforces the default behavior of R (as used for
# generating aov0). NOTE: contr.treatment is an R function. The fact that one of the variables
# in the analysis is called treatment is coincidental!

data4 <- data2
contrasts(data4$cars) <- contr.treatment(4)
contrasts(data4$drivers) <- contr.treatment(4)
data4$treatment <- as.factor(data4$treatment)
contrasts(data4$treatment) <- contr.treatment(4)

aov2 <- aov(data = data4, formula = response ~ treatment + cars + drivers)
summary(aov2)
coef(aov2)

# Interpretation of these coefficients:
# 
# intercept -- overall mean (as before!)

# treatment2 now represents the difference from what is seen under treatment1 (treatment = 'A'):
mean_reference <- tapply(data2, INDEX = data2$treatment, FUN = function(z) mean(z$response))[1]
tapply(data2, INDEX = data2$treatment, FUN = function(z) mean(z$response))[-1] -  mean_reference
#  (analogous for treatment3 and treatment4)

# checking that this coincides with the R default
mean(abs(coef(aov0) - coef(aov2)))

# the two choices contr.sum and contr.treatment lead to different design matrices, but the 
# range spaces of these matrix are equivalent:
X0 <- model.matrix(aov0)
X1 <- model.matrix(aov1)
# print(X0)
# print(X1)

# to check that the range spaces are the same
# tcrossprod(svd(X0)$u) - tcrossprod(svd(X1)$u) (these compares the two Hat [aka projection matrices] of the fit). 

# Accordingly, there is no difference in model fit: 
mean(abs(fitted(aov0) - fitted(aov1))
