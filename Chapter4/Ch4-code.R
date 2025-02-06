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
