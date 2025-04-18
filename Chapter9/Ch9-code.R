##########################################################################
#
#
# R-code for Ch.~9: Analysis of Covariance (ANCOVA)
#
#
##########################################################################
ls()
rm(list = ls())
##########################################################################
# Blood Pressure Data
##########################################################################

n <- 10 # number
t <- 2  # number of treatments
N <- n*t # replicates x t
data <- data.frame(x = c(135, 125, 125, 130, 105, 130, 140, 93, 110, 100, 
                         90, 135, 130, 115, 110, 140, 130, 95, 90, 105),
                   y = c(45, 45, 20, 50, 25, 37, 50, 20, 25, 15,
                         34, 55, 50, 45, 30, 45, 45, 23, 40, 35),
                   treatment = as.factor(rep(c("A", "B"), each = n))) 



# standard ANOVA:
aov0 <- aov(y ~ treatment, contrasts = list(treatment = "contr.sum"), data = data)
summary(aov0)

# Normal Q-Q plot
plot(aov0, which = 2)

# scatter-plot Y vs. X

plot(data$x, data$y, xlab = "x", ylab = "y",
    col = ifelse(data$treatment == "A", "red", "blue"),
    pch = 16)

# fit ANCOVA model II
ancova2 <- lm(y ~ I(x - mean(x)) + treatment, contrasts = list(treatment = "contr.sum"), data = data) 
summary(ancova2)
anova(ancova2)

# ANCOVA model I
ancova1 <- lm(y ~ x + treatment, contrasts = list(treatment = "contr.sum"), data = data) 
summary(ancova1)
anova(ancova1)

# visualizing the least squares fit:

xgrid <- seq(from = min(data$x), to = max(data$x), by = 0.1)
pred_grid <- predict(ancova1, newdata = data.frame(x = c(xgrid, xgrid), treatment = rep(c("A", "B"), times = length(xgrid))))

plot(c(xgrid, xgrid), pred_grid, cex = 0.25, col = rep(c("red", "blue"), times = length(xgrid)),
     xlab = "x", ylab = expression(hat(y)), ylim = range(data$y) * c(.98, 1.02))
points(data$x, data$y, col = ifelse(data$treatment == "A", "red", "blue"),
    pch = 16) 

abline(coef(lm(y ~ x, data = data)), lwd = 2, lty = "dashed")
text(105, 55, labels = "dashed line:\n least squares regression \n            line assuming no treatment effect")


### Error decomposition, AOV table calculations

attach(data)
SSy <- sum((y - mean(y))^2)
SSx <- sum((x - mean(x))^2)           
SP_xy <- sum((y - mean(y))*(x - mean(x)))

SS_covar <- SP_xy^2 / SSx

SSy_star <- sum(tapply(y, INDEX = treatment, function(z) sum((z - mean(z))^2)))
SSx_star <- sum(tapply(x, INDEX = treatment, function(z) sum((z - mean(z))^2)))
SSxy_star <- sum(unlist(tapply(x, INDEX = treatment, function(z) z - mean(z))) *
             unlist(tapply(y, INDEX = treatment, function(z) z - mean(z)))) 
SS_Tradj <- (SSy - SS_covar) - (SSy_star - SSxy_star^2 / SSx_star)

SSE <- SSy_star - SSxy_star^2 / SSx_star

print(SS_covar); print(SS_Tradj); print(SSE)
anova(ancova2)

### ANCOVA model III -- baseline constraint 

ancova3 <- lm(y ~ x + treatment, contrasts = list(treatment = "contr.treatment"), data = data) 
summary(ancova3)
anova(ancova3)

### ANCOVA model IV -- models with different regression slopes (aka interaction model)

xA <- x * (treatment == "A")
xB <- x * (treatment == "B")

#
ancova4 <- lm(y ~ xA + xB  + treatment, contrasts = list(treatment = "contr.treatment"), data = data) 
summary(ancova4)
anova(ancova4)

### ANCOVA model V -- "proper" interaction model with referece coding 
ancova5 <- lm(y ~ x*treatment, contrasts = list(treatment = "contr.treatment"), data = data) 
summary(ancova5)

### Visualization

pred_grid_int <- predict(ancova5, newdata = data.frame(x = c(xgrid, xgrid), treatment = rep(c("A", "B"), times = length(xgrid))))

plot(c(xgrid, xgrid), pred_grid_int, cex = 0.25, col = rep(c("red", "blue"), times = length(xgrid)),
     xlab = "x", ylab = expression(hat(y)), ylim = range(data$y) * c(.98, 1.02))
points(data$x, data$y, col = ifelse(data$treatment == "A", "red", "blue"),
    pch = 16) 

##########################################################################################################
