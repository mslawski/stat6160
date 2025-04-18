##########################################################################
#
#
# R-code for Ch.~3b: One-way Analysis of Variance
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


# generate ANOVA table
anova_result <- aov(data = data, formula = etchrate ~ power_setting) 
print(anova_result)
summary(anova_result)

# generate individual contrasts for each mean comparison

# extract treatment effects
taus <- contrasts(data$power_setting) %*% coef(anova_result)[-1]

# construct matrix for evaluating pairwise differences 
C <- matrix(nrow = 4, ncol = choose(4, 2), data = 0)
l <- 0
for(i in 1:3){
    for(j in (i+1):4){
        l <- l + 1
        C[i,l] <- 1
        C[j,l] <- -1
    }
}

# contrasts
Gamma <- t(taus) %*% C

# associated t-statistic --- note: 5 observation in each group
Ts <- Gamma / sqrt(sigma(anova_result)^2 * (1/5 + 1/5))
# p-values (in percent)
k <- length(unique(power_setting))
pvals <- 2*(1-pt(abs(Ts), df = length(etchrate) - k)) * 100
# compute upper and lower confidence limits 
se <- sqrt(sigma(anova_result)^2 * (1/5 + 1/5))
crit <- qt(1-0.05/2, df = length(etchrate) - k)
lower <- Gamma - se * crit
upper <- Gamma + se * crit

lev <- levels(data$power_setting)
grp1 <- lev[apply(C, 2, function(z) which(z == -1))]
grp2 <- lev[apply(C, 2, function(z) which(z == 1))]

# create table 
data.frame(grp1, grp2, diff = as.numeric(Gamma), lower = as.numeric(lower), upper = as.numeric(upper),
           T = as.numeric(Ts), pvalue_perc = zapsmall(as.numeric(pvals))) 

# compute critical value when using Bonferroni correction
crit_bon <- qt(1-0.05/(2*choose(k,2)), df = length(etchrate) - k)
print(crit_bon)

# compute critical value when using Tukey's method 
qtukey(0.95, nmeans = k, df = length(etchrate) - k)/sqrt(2)


### Matrix representations 

# design matrix for a four-level factor, using a sum-to-zero constraint
# on the treatment effects.  
M <- matrix(nrow = 4, ncol = 4, data = c(1,1,0,0,
                                         1,0,1,0,
                                         1,0,0,1,
                                         1,-1,-1,-1), byrow = TRUE)
Minv <- solve(M) # returns matrix inverse
Minv %*% M # check if result is identity 
M %*% Minv

# returns linear model design matrix 
X <- model.matrix(anova_result)
betahat <- coef(anova_result)

# reproduce bethat using matrix algebra
y <- anova_result$model$etchrate
XtX <- crossprod(X)
Xty <- crossprod(X, y)
betahat_ma <- solve(XtX, Xty) # preferred over: solve(XtX) %*% Xty because of better numerical properties 
print(betahat)
print(betahat_ma) # identical. 

# To find the variance covariance matrix:
sigma(anova_result)^2 * solve(XtX) # according to formula on slide
zapsmall(vcov(anova_result)) # same result using the 'vcov' function in R
