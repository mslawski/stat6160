##########################################################################
#
#
# R-code for Ch.~6: Factorial Designs for two/three factors 
#
#
##########################################################################
ls()
rm(list = ls())
##########################################################################
# Battery Example
##########################################################################

n = 4 # number of replications per treatment combination
I <- 3
J <- 3
data <- data.frame(type = as.factor(rep(1:3, each = J*n)),
                   temperature = as.factor(rep(c(rep(15, n), rep(70, n), rep(125, n)), times = I)),
                   response = c(130, 155, 74, 180, 34, 40, 80, 75, 20, 70, 82, 58,
                                150, 188, 159, 126, 136, 122, 106, 115, 25, 70, 58, 45,                                                                        138, 110, 168, 160, 174, 120, 150, 139, 96, 104, 82, 60)) 

# check data entry: 

tapply(data, INDEX = list(data$type, data$temperature), FUN = function(z) z$response[1])
tapply(data, INDEX = list(data$type, data$temperature), FUN = function(z) z$response[2])
tapply(data, INDEX = list(data$type, data$temperature), FUN = function(z) z$response[3])
tapply(data, INDEX = list(data$type, data$temperature), FUN = function(z) z$response[4])


### ANOVA table

# this one fits the plain main effect model 
anova_result0 <- aov(response ~ type + temperature, data = data,
                    contrasts = list(type = "contr.sum", temperature = "contr.sum"))

print(anova_result0)

# this one fits the interaction model: 
anova_result1 <- aov(response ~ type * temperature, data = data,
                    contrasts = list(type = "contr.sum", temperature = "contr.sum"))

print(anova_result1)
# note that the SS for the main effects remain unchanged: 
# this is because the corresponding columns in the design matrix are orthogonal: 
X1 <- model.matrix(anova_result1)
X1tX1 <- crossprod(X1) # computes the crossproduct t(X1) %*% X1 

### compare parameter estimates  

# intercept: 
muhat <- mean(data$response)

# main effects for variable 'type':
alphahat <- tapply(data, INDEX = data$type, function(z) mean(z$response)) - muhat 

# main effects for variable 'temperature': 
betahat <- tapply(data, INDEX = data$temperature, function(z) mean(z$response)) - muhat 

# interaction effects (here organized as a 3-by-3 matrix) 
alpha_beta_hat <- tapply(data, INDEX = list(data$type, data$temperature), function(z) mean(z$response)) - (outer(alphahat, betahat, FUN = "+")  + muhat)


# interaction plots   
interaction.plot(x.factor = data$temperature, trace.factor = data$type,
                 response = data$response, fixed = TRUE)

# Note that we can choose the other variable to be on the x-axis: 
interaction.plot(x.factor = data$type, trace.factor = data$temperature,
                 response = data$response, fixed = TRUE)

### The interaction plots confirm that the use of interaction terms
### is adequate. Note that the F-test also confirms the addition 
### of interaction terms.

### Finally, we can perform Tukey's non-additivity test:
library(agricolae)

nonadditivity(y = data$response, factor1 = data$type, factor2 = data$temperature,
              df = anova_result0$df.residual, MSerror = sigma(anova_result0)^2)

# Interestingly, the non-additive term is not significant. Note that this
# test is not equivalent to the F-test for the significance of the interaction terms.
# Tukey's non-additivity is based on a significantly constrained, 1-parameter model
# for the interactions. 

### Model Adequacy Checking:

# first: the main effect model 
plot(anova_result0, 0) # residual vs. fitted
plot(anova_result0, 1) # Q-Q plot
plot(as.numeric(data$type), anova_result0$residuals, xlab = "type") # residuals vs. factor 1 
plot(as.numeric(data$temperature), anova_result0$residuals, xlab = "temperature") # residuals vs. factor 2

# the residual vs. fitted plot for the main effect model shows
# a curvilinear up-and-down pattern. This in alignment with 
# the significance of the interaction terms. 

# now: the interaction model 
pdf("diagnostics_battery.pdf")
plot(anova_result1, 1) # residuals vs. fitted 
plot(anova_result1, 2) # Q-Q plot (slight departure from Normality) 
plot(as.numeric(data$type), anova_result1$residuals, xlab = "type") # residuals vs. factor1 
plot(as.numeric(data$temperature), anova_result1$residuals, xlab = "temperature") # residuals vs. factor2
dev.off()


##########################################################################
# Response Surfaces: Modeling temperature as a quadratic function 
##########################################################################

data2 <- data
data2$temperature <- as.numeric(levels(data$temperature)[data$temperature])

lmq <- lm(response ~ type * (temperature + I(temperature**2)), data = data,
          contrasts = list(type = "contr.sum", temperature = "contr.sum"))

print(anova_result1)
# note that the SS for the main effects remain unchanged: 
# this is because the corresponding columns in the design matrix are orthogonal: 
X1 <- model.matrix(anova_result1)
X1tX1 <- crossprod(X1) # computes the crossproduct t(X1) %*% X1 

### compare parameter estimates  

# intercept: 
muhat <- mean(data$response)

# main effects for variable 'type':
alphahat <- tapply(data, INDEX = data$type, function(z) mean(z$response)) - muhat 

# main effects for variable 'temperature': 
betahat <- tapply(data, INDEX = data$temperature, function(z) mean(z$response)) - muhat 

# interaction effects (here organized as a 3-by-3 matrix) 
alpha_beta_hat <- tapply(data, INDEX = list(data$type, data$temperature), function(z) mean(z$response)) - (outer(alphahat, betahat, FUN = "+")  + muhat)


# interaction plots   
interaction.plot(x.factor = data$temperature, trace.factor = data$type,
                 response = data$response, fixed = TRUE)

# Note that we can choose the other variable to be on the x-axis: 
interaction.plot(x.factor = data$type, trace.factor = data$temperature,
                 response = data$response, fixed = TRUE)

### The interaction plots confirm that the use of interaction terms
### is adequate. Note that the F-test also confirms the addition 
### of interaction terms.

### Finally, we can perform Tukey's non-additivity test:
library(agricolae)

nonadditivity(y = data$response, factor1 = data$type, factor2 = data$temperature,
              df = anova_result0$df.residual, MSerror = sigma(anova_result0)^2)

# Interestingly, the non-additive term is not significant. Note that this
# test is not equivalent to the F-test for the significance of the interaction terms.
# Tukey's non-additivity is based on a significantly constrained, 1-parameter model
# for the interactions. 

### Model Adequacy Checking:

# first: the main effect model 
plot(anova_result0, 0) # residual vs. fitted
plot(anova_result0, 1) # Q-Q plot
plot(as.numeric(data$type), anova_result0$residuals, xlab = "type") # residuals vs. factor 1 
plot(as.numeric(data$temperature), anova_result0$residuals, xlab = "temperature") # residuals vs. factor 2

# the residual vs. fitted plot for the main effect model shows
# a curvilinear up-and-down pattern. This in alignment with 
# the significance of the interaction terms. 

# now: the interaction model 
pdf("diagnostics_battery.pdf")
plot(anova_result1, 1) # residuals vs. fitted 
plot(anova_result1, 2) # Q-Q plot (slight departure from Normality) 
plot(as.numeric(data$type), anova_result1$residuals, xlab = "type") # residuals vs. factor1 
plot(as.numeric(data$temperature), anova_result1$residuals, xlab = "temperature") # residuals vs. factor2
dev.off()


##########################################################################
# Response Surfaces: Modeling temperature as a quadratic function 
##########################################################################

data2 <- data
data2$temperature <- as.numeric(levels(data$temperature)[data$temperature])

lmq <- lm(response ~ type * (temperature + I(temperature**2)), data = data2,
          contrasts = list(type = "contr.sum"))
anova(lmq)

sum(residuals(lmq)^2)
sum(residuals(anova_result1)^2)

### Residuals Sums of Squares are identical, and so are the number of parameters
### of the two models. An advantage of the quadratic model is increased interpretability
### via explicit response curves. 

temperature_grid <- seq(from = min(data2$temperature), to = max(data2$temperature), by = 1) 
newdata <- data.frame(type = factor(rep(1:3, times = length(temperature_grid))), temperature = rep(temperature_grid, each = 3))

predictions <- predict(lmq, newdata = newdata)

# plot response curves
for(i in 1:nlevels(newdata$type)){
    subs <- newdata$type == i
    if(i == 1){    
        plot(newdata$temperature[subs], predictions[subs], col = i, type = "l",
             xlab = "temperature", ylab = "predicted response", xlim = range(temperature_grid) * c(.97, 1.03),
             ylim = range(predictions) * c(.97, 1.03))
        points(data2$temperature[data$type == i], data2$response[data2$type == i], col = i, pch = 16) 
    }
    else{
        lines(newdata$temperature[subs], predictions[subs], col = i)
        points(data2$temperature[data$type == i], data2$response[data2$type == i], col = i, pch = 16) 
    }
}

legend("topright", col = 1:3, lwd = 1.8,
       legend = paste("type", 1:3, sep = " "), cex = 1.5)


