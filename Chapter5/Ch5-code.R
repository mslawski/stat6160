##########################################################################
#
#
# R-code for Ch.~5: Balanced Incomplete Block Design
#
#
##########################################################################
ls()
rm(list = ls())
##########################################################################
# Tire Experiment Data 
##########################################################################

data <- data.frame(tire = as.factor(rep(1:4, each = 3)), compound = as.factor(c("A", "B", "C", "A", "B", "D", "A", "C", "D", "B", "C", "D")),
                   response = c(238, 238, 279, 196, 213, 308, 254, 334, 367, 312, 421, 412))

# create table as on the slides: 

tapply(data, INDEX = list(data$tire, data$compound), FUN = function(z) z$response)


### enforce sum-to-zero constraints on the coefficients prior to running anova 

contrasts(data$tire) <- contr.sum(4)
contrasts(data$compound) <- contr.sum(4)


# generate ANOVA table
anova_result <- aov(response ~ tire + compound, data = data) 
print(anova_result)
summary(anova_result)

# align the above anova table with the expressions shown on the slides

Ti_ <- tapply(data, INDEX = list(data$compound), FUN = function(z) sum(z$response))
B_j <- tapply(data, INDEX = list(data$tire), FUN = function(z) sum(z$response))
k <- tapply(data, INDEX = list(data$tire), function(z) sum(table(z$tire)))[1] 

aux <- tapply(data, INDEX = list(data$compound), FUN = function(z) table(z$tire)) 
Q_i <- Ti_ -  (1/k) * unlist(lapply(aux, function(z) sum(z * B_j)))
lambda <- 2
t <- nlevels(data$compound)
tauhat_i <- k * Q_i / (lambda * t)    

k * sum(Q_i^2)/(lambda * t)

tau_i <- k * Q_i / (lambda * t)
print(tau_i)
coef(lm_res)


# Non-orthogonality of treatment and blocks

X <- model.matrix(anova_result)

# we see that the columns corresponding to tire and compound are not orthogonal (their inner product is not zero):
crossprod(X)

# by comparison in case of a complete block design, treatments and block factors are orthogonal. As a result,
# the treatment effects and block effects can be estimated/tested independently.

# !!! Because of lack of orthogonality, there is not a single canonical way of generating SS decompositions.
# For example, if we change the order of the terms in the LS decomposition, the SS decomposition becomes different:

anova_result2 <- aov(response ~ compound + tire, data = data) 
print(anova_result2) #! note the changes in the output
summary(anova_result2) #! note the changes in the output

# This is a useful package for generating various standard designs:  
library(AlgDesign)

design_complete <- gen.factorial(4, nVars=2, center=FALSE,  factors="all",varNames=c("compound", "tire"))
contrasts(design_complete$compound) <- contr.sum(4)
contrasts(design_complete$tire) <- contr.sum(4)
Xcomplete <- model.matrix(~compound + tire, data = design_complete)

crossprod(Xcomplete)

# we see that the the X' X matrix is block diagonal
# (and thus the covariance of the coefficients (X' X)^{-1} is block diagonal as well). 


### pair-wise comparisons:  

gr <- expand.grid(1:4, 1:4)
gr <- gr[gr[,1] > gr[,2],]
diffs <- tau_i[gr[,1]] - tau_i[gr[,2]] 
names(diffs) <- NULL

se <- sqrt(sigma(anova_result)^2 * (2*k)/(lambda * t))
Ts <- diffs / se
pvals <- 2*(1 - pt(abs(Ts), df = anova_result$df.residual))
pvals_format <- format(pvals, nsmall = 4, scientific = FALSE)
lev <- levels(data$compound)
labels <- paste(lev[gr[,1]], lev[gr[,2]], sep = " vs. ") 
outp <- data.frame(pair = labels, difference = diffs, 
			 stand_err = se, tstatistic = Ts, pval = pvals_format)   

print(outp)


### ASIDE --- functionality for generating specific designs:

library(agricolae)

### [1] Generate LSD's (Latin Square Designs)

p <- 4
ls4design <- design.lsd(letters[1:4], seed = 23)
print(ls4design)

### [2] Incomplete Block Designs 

# generate BIBD as in the tire example 
b <- 4
t <- 4
k <- 3 
BIB <- optBlock( ~ ., withinData = factor(1:t), blocksizes = rep(k, b))
block_assignments <- BIB$Blocks
block_assignments_m <- matrix(nrow = b, ncol = t, data = 0)

for(i in 1:length(block_assignments)){
block_assignments_m[i,unlist(block_assignments[[i]])] <- 1
}

print(block_assignments_m)

# check that this is a balanced designs in terms of group comparisons:

pairs <- expand.grid(1:t, 1:t)
pairs <- pairs[pairs[,1] > pairs[,2], ]

# every pair occurs exactly twice
apply(pairs, 1, function(z) sum(rowSums(block_assignments_m[,c(z[1], z[2])]) == 2))

    
# generate BIBD with 7 blocks, 7 treatments

b <- 7
t <- 7
k <- 3 
BIB <- optBlock( ~ ., withinData = factor(1:t), blocksizes = rep(k, b))
block_assignments <- BIB$Blocks
block_assignments_m <- matrix(nrow = b, ncol = t, data = 0)

for(i in 1:length(block_assignments)){
block_assignments_m[i,unlist(block_assignments[[i]])] <- 1
}

print(block_assignments_m)

# every pair occurs exactly once
apply(pairs, 1, function(z) sum(rowSums(block_assignments_m[,c(z[1], z[2])]) == 2))
