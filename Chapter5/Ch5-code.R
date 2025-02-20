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
coef(anova_result)


# Non-orthogonality of treatment and blocks

X <- model.matrix(anova_result)

# we see that the columns corresponding to tire and compound are not orthogonal (their inner product is not zero):
crossprod(X)

# by comparison in case of a complete block design, treatments and block factors are orthogonal. As a result,
# the treatment effects and block effects can be estimated/tested independently.

# This is a useful package for generating various standard designs:  
library(AlgDesign)

design_complete <- gen.factorial(4, nVars=2, center=FALSE,  factors="all",varNames=c("compound", "tire"))
contrasts(design_complete$compound) <- contr.sum(4)
contrasts(design_complete$tire) <- contr.sum(4)
Xcomplete <- model.matrix(~compound + tire, data = design_complete)

crossprod(Xcomplete)

# we see that the the X' X matrix is block diagonal
# (and thus the covariance of the coefficients (X' X)^{-1} is block diagonal as well). 

				    
# !!! Because of lack of orthogonality, there is not a single canonical way of generating SS decompositions.
# For example, if we change the order of the terms in the LS decomposition, the SS decomposition becomes different:

anova_result2 <- aov(response ~ compound + tire, data = data) 
print(anova_result2) #! note the changes in the output
summary(anova_result2) #! note the changes in the output

# ASIDE -- by default, the ANOVA function in R computes sequential 
# sums of squares and associated test, i.e., term-by-term in the sequence they are entered.  
# For example, in the output of anova_result2, we compare
#
# SS(1) | SS(1 + compound)
# SS(1 + compound + tire) | SS(compound) 
#
# This is not of interest in our case --- we compare the reduction in 
# SS achieved by the variable "compound" relative to the intercept model), then the reduction in SS achieved by tire relative to the compound + intercept model. 
# Instead, we are interested in the opposite sequence 
# SS(1) | SS(1 + tire)
# SS(1 + tire + compound) | SS(tire)
#
# since "tire" is not the term of primary interest but a simply a term we would to account for.
#
# The use of the drop1() function is more symmetric, it compares the increase in SS after 
# dropping each term
drop1(anova_result, test = "F")
drop1(anova_result2, test = "F")

# Alternatively, we may use the 'Anova' function in the 'car' package.
library(car)
Anova(anova_result, type = 2) # type 2 means "type 2" (i.e., *partial* sum of squares)
Anova(anova_result, type = "II") # equivalent to previous
Anova(anova_result2, type = "II") # again, equivalent.

# Note that for balanced designs (such as the complete block design and LSD), 
# the sequence in which terms are added into the model does *NOT* matter:  
# Sequential (type1) SS and partial (type2) SS will coincide. For unbalanced 
# designs (such as BIBD), we need to be more careful since by default R 
# return type 1 SS, which may not be of interest (depending on the sequence 
# in which terms are entered). 

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

### [2] Incomplete Block Designs (using the 'Algdesign' package)

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
pairs <- expand.grid(1:t, 1:t)
pairs <- pairs[pairs[,1] > pairs[,2], ]

apply(pairs, 1, function(z) sum(rowSums(block_assignments_m[,c(z[1], z[2])]) == 2))
