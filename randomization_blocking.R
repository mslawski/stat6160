#
# Simulation Study showing the benefits of randomization and blocking 
#

# t --  confounding variable that we want to block on 
t <- rep(1:5, each = 2)

# d --- treatment variable (binary): to be allocated below

# data-generating process depending on t and d
generate_response <- function(d, t){

beta0 <- 3 # intercept
alpha <- 2 # treatment effect
beta_t <- 1 # effect of confounder

X <- cbind(1, d, t)     
y <- X %*% c(beta0, alpha, beta_t)

y    
}    

### SETUP 1 --- we assign the first 5 subjects to treatment
d <- rep(c(1,0), each = 5)

# resulting "treatment effect" estimate 

y1 <- generate_response(d,t)

mean(y1[d == 1]) - mean(y1[d == 0]) 
# treatment effect estimate is heavily biased (-0.4 vs. 2[= truth]) 

### SETUP 2 --- full randomization 
    
nrepl <- 1E4 # we consider 10,000 randomizations, take avg.
est <- numeric(nrepl) 

for(i in 1:nrepl){
    d <- numeric(10)
    d[sample(1:10, 5, replace = FALSE)] <- 1
    y_i <- generate_response(d,t)
    est[i] <- mean(y_i[d == 1]) - mean(y_i[d == 0]) 

}

# mean estimator is unbiased, but variance of the estimator
# is huge. 

hist(est, nclass = 100, prob = T)
mean(est)
sd(est)


### SETUP 3 --- randomization within the following blocks
### for the variable t: [1 -- 2][3 -- 4][5]

blocks <- list(c(1,2),c(3,4),5)
block_sizes <- c(4,4,2)

nrepl <- 1E4 # we consider 10,000 randomizations, take avg.
est_block <- numeric(nrepl) 

for(i in 1:nrepl){
    d <- numeric(10)
    for(j in 1:length(blocks)){
        block_j <- which(t %in% blocks[[j]])
        s_j <- sample(block_j, block_sizes[j]/2, replace = FALSE) 
        d[s_j] <- 1
    }
    y_i <- generate_response(d,t)
    est_block[i] <- mean(y_i[d == 1]) - mean(y_i[d == 0]) 
}

mean(est_block)
sd(est_block)

# the variability of the estimator is almost reduced by
# a factor of 3!

### SETUP 3 --- randomization within the following blocks
### for the variable t: [1][2][3][4][5] (complete paired design).

blocks <- list(1,2,3,4,5)
block_sizes <- c(2,2,2,2,2)

nrepl <- 1E4 # we consider 10,000 randomizations, take avg.
est_block <- numeric(nrepl) 

for(i in 1:nrepl){
    d <- numeric(10)
    for(j in 1:length(blocks)){
        block_j <- which(t %in% blocks[[j]])
        s_j <- sample(block_j, block_sizes[j]/2, replace = FALSE) 
        d[s_j] <- 1
    }
    y_i <- generate_response(d,t)
    est_block[i] <- mean(y_i[d == 1]) - mean(y_i[d == 0]) 
}

mean(est_block)
sd(est_block)

# the variability of the estimator is zero!
