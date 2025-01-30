##########################################################################
#
#
# R-code demonstrating the use of Scheffe's method to safeguard 
# against false positive findings from ``data snooping" (i.e., adaptive
# selection of hypothesis to be tested -- here 'adaptive' means 
# that the hypothesis is selected *after* inspecting the data).
#
##########################################################################
ls()
##########################################################################

# Here, we simulate (normally distributed) data from k=5 treatment groups.
# There is *no difference* among the five groups, i.e., the data is 
# pure noise. After seeing the data, we identify the pair of groups
# with maximum difference, and then test whether that difference is
# statistically significant. We notice a significantly enhanced type I
# error.  
##########################################################################

n <- 10 # number of samples per treatment 
k <- 5  # number of treatments 
N <- 50 # total number of samples 
alpha <- .05 # significance level 

nrepl <- 1E4 # number of replications in the simulation
rejected_unadjusted <- numeric(nrepl)
rejected_tukey <- numeric(nrepl) 
rejected_scheffe <- numeric(nrepl)

for(r in 1:nrepl){

    y <- rnorm(N)
    group <- rep(1:k, each = n)
    group_means <- tapply(y, INDEX = group, FUN = mean)
    index_max <- which.max(group_means)
    index_min <- which.min(group_means)

    aov_res <- aov(y ~ group)

    # compute contrast associated with index_max and index_min
    C <- group_means[index_max] - group_means[index_min]
    S <- sqrt(sigma(aov_res)^2 *  (1/n + 1/n))
    F <- C^2 / S^2
    rejected_unadjusted[r] <- (F > qf(1-alpha, df1 = 1, df2 = N-k))

    # Tukey's studentized range statistic and Scheffe's method are expected
    # to safeguard against false discoveries from 'data snooping'.
    # Note that the t-statistic equals the root of the F-statistic.
    rejected_tukey[r] <- sqrt(F) > qtukey(1-alpha, nmeans = k, df = N-k)/sqrt(2)
    rejected_scheffe[r] <- F > ((k-1) * qf(1-alpha, df1 = k-1, N-k))
}

# without adjustment, the type I error is massively inflated. 
mean(rejected_unadjusted)
mean(rejected_tukey) # slightly below expected type I error rate
mean(rejected_scheffe) # much below expected type I error rate (conservative, but safe -- works regardless of the type of contrasts and numer of contrasts) 

