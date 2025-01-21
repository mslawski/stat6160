##########################################################################
#
#
# R-code for Ch.~2: Simple Comparative Experiments
#
#
##########################################################################

##########################################################################
# Tomato Experiment Box, Hunter, Hunter (2005)
##########################################################################
fertilizer <- as.factor(c("A", "A", "B", "B", "A", "B", "B", "B", "A", "A", "B"))
weights <- c(29.9, 11.4, 26.6, 23.7, 25.3, 28.5, 14.2, 17.9, 16.5, 21.1, 24.3)

# use the t-test function in R
# NOTE: one-sided test H_1: \mu_B > \mu_A
t.test(x = weights[fertilizer == "B"], y = weights[fertilizer == "A"],
       alternative = "greater",
       var.equal = TRUE)

# extract two-sided confidence interval for difference of means 
t.test(x = weights[fertilizer == "B"], y = weights[fertilizer == "A"],
       alternative = "two.sided", var.equal = TRUE)$conf.int

# test equality of variances 

nA <- sum(fertilizer == "A")
nB <- sum(fertilizer == "B")

sA <- var(weights[fertilizer == "A"]) 
sB <- var(weights[fertilizer == "B"]) 

F <- sA/sB

# two-sided p-value (one-sided p-value coincides with the slide)
2 * (1 - pf(F, df1 = nA - 1, df2 = nB - 1))

# single-line command in R
var.test(weights[fertilizer == "A"], weights[fertilizer == "B"])

# assessing Normality visually (note that this is somewhat a stretch, given small sample sizes)
par(mfrow = c(2,1))
qqnorm((weights[fertilizer == "A"] - mean(weights[fertilizer == "A"]))/sqrt(sA), main = "Level A")
abline(0,1)


qqnorm((weights[fertilizer == "B"] - mean(weights[fertilizer == "B"]))/sqrt(sB), main = "Level B")
abline(0,1)


## Randomization Test 

# generate grid of all possible assignments of nA + nB samples (2^11 = 2048)
gr <- expand.grid(0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1)

# filter out all assignments with six ones in total (corresponding to nB)

gr_red <- gr[rowSums(gr) == nB,]

# number of admissible re-assignments
nrow(gr_red)

# we can also find this number by calculating a binomial coefficient 
choose(nA + nB, nB)

# obtain randomization distribution of difference in means
dist_rand <- apply(gr_red, 1, function(z) mean(weights[z == 1]) - mean(weights[z == 0])) 

par(mfrow = c(1,1))
hist(dist_rand, 13, main = "Histogram of randomization distribution",
     xlab = "statistic (mean difference)")

abline(v = mean(weights[fertilizer == "B"]) - mean(weights[fertilizer == "A"]), col = "red", lwd = 2) 

# obtain randomization p-value (one-sided) -- very consistent with the p-value of the t-test 
mean(dist_rand > (mean(weights[fertilizer == "B"]) - mean(weights[fertilizer == "A"])))


## Power calculations

# suppose that the mean difference between the two groups was two,
# what power did the experiment have in detecting this difference

# the pooled standard deviation is around 6
t.test(x = weights[fertilizer == "B"], y = weights[fertilizer == "A"],
       alternative = "greater",
       var.equal = TRUE)$stderr * 1/sqrt(1/nA + 1/nB)

power.t.test(n = min(nA, nB), delta = 2, sd = 6, sig.level = 0.05,
             type = "two.sample", alternative = "one.sided")

# with a power of .1227, the experiment was significantly underpowered 

# how many samples per group would we need to achieve a power of 80%: 

power.t.test(power = .8, delta = 2, sd = 6, sig.level = 0.05,
             type = "two.sample", alternative = "one.sided")


##########################################################################
# Paired Design: Shoe experiment 
##########################################################################

shoes <- data.frame(rbind(c("LR", 13.2, 14),
               c("LR", 8.2, 8.8),
               c("RL", 10.9, 11.2),
               c("LR", 14.3, 14.2),
               c("RL", 10.7, 11.8),
               c("LR", 6.6, 6.4),
               c("LR", 9.5, 9.8),
               c("LR", 10.8, 11.3),
               c("RL", 8.8, 9.3),
               c("LR", 13.3, 13.6)))
colnames(shoes) <- c("L/R", "A", "B")
shoes[,"A"] <- as.numeric(shoes[,"A"])
shoes[,"B"] <- as.numeric(shoes[,"B"])


# paired t-test 
t.test(x = shoes[,"B"], y = shoes[,"A"], paired = TRUE, alternative = "greater")

# corresponding two-sided confindence interval 
t.test(x = shoes[,"B"], y = shoes[,"A"], paired = TRUE, alternative = "two.sided")$conf.int

# this would be the two-sample test, ignoring paired design: 
t.test(x = shoes[,"B"], y = shoes[,"A"], paired = FALSE, alternative = "greater",
       var.equal = TRUE)

## linear model (last slide) 

# convert data set to long format 

library(tidyr)

shoes_long <- shoes %>% 
  pivot_longer(
    cols = c("A", "B"), 
    names_to = "material",
    values_to = "outcome")

shoes_long <- data.frame(shoes_long, person_ID = rep(1:10, each = 2))
shoes_long$material <- as.factor(shoes_long$material)
shoes_long$person_ID <- as.factor(shoes_long$person_ID)

contrasts(shoes_long$material) <-  contr.sum(2)
contrasts(shoes_long$person_ID) <- contr.sum(10)

summary(lm(outcome ~ material + person_ID, data = shoes_long))

# note that p-value for 'material' equals to the two-sided p-value from the paired test.  
