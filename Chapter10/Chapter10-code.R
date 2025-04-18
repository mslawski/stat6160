##########################################################################
#
#
# R-code for Ch.~10: Nested and Split-Plot Designs
#
#
##########################################################################
ls()
rm(list = ls())
##########################################################################
# Two-Stage Nested Design: Suppliers and Batches 
##########################################################################

n <- 3 # number of replicates 
a <- 3  # number of suppliers
b <- 4 #  number of batches (within each supplier) 

N <- n * a * b #total number of measurements  
supplier <- rep(rep(1:a, each = b), each = n)
batches <- rep(rep(1:b, times = a), each = n)
response <- c(1, -1, 0, -2, -3, -4, -2, 0, 1, 1, 4, 0,
              1, -2, -3, 0, 4, 2, -1, 0, -2, 0, 3, 2,
              2, 4, 0, -2, 0, 2, 1, -1, 2, 3, 2, 1)

data <- data.frame(supplier = factor(supplier),
                   batches  = factor(batches),
                   response = response)

# ANOVA fit
# the forward slash "/" is used to express that batches are nested within supplier 
aov0 <- aov(response ~ supplier/batches, contrasts = list(supplier = "contr.sum",
                                                   batches = "contr.sum"), data = data)
summary(aov0)

# Note that the F-tests assume that both factors are fixed.

# If we assume that "batches" is a random factor while  
# "supplier" is fixed, we can keep the F-statistic
# from the above output, but need to change the  
# the F-statistic for the fixed factor as follows: 

aov0_sum <- unlist(summary(aov0))
FA <- aov0_sum["Mean Sq1"] / aov0_sum["Mean Sq2"]  
pA <- 1-pf(FA, df1 = aov0_sum["Df1"], df2 = aov0_sum["Df2"]) 
print(FA); print(pA)


##########################################################################
# Split-Plot Design 
##########################################################################


r <- 3 # 3 replicates
a <- 3 # 3 preparation methods
b <- 4 # 4 temperatures

PrepMethod <- rep(rep(c(1,2,3), each = b), times = r)
Temperature <- rep(rep(c(200, 225, 250, 275), times = 3), times = r)
Replicate <- rep(1:3, each = a * b) 
response <- c(30, 35, 37, 36, 34, 41, 38, 42, 29, 26, 33, 36,
              28, 32, 40, 41, 31, 36, 42, 40, 31, 30, 32, 40,
              31, 37, 41, 40, 35, 40, 39, 44, 32, 34, 39, 45)

data <- data.frame(PrepMethod = factor(PrepMethod),
                   Temperature = factor(Temperature),
                   Replicate = factor(Replicate),
                   response = response)

aov0 <- aov(response ~ PrepMethod*Temperature*Replicate,
            contrasts = list(PrepMethod = "contr.sum", Temperature = "contr.sum",
                             Replicate = "contr.sum"), data = data)

summary(aov0)

# Compute F-statistic and p-values, according to the table
# provided on the last slide:

aov0_sum <- unlist(summary(aov0))
MS <- aov0_sum[grep("Mean Sq", names(aov0_sum))]
DF <- aov0_sum[grep("Df", names(aov0_sum))]

# 1 --- PrepMethod (A)
# 3 --- Replicates
# 5 --- Replicates x A (Whole Plot Error) 
# 2 --- Temperature (B)
# 6 --- Replicates x B
# 4 --- AB interaction
# 7 --- Replicates x AB (Subplot error)

# F-test for A: MS_A / MS_{Replicates x A}

FA <- MS[1]/MS[5]
pA <- 1 - pf(FA, df1 = DF[1], df2 = DF[5])
print(FA); print(pA)

# F-test for B: MS_B / MS_{Replicates x B}
FB <- MS[2]/MS[6]
pB <- 1 - pf(FB, df1 = DF[1], df2 = DF[5])
print(FB); print(pB)

# F-test for AB: MS_AB / MS_{Replicates x AB}
FAB <- MS[4]/MS[7]
pAB <- 1 - pf(FAB, df1 = DF[4], df2 = DF[7])
print(FAB); print(pAB)

##########################################################################################################

