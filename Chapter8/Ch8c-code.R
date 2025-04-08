##########################################################################
#
#
# R-code for Ch.~8: Random Effects and Mixed Effects Models
#                   for factor experiments with three or more factors,
#                   Sattherswaite approximation
#
##########################################################################
ls()
rm(list = ls())
##########################################################################
# Pressure Drop Data set 
##########################################################################

data0 <- read.table("PressureDrop.csv", header = TRUE, sep = ",")
data <- data0
factornames <- c("Temp", "Operator", "Gauge")
for(f in factornames)
    data[,f] <- as.factor(data[,f])

aov0 <- aov(Drop ~ Temp*Operator*Gauge, data = data,
            contrasts = list(Temp= "contr.sum",
                             Operator = "contr.sum", Gauge = "contr.sum"))

summary(aov0)


# calculate F-statistics for each of the terms 
aov0_sum <- unlist(summary(aov0))

MSEs <- aov0_sum[grep("Mean Sq", names(aov0_sum))]
DFs <- aov0_sum[grep("Df", names(aov0_sum))]

# [1] - A, [2] - B, [3] - C
# [4] - AB, [5] - AC, [6] - BC, [7] - ABC, [8] - E.

# (Operator)
F_B <- MSEs[2] / MSEs[6]
p_B <- 1 - pf(F_B, df1 = DFs[2], df2 = DFs[6])

F_C <- MSEs[3] / MSEs[6]
p_C <- 1 - pf(F_C, df1 = DFs[3], df2 = DFs[6])

F_AB <- MSEs[4]/MSEs[7] 
p_AB <- 1 - pf(F_AB, df1 = DFs[4], df2 = DFs[7])
        
F_AC <- MSEs[5]/MSEs[7] 
p_AC <- 1 - pf(F_AC, df1 = DFs[5], df2 = DFs[7])

### Note that denominator for this interaction term becomes MS_E! 
F_BC <- MSEs[6]/MSEs[8] 
p_BC <- 1 - pf(F_BC, df1 = DFs[6], df2 = DFs[8])

F_ABC <- MSEs[7]/MSEs[8] 
p_ABC <- 1 - pf(F_ABC, df1 = DFs[7], df2 = DFs[8])
        
# approximate F-statistic (Sattherswaite approximation)
# for the A main effect

MSnum <- MSEs[1] + MSEs[7]
MSdenom <- MSEs[4] + MSEs[5]

F_A <- MSnum / MSdenom 

DFnum <- MSnum^2 / (MSEs[1]^2/DFs[1] + MSEs[7]^2/DFs[7])
DFdenom <- MSdenom^2 / (MSEs[4]^2/DFs[4] + MSEs[5]^2/DFs[5])    

# here, we simply round the DFs to the nearest integer:
p_A <- 1 - pf(F_A, df1 = DFnum, df2 = DFdenom)
