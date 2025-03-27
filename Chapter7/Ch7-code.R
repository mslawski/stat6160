##########################################################################
#
#
# R-code for Ch.~7: 2^k factorial designs
#
#
##########################################################################
ls()
rm(list = ls())
##########################################################################
# Spring Experiment
##########################################################################
k <- 3
N <- 2^k
data <- matrix(nrow = N, ncol = k+1)
data[,1] <- rep(c(70, 120), each = 2^(k-1))
data[,2] <- rep(rep(c(0.5, 0.7), each = 2), times = 2)
data[,3] <- rep(c(1450, 1600), times = 2^(k-1))
data[,4] <- c(67, 79, 61, 75, 59, 90, 52, 87)
colnames(data) <- c("OT", "CP", "ST", "Y")
data <- as.data.frame(data)
data[, "OT"] <- as.factor(data[, "OT"])
data[,"CP"] <- as.factor(data[, "CP"])
data[,"ST"] <- as.factor(data[, "ST"])

# design matrix
X <- model.matrix(~OT*CP*ST, data = data,
                  contrasts = list(OT = "contr.sum", CP = "contr.sum", ST = "contr.sum"))

# note that the columns for the interaction columns are obtained by
# multiplying the columns of the main effect
all.equal(X[,"OT1"] * X[,"CP1"], X[,"OT1:CP1"])
all.equal(X[,"CP1"] * X[,"ST1"], X[,"CP1:ST1"])
all.equal(X[,"OT1"] * X[,"ST1"], X[,"OT1:ST1"])
all.equal(X[,"OT1"] * X[,"CP1"] * X[,"ST1"], X[,"OT1:CP1:ST1"])

# check orthogonality:
crossprod(X)

# effects:
aov1 <- aov(Y ~ OT*CP*ST, data = data,
            contrasts =list(OT = "contr.sum", CP = "contr.sum", ST = "contr.sum"))

effects <- coef(aov1)[-1] * 2 # [-1] for dropping the intercept
model.tables(aov1, type = "effects", se = F) # built-in R function

all_means <- model.tables(aov1, type = "means", se = F)

# main effects plot: 
plot(c(1,2), all_means$tables$OT, type = "b", ylim = c(55, 85), xlim = c(0.5, 9.5),
     cex.axis = 0.01, ylab = "mean difference", xlab = "", pch = 16 ,cex.lab = 1.7) 
points(c(4,5), all_means$tables$CP, type = "b", ylim = c(55, 85), xlim = c(0.5, 9.5),
     cex.axis = 0.01, ylab = "mean difference", xlab = "", pch = 16) 
points(c(7,8), all_means$tables$ST, type = "b", ylim = c(55, 85), xlim = c(0.5, 9.5),
       cex.axis = 0.01, ylab = "mean difference", xlab = "", pch = 16)
text(1.5, 75, "OT", cex = 1.5)
text(4.5, 75, "CP", cex = 1.5)
text(7.5, 75, "ST", cex = 1.5)

# interaction plot (manual):
plot(c(1,2), all_means$tables["OT:CP"][[1]][1,], type = "b", ylim = c(50, 95), xlim = c(0.5, 9.5),
     cex.axis = 0.01, ylab = "mean difference", xlab = "", pch = 16 ,cex.lab = 1.7) 
points(c(1,2), all_means$tables["OT:CP"][[1]][2,], type = "b",
       cex.axis = 0.01, ylab = "mean difference", xlab = "", pch = 16 ,cex.lab = 1.7, col = "blue")


points(c(4,5), all_means$tables["CP:ST"][[1]][1,], type = "b", 
       cex.axis = 0.01, ylab = "mean difference", xlab = "", pch = 16 ,cex.lab = 1.7)
points(c(4,5), all_means$tables["CP:ST"][[1]][2,], type = "b", 
     cex.axis = 0.01, ylab = "mean difference", xlab = "", pch = 16 ,cex.lab = 1.7, col = "blue") 


points(c(7,8), all_means$tables["OT:ST"][[1]][1,], type = "b", 
       cex.axis = 0.01, ylab = "mean difference", xlab = "", pch = 16 ,cex.lab = 1.7)
points(c(7,8), all_means$tables["OT:ST"][[1]][2,], type = "b", 
     cex.axis = 0.01, ylab = "mean difference", xlab = "", pch = 16 ,cex.lab = 1.7, col = "blue") 

text(1.5, 80, "OT:CP", cex = 1.5)
text(4.5, 80, "CP:ST", cex = 1.5)
text(7.5, 80, "OT:ST", cex = 1.5)

# Normal Plot and Half-Normal Plot

# Normal Plot
ord <- order(coef(aov1)[-1])
plot(qnorm((1:(N-1))/N), coef(aov1)[-1][ord], pch = 16,
     xlab = "Normal Quantiles", ylab = "Effects", cex.lab = 1.7, ylim = c(-12,8), xlim = c(-1.5, 1.5))
text(qnorm((1:(N-1))/N), coef(aov1)[-1][ord] + 0.5,
     labels = names(coef(aov1)[-1])[ord])

# Half-Normal Plot 

ord <- order(abs(coef(aov1)[-1]))
plot(qnorm(((N - 1) + (1:(N-1)))/(2*(N-1) + 1)), abs(coef(aov1)[-1])[ord], pch = 16,
     xlab = "Half-Normal Quantiles", ylab = "Effects", cex.lab = 1.7, ylim = c(0, 12), xlim = c(0, 1.5))

text(qnorm(((N - 1) + (1:(N-1)))/(2*(N-1) + 1)), abs(coef(aov1)[-1][ord]) + 0.5,
     labels = names(coef(aov1)[-1])[ord])

# How would go about finding your final model (?)
