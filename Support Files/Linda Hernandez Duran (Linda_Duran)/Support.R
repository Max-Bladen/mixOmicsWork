library(mixOmics)
#library(devtools)
#setwd("~/GitHub/mixOmics/R")
#load_all()

# ============================================================================ #
# Set up data

# read in csv's
Xv<- read.csv("Val_sum_venomt.csv") # 30 x 98
Yb<- read.csv("Behaviours_val.csv") # 23 x 7

rownames(Xv) <- Xv$ID # set ID feature as rownames and remove ID feature
Xv <- Xv[,-1]

rownames(Yb) <- Yb$ID # set ID feature as rownames and remove ID feature
Yb <- Yb[,-1]

# identify which rows of Xv are also in Yb
# remove any rows from Xv that don't have equivalent sample in Yb
usable.idx <- which(rownames(Xv) %in% rownames(Yb))
Xv <- Xv[usable.idx, ]

# ============================================================================ #
# Using the ridge methodology without removing any features

# set the grid
grid <- seq(0.001,0.2,length=10)

# tune the rcc object to yield the optimal lambda values
cv.tunercc.val<-tune.rcc(Xv,Yb,
                         grid1 = grid,
                         grid2 = grid,
                         validation = "loo")

# extract these optimal lambda values
opt.l1_val<-cv.tunercc.val$opt.lambda1
opt.l2_val<-cv.tunercc.val$opt.lambda2


CV.rcc.valida.ridge <- rcc(Xv,Yb, method = "ridge",
                           lambda1 = opt.l1_val, lambda2 = opt.l2_val)
### produces no errors

plot(CV.rcc.valida.ridge, type = "barplot") 
### produces no errors

plotVar(CV.rcc.valida.ridge)
### PRODUCES THE FOLLOWING WARNINGS (but still plots):
# Warning messages:
# 1: In cor(object$X, object$variates$X[, c(comp1, comp2)] + object$variates$Y[,  : 
#       the standard deviation is zero
# 2: Removed 14 rows containing missing values (geom_point). 
# 3: Removed 14 rows containing missing values (geom_text).

# ============================================================================ #
# Using the shrinkage methodology without removing any features

CV.rcc.valida.shrink <- rcc(Xv,Yb, method = "shrinkage")
### PRODUCES THE FOLLOWING WARNINGS:
# Warning messages:
#   1: 14 instances of variables with zero scale detected! 
#   2: 14 instances of variables with zero scale detected! 
#   3: 14 instances of variables with zero scale detected! 

plot(CV.rcc.valida.shrink, type = "barplot")
### produces no errors

plotVar(CV.rcc.valida.shrink)
### PRODUCES THE FOLLOWING WARNINGS (but still plots):
# Warning messages:
# 1: In cor(object$X, object$variates$X[, c(comp1, comp2)] + object$variates$Y[,  : 
#       the standard deviation is zero
# 2: Removed 14 rows containing missing values (geom_point). 
# 3: Removed 14 rows containing missing values (geom_text).


# ============================================================================ #
# Let's now try removing rows which have all 0s or all 1s, meaning that they
# have no variation (and hence no standard deviation)

# which features have zero variance (all values are the same)
zero.var.feats <- as.vector(which(apply(Xv, 2, sd) == 0))
Xv <- Xv[, -zero.var.feats] # remove them from the Xv dataframe

# ============================================================================ #
# Using the ridge methodology AFTER removing zero-variance features

CV.rcc.valida2.ridge <- rcc(Xv,Yb, method = "ridge",
                            lambda1 = opt.l1_val, lambda2 = opt.l2_val)
### produces no errors

plot(CV.rcc.valida2.ridge, type = "barplot")
### produces no errors

plotVar(CV.rcc.valida2.ridge)
### produces no errors

# ============================================================================ #
# Using the shrinkage methodology AFTER removing zero-variance features

CV.rcc.valida2.shrink <- rcc(Yb,Xv, method = 'shrinkage')
### produces no errors

plot(CV.rcc.valida2.shrink, type = "barplot")
### produces no errors

plotVar(CV.rcc.valida2.shrink)
### produces no errors

plotVar(CV.rcc.valida2.shrink, var.names = c(T,F),
        cex = c(4, 4), cutoff = 0.5,
        title = "(b) H. valida, rCCA shrinkage comp 1 - 2")
### produces no errors

#cim(CV.rcc.valida2.shrink, comp = 1:2, xlab = "Venom masses", ylab = "behaviours")


