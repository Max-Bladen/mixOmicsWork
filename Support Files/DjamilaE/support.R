suppressMessages(library(mixOmics))

data("liver.toxicity")

X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
group <- liver.toxicity$treatment$Dose.Group

# taken from the sPLS Case Study (mixOmics.org)
optimal.ncomp <- 2
optimal.keepX <- c(35, 45)
optimal.keepY <- c(4,4)

liver.spls <- spls(X, Y, mode = "canonical",
                   ncomp = optimal.ncomp,
                   keepX = optimal.keepX,
                   keepY = optimal.keepY)

par(mfrow=c(2,2))
gene.splsda <- splsda(X, group,
                      ncomp = optimal.ncomp,
                      keepX = optimal.keepX)
plotLoadings(gene.splsda, contrib = "max", method = "median",
             title = "Figure 1a, comp1")
plotLoadings(gene.splsda, contrib = "max", method = "median", comp = 2,
             title = "Figure 1b, comp2") #





par(mfrow=c(2,2))
treatment.splsda <- splsda(Y, group,
                           ncomp = optimal.ncomp,
                           keepX = optimal.keepY)
plotLoadings(treatment.splsda, contrib = "max", method = "median",
             title = "Figure 2a, comp1")
plotLoadings(treatment.splsda, contrib = "max", method = "median", comp = 2,
             title = "Figure 2b, comp2") #



selected.genes <- rownames(which(gene.splsda$loadings$X!=0, arr.ind = T))
selected.treaments <- rownames(which(treatment.splsda$loadings$X!=0, arr.ind = T))

X.s <- X[, selected.genes]
Y.s <- Y[, selected.treaments]

heatmap(cor(X.s, Y.s))

subset.spls <- spls(X.s, Y.s, mode = "canonical", ncomp = 5)
sub.spls.perf <- perf(subset.spls, folds = 5, nrepeat = 5)
plot(sub.spls.perf, criterion = "cor.tpred")

sub.spls.tune <- tune.spls(X, Y, test.keepX = c(1:10), folds = 5, ncomp = 5)
plot(sub.spls.tune)
