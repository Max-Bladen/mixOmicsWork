library(mixOmics)

data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
nutrimouse.sgccda <- wrapper.sgccda(X=data,
                                    Y = Y,
                                    design = design,
                                    keepX = list(gene=c(10,10), lipid=c(15,15)),
                                    ncomp = 2,
                                    scheme = "horst")
circosPlot(nutrimouse.sgccda, cutoff = 0.7)