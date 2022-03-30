# ---------------------------------------------------------------------------- #
# ============================= PACKAGE VERSIONS ============================= #
# ---------------------------------------------------------------------------- #

library(mixOmics)

# ---------------------------------------------------------------------------- #

library(devtools)
setwd("C:/Users/Work/Desktop/mixOmics/R")

rm(list = c("network"))
devtools::load_all()

# ---------------------------------------------------------------------------- #
# =============================== LOAD IN DATA =============================== #
# ---------------------------------------------------------------------------- #

data("nutrimouse")
X <- nutrimouse$gene
Y <- nutrimouse$lipid

# ---------------------------------------------------------------------------- #
# ================================== TESTS =================================== #
# ---------------------------------------------------------------------------- #

pls.obj <- pls(X, Y)


setwd("C:/Users/Work/Desktop/mixOmics/R")
devtools::load_all()
#dev.new()

data("nutrimouse")
X <- nutrimouse$gene
Y <- nutrimouse$lipid

pls.obj <- pls(X, Y)

hist(1:10)
network.obj <- network(pls.obj, cutoff = 0.7, plot.graph = F)
plotIndiv(pls.obj)



network.obj <- network(pls.obj, cutoff = 0.75, save = 'jpeg', name.save = "asdf")













#-----------------------------------#
# construction of the initial graph #
#-----------------------------------#

# if (!plot.graph) {
#     ff <- tempfile()
#     png(filename=ff)
# }

nn = vcount(gR)
V(gR)$label.cex = min(2.5 * cex.node.name/log(nn), 1)
E(gR)$label.cex = min(2.25 * cex.edge.label/log(nn), 1)
cex0 = 2 * V(gR)$label.cex

def.par = par(no.readonly = TRUE)
dev.new()
par(pty = "s", mar = c(0, 0, 0, 0),mfrow=c(1,1))
plot(1:100, 1:100, type = "n", axes = FALSE, xlab = "", ylab = "")
cha = V(gR)$label
cha = paste("", cha, "")
xh = strwidth(cha, cex = cex0) * 1.5
yh = strheight(cha, cex = cex0) * 3

V(gR)$size = xh
V(gR)$size2 = yh

dev.off()

#return(NA)

if (is.null(layout.fun))
{
  l = layout.fruchterman.reingold(gR, weights = (1 - abs(E(gR)$weight)))
} else {
  l = layout.fun(gR)
}

if (isTRUE(!interactive))
{
  if (isTRUE(show.color.key))
  {
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    image(z.mat, col = col, xaxt = "n", yaxt = "n")
    box()
    par(usr = c(0, 1, 0, 1))
    axis(1, at = xv, labels = lv, cex.axis = keysize.label)
    title("Color key", font.main = 1, cex.main = keysize.label)
    par(def.par)
    par(new = TRUE)
  }
  
  par(pty = "s", mar = c(0, 0, 0, 0),mfrow=c(1,1))
  plot(gR, layout = l)
  par(def.par)
}

#-----------------------#
# procedure interactive #
#-----------------------#
gE.none = FALSE
if (isTRUE(interactive))
{
  
  # cutoff control bar #
  #-----------------------#
  min.cut = cutoff
  max.cut = max(abs(w))
  
  cutoff.old = cutoff
  
  dev.new("width" = 5, "height" = 2.7, xpos = -1)
  def.par = par(no.readonly = TRUE)
  
  cuts = seq(0, 1, length = 21)
  par(mai = c(0.25, 0.15, 0.3, 0.15), bg = gray(0.95))
  layout(matrix(c(0, 1, 0), ncol = 1, nrow = 3),
         widths = 1, heights = c(0.25, 1, 0.25))
  
  plot(cuts, type = "n", rep(0, 21), xlab = "", ylab = "",
       xlim = c(-0.10, 1.10), axes = FALSE)
  title("cutoff control", cex.main = 1.9, font.main = 1)
  text(0.5, -0.6, "value", cex = 1.5)
  text(0, -0.6, round(min.cut, 2), cex = 1.4)
  text(1, -0.6, round(max.cut, 2), cex = 1.4)
  mtext(min.cut, side = 1, line = -1, outer = FALSE, cex = 0.95)
  
  rect(-0.1, -0.3, -0.02, 0.3, col = "white")
  rect(1.02, -0.3, 1.1, 0.3, col = "white")
  points(1.06, 0, pch = 3, cex = 2.4)
  lines(c(-0.085, -0.035), c(0, 0))
  
  for (i in seq(0, 1, length = 21))
    lines(c(i, i), c(-0.22, 0.2))
  
  x = pos = 0
  rect(-0.01, -0.045, x, 0.04, col = "red")
  rect(x, -0.045, 1.01, 0.04, col = "white")
  
  bar.dev = dev.cur()
  
  # construction of the initial graph #
  #-----------------------------------#
  dev.new()
  net.dev = dev.cur()
  
  if (isTRUE(show.color.key))
  {
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    image(z.mat, col = col, xaxt = "n", yaxt = "n")
    box()
    par(usr = c(0, 1, 0, 1))
    axis(1, at = xv, labels = lv, cex.axis = keysize.label)
    title("Color key", font.main = 1, cex.main = keysize.label)
    par(def.par)
    par(new = TRUE)
  }
  
  par(pty = "s", mar = c(0, 0, 0, 0))
  plot(gR, layout = l)
  par(def.par)
  
  old.pos = -1
  
  repeat {
    dev.set(bar.dev)
    
    z = locator(1, type = "n")
    x = z[[1]]
    y = z[[2]]
    
    if (is.null(z)) break
    
    if (0 <= x & x <= 1 & -0.22 <= y & y <= 0.22)
    {
      rect(0, -0.045, x, 0.04, col = "red")
      rect(x, -0.045, 1.01, 0.04, col = "white")
      pos = x
    }
    
    if (1.02 <= x & x <= 1.1 & -0.3 <= y & y <= 0.3)
    {
      x = pos + 0.05
      idx = which.min(abs(cuts - x))
      x = cuts[idx]
      pos = x
      rect(0, -0.045, x, 0.04, col = "red")
      rect(x, -0.045, 1.01, 0.04, col = "white")
    }
    
    if (-0.1 <= x & x <= -0.02 & -0.3 <= y & y <= 0.3)
    {
      x = pos - 0.05
      idx = which.min(abs(cuts - x))
      x = cuts[idx]
      pos = x
      rect(0, -0.045, x, 0.04, col = "red")
      rect(x, -0.045, 1.01, 0.04, col = "white")
    }
    
    if (old.pos != pos)
    {
      old.pos = pos
      rect(0.4, -0.8, 0.6, -1.5, col = gray(0.95), border = NA)
      cutoff = (max.cut - min.cut) * pos + min.cut
      mtext(round(cutoff, 3), side = 1, line = -1, cex = 0.9)
      
      
      # new graph plot #
      #----------------#
      dev.set(net.dev)
      
      if (cutoff >= cutoff.old)
      {
        
        # selection of the edges to remove of the network #
        #-------------------------------------------------#
        supp.edge = E(gR)[abs(E(gR)$weight) < cutoff]
        
        # Generation of the graph with all the significant edges #
        #--------------------------------------------------------#
        gE = delete.edges(gR, supp.edge)
        gE = delete.vertices(gE, which(degree(gE) == 0))
        
        # graph plot #
        #------------#
        nn = vcount(gE)
        V(gR)$label.cex = min(2.5 * cex.node.name/log(nn), 1)
        E(gR)$label.cex = min(2.25 * cex.edge.label/log(nn), 1)
        cex0 = 2 * V(gE)$label.cex
        
        def.par = par(no.readonly = TRUE)
        
        par(pty = "s", mar = c(0, 0, 0, 0))
        plot(1:100, 1:100, type = "n", xaxt = "n")
        cha = V(gE)$label
        cha = paste("", cha, "")
        xh = strwidth(cha, cex = cex0) * 1.5
        yh = strheight(cha, cex = cex0) * 3
        
        V(gE)$size = xh
        V(gE)$size2 = yh
        
        par(def.par)
        
        if (is.null(layout.fun))
        {
          l = layout.fruchterman.reingold(gE, weights = (1 - abs(E(gE)$weight)))
        } else {
          l = layout.fun(gE)
        }
        
        if (isTRUE(show.color.key))
        {
          layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
          par(mar = c(5, 4, 2, 1), cex = 0.75)
          image(z.mat, col = col, xaxt = "n", yaxt = "n")
          box()
          par(usr = c(0, 1, 0, 1))
          axis(1, at = xv, labels = lv, cex.axis = keysize.label)
          title("Color key", font.main = 1, cex.main = keysize.label)
          par(def.par)
          par(new = TRUE)
        }
        
        par(pty = "s", mar = c(0, 0, 0, 0))
        plot(gE, layout = l)
        par(def.par)
        
        cutoff.old = cutoff
      } else {
        # selection of the edges to incluir in the network #
        #--------------------------------------------------#
        supp.edge = E(gR)[abs(E(gR)$weight) < cutoff]
        
        # generation of the graph with all the significant edges #
        #--------------------------------------------------------#
        gE = delete.edges(gR, supp.edge)
        gE = delete.vertices(gE, which(degree(gE) == 0))
        
        # graph plot #
        #------------#
        nn = vcount(gE)
        V(gR)$label.cex = min(2.5 * cex.node.name/log(nn), 1)
        E(gR)$label.cex = min(2.25 * cex.edge.label/log(nn), 1)
        cex0 = 2 * V(gE)$label.cex
        
        def.par = par(no.readonly = TRUE)
        
        par(pty = "s", mar = c(0, 0, 0, 0))
        plot(1:100, 1:100, type = "n", xaxt = "n")
        cha = V(gE)$label
        cha = paste("", cha, "")
        xh = strwidth(cha, cex = cex0) * 1.5
        yh = strheight(cha, cex = cex0) * 3
        
        V(gE)$size = xh
        V(gE)$size2 = yh
        
        par(def.par)	
        
        if (is.null(layout.fun))
        {
          l = layout.fruchterman.reingold(gE, weights = (1 - abs(E(gE)$weight)))
        } else {
          l = layout.fun(gE)
        }
        
        if (isTRUE(show.color.key))
        {
          layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
          par(mar = c(5, 4, 2, 1), cex = 0.75)
          image(z.mat, col = col, xaxt = "n", yaxt = "n")
          box()
          par(usr = c(0, 1, 0, 1))						
          axis(1, at = xv, labels = lv, cex.axis = keysize.label)
          title("Color key", font.main = 1, cex.main = keysize.label)
          par(def.par)
          par(new = TRUE)
        }
        
        par(pty = "s", mar = c(0, 0, 0, 0))
        plot(gE, layout = l)
        par(def.par)
        
        cutoff.old = cutoff
      }
      
      gE.none = TRUE
    }
    
  } # end loop
  
  if (gE.none != FALSE)
    gR = gE
}
res=list(gR = gR)


if(any(class.object %in% object.blocks))
{
  l = 1
  for (i in 1:(length(blocks)-1))
  {
    for (j in (i + 1):length(blocks))
    {
      M_block[[l]][abs(M_block[[l]]) < cutoff] = 0
      res[paste("M",blocks[i],blocks[j],sep="_")] = list(M_block[[l]])
      l = l + 1
    }
  }
} else {
  mat[abs(mat) < cutoff] = 0
  res$M=mat
}

res$cutoff = cutoff

if (!is.null(save))
  dev.off()

if (!plot.graph) {
  #unlink(ff)
  dev.off()
  dev.new
}

return(invisible(res))
}