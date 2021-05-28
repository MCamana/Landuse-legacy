plotTaxa2<-function (titan.out, z1 = T, z2 = T, interval = T, prob95 = F, 
                     z.med = F, xlabel = "Environmental Gradient", log = "", at = NULL, 
                     xmin = min(titan.out$sppmax[, 8]), xmax = max(titan.out$sppmax[, 
                                                                                    12]) * 1.05, tck = 0.025, bty = "u", ntick = 6, prtty = T, 
                     dig = 1, leg.x = 0.8, leg.y = 10, cex.taxa = 0.75, cex = 1.25, 
                     cex.axis = 1.25, cex.leg = 1.25, cex.lab = 1.25, legend = T, 
                     col1 = "black", fil1 = "black", col2 = "black", fil2 = "white", 
                     write = F, all = F, ...) 
{
  imax = titan.out$arguments[[5]]
  boot = titan.out$arguments[[3]] > 0.5
  if (all) {
    boot = F
  }
  if (boot) {
    if (z1) {
      sppsub1 <- subset(titan.out$sppmax, titan.out$sppmax[, 
                                                           16] == 1)
    }
    if (z2) {
      sppsub2 <- subset(titan.out$sppmax, titan.out$sppmax[, 
                                                           16] == 2)
    }
  }
  else {
    if (z1) {
      sppsub1 <- subset(titan.out$sppmax, titan.out$sppmax[, 
                                                           4] == 1 & titan.out$sppmax[, 6] <= 0.05)
    }
    if (z2) {
      sppsub2 <- subset(titan.out$sppmax, titan.out$sppmax[, 
                                                           4] == 2 & titan.out$sppmax[, 6] <= 0.05)
    }
  }
  if (z1) {
    if (nrow(sppsub1) < 1) {
      stop("z1 is empty, set z1=FALSE, change significance criteria, or set boot=FALSE if bootstrapping was not used to generate the titan object")
    }
  }
  if (z2) {
    if (nrow(sppsub2) < 1) {
      stop("z2 is empty, set z2=FALSE, change significance criteria, or set boot=FALSE if bootstrapping was not used to generate the titan object")
    }
  }
  #par(mar = c(8, 8, 1, 8), oma = c(0, 3, 0, 3))
  if (z1 & z2) {
    if (nrow(sppsub1) >= nrow(sppsub2)) {
      sppsub.gt <- sppsub1
    }
    else {
      sppsub.gt <- sppsub2
    }
  }
  else {
    if (z1) {
      sppsub.gt <- sppsub1
    }
    else {
      sppsub.gt <- sppsub2
    }
  }
  plot(sppsub.gt[, 1], ((max(rank((sppsub.gt[, 1]), ties.method = "first")) + 
                           1) - rank((sppsub.gt[, 1]), ties.method = "first")), 
       xlim = c(xmin, xmax), ylim = c(0.5, max(rank(sppsub.gt[, 
                                                              1], ties.method = "first") + 1)), cex = 0, tck = tck, 
       log = log, axes = FALSE, ylab = "", xlab = "")
  if (boot) {
    if (prob95) {
      if (z1) {
        yvalues1 = ((max(rank((sppsub1[, 12]), ties.method = "first")) + 
                       1) - rank((sppsub1[, 12]), ties.method = "first"))
      }
      if (z2) {
        yvalues2 = rank((sppsub2[, 8]), ties.method = "first") + 
          0.5
      }
    }
    else {
      if (z.med) {
        if (z1) {
          yvalues1 = ((max(rank((sppsub1[, 10]), ties.method = "first")) + 
                         1) - rank((sppsub1[, 10]), ties.method = "first"))
        }
        if (z2) {
          yvalues2 = rank((sppsub2[, 10]), ties.method = "first") + 
            0.5
        }
      }
      else {
        if (imax) {
          if (z1) {
            yvalues1 = ((max(rank((sppsub1[, 1]), ties.method = "first")) + 
                           1) - rank((sppsub1[, 1]), ties.method = "first"))
          }
          if (z2) {
            yvalues2 = rank((sppsub2[, 1]), ties.method = "first") + 
              0.5
          }
        }
        else {
          if (z1) {
            yvalues1 = ((max(rank((sppsub1[, 2]), ties.method = "first")) + 
                           1) - rank((sppsub1[, 2]), ties.method = "first"))
          }
          if (z2) {
            yvalues2 = rank((sppsub2[, 2]), ties.method = "first") + 
              0.5
          }
        }
      }
    }
  }
  else {
    if (imax) {
      if (z1) {
        yvalues1 = ((max(rank((sppsub1[, 1]), ties.method = "first")) + 
                       1) - rank((sppsub1[, 1]), ties.method = "first"))
      }
      if (z2) {
        yvalues2 = rank((sppsub2[, 1]), ties.method = "first") + 
          0.5
      }
    }
    else {
      if (z1) {
        yvalues1 = ((max(rank((sppsub1[, 2]), ties.method = "first")) + 
                       1) - rank((sppsub1[, 2]), ties.method = "first"))
      }
      if (z2) {
        yvalues2 = rank((sppsub2[, 2]), ties.method = "first") + 
          0.5
      }
    }
  }
  if (boot & interval) {
    if (z1) {
      segments(sppsub1[, 8], yvalues1, sppsub1[, 12], yvalues1, 
               col = col1, lwd = 2)
    }
    if (z2) {
      segments(sppsub2[, 8], yvalues2, sppsub2[, 12], yvalues2, 
               col = col2, lwd = 2, lty = 3)
    }
  }
  if (z1) {
    grpcol = rep(NA, nrow(sppsub1))
  }
  if (z2) {
    grpcol2 = rep(NA, nrow(sppsub2))
  }
  if (z1) {
    for (i in 1:nrow(sppsub1)) {
      if (sppsub1[i, 4] > 1.5) {
        grpcol[i] = col2
      }
      else {
        grpcol[i] = fil1
      }
    }
  }
  if (z2) {
    for (i in 1:nrow(sppsub2)) {
      if (sppsub2[i, 4] > 1.5) {
        grpcol2[i] = fil2
      }
      else {
        grpcol2[i] = fil1
      }
    }
  }
  if (boot) {
    if (prob95) {
      if (z.med) {
        if (z1) {
          symbols(sppsub1[, 12], yvalues1, circles = sppsub1[, 
                                                             15], inches = 0.1, add = TRUE, xlim = c(0, 
                                                                                                     5), fg = col1, bg = grpcol, lwd = 2)
        }
        if (z2) {
          symbols(sppsub2[, 8], yvalues2, circles = sppsub2[, 
                                                            15], inches = 0.1, add = TRUE, xlim = c(0, 
                                                                                                    5), fg = col2, bg = grpcol2, lwd = 2)
        }
      }
      else {
        if (z1) {
          symbols(sppsub1[, 12], yvalues1, circles = sppsub1[, 
                                                             7], inches = 0.1, add = TRUE, xlim = c(0, 
                                                                                                    5), fg = col1, bg = grpcol, lwd = 2)
        }
        if (z2) {
          symbols(sppsub2[, 8], yvalues2, circles = sppsub2[, 
                                                            7], inches = 0.1, add = TRUE, xlim = c(0, 
                                                                                                   5), fg = col2, bg = grpcol2, lwd = 2)
        }
      }
    }
    else {
      if (z.med) {
        if (z1) {
          symbols(sppsub1[, 10], yvalues1, circles = sppsub1[, 
                                                             15], inches = 0.1, add = TRUE, xlim = c(0, 
                                                                                                     5), fg = col1, bg = grpcol, lwd = 2)
        }
        if (z2) {
          symbols(sppsub2[, 10], yvalues2, circles = sppsub2[, 
                                                             15], inches = 0.1, add = TRUE, xlim = c(0, 
                                                                                                     5), fg = col2, bg = grpcol2, lwd = 2)
        }
      }
      else {
        if (imax) {
          if (z1) {
            symbols(sppsub1[, 1], yvalues1, circles = sppsub1[, 
                                                              7], inches = 0.1, add = TRUE, xlim = c(0, 
                                                                                                     5), fg = col1, bg = grpcol, lwd = 2)
          }
          if (z2) {
            symbols(sppsub2[, 1], yvalues2, circles = sppsub2[, 
                                                              7], inches = 0.1, add = TRUE, xlim = c(0, 
                                                                                                     5), fg = col2, bg = grpcol2, lwd = 2)
          }
        }
        else {
          if (z1) {
            symbols(sppsub1[, 2], yvalues1, circles = sppsub1[, 
                                                              7], inches = 0.1, add = TRUE, xlim = c(0, 
                                                                                                     5), fg = col1, bg = grpcol, lwd = 2)
          }
          if (z2) {
            symbols(sppsub2[, 2], yvalues2, circles = sppsub2[, 
                                                              7], inches = 0.1, add = TRUE, xlim = c(0, 
                                                                                                     5), fg = col2, bg = grpcol2, lwd = 2)
          }
        }
      }
    }
  }
  else {
    if (imax) {
      if (z1) {
        symbols(sppsub1[, 1], yvalues1, circles = sppsub1[, 
                                                          7], inches = 0.1, add = TRUE, xlim = c(0, 5), 
                fg = col1, bg = grpcol, lwd = 2)
      }
      if (z2) {
        symbols(sppsub2[, 1], yvalues2, circles = sppsub2[, 
                                                          7], inches = 0.1, add = TRUE, xlim = c(0, 5), 
                fg = col2, bg = grpcol2, lwd = 2)
      }
    }
    else {
      if (z1) {
        symbols(sppsub1[, 2], yvalues1, circles = sppsub1[, 
                                                          7], inches = 0.1, add = TRUE, xlim = c(0, 5), 
                fg = col1, bg = grpcol, lwd = 2)
      }
      if (z2) {
        symbols(sppsub2[, 2], yvalues2, circles = sppsub2[, 
                                                          7], inches = 0.1, add = TRUE, xlim = c(0, 5), 
                fg = col2, bg = grpcol2, lwd = 2)
      }
    }
  }
  if (z1) {
    axis(2, at = yvalues1, labels = rownames(sppsub1), las = 1, 
         mgp = c(1, 0.1, 0), cex.axis = cex.taxa, tck = tck)
  }
  if (z2) {
    axis(4, at = yvalues2, labels = rownames(sppsub2), mgp = c(1, 
                                                               0.1, 0), las = 1, cex.axis = cex.taxa, tck = tck)
  }
  if (log == "x") {
    axis(1, at = at, mgp = c(1, 0.1, 0), cex.axis = cex.axis, 
         tck = tck)
  }
  else {
    if (prtty) {
      axis(1, pretty(xmin:xmax, ntick), mgp = c(1, 0.1, 
                                                0), cex.axis = cex.axis, tck = tck)
    }
    else {
      axis(1, at = seq(from = round(xmin, digits = dig), 
                       to = round(xmax, digits = dig), by = round((xmax - 
                                                                     xmin)/4, digits = dig)), mgp = c(1, 0.1, 0), 
           cex.axis = cex.axis, tck = tck)
    }
  }
  mtext(xlabel, side = 1, line = 1.5, cex = cex)
  
  box(which = "plot")
  if (z1 & z2 & write) {
    sppsub <- list(sppsub1, sppsub2)
    names(sppsub) <- c("sppsub1", "sppsub2")
    return(sppsub)
  }
  if (z1 & write) {
    return(sppsub1)
  }
  if (z2 & write) {
    return(sppsub2)
  }
}
