envelopeD_jmn <- function(cases, contr, nsims, rsplit, Rmax) {
  Dsim <- c() # we'll store our result in this object
  win <- cases$window
  x <- c(cases$x, contr$x)
  y <- c(cases$y, contr$y)
  cc <- c(rep("case", cases$n), rep("control", contr$n))
  for (i in 1:nsims) {
    cc <- sample(cc)
    simcases <- ppp(x = x[cc == "case"], 
                    y = y[cc =="case"], 
                    window = win)
    simcontr <- ppp(x = x[cc == "control"], 
                    y = y[cc == "control"], 
                    window = win)
    # rescaling each to km
    # simcases.km <- rescale(simcases, 1000, "km")
    # simcontr.km <- rescale(simcases, 1000, "km")
    
    Kcases <- Kest(simcases, 
                   correction = "isotropic", 
                   r = rsplit,
                   rmax = Rmax)$iso 
    Kcontrols <- Kest(simcontr, 
                      correction = "isotropic", 
                      r = rsplit,
                      rmax = Rmax)$iso
    dsim <- Kcases - Kcontrols 
    Dsim <- cbind(Dsim, dsim)
  }
  qts <- apply(Dsim, 1, quantile, probs = c(0.025, 0.975))
  
  Kcases <- Kest(cases, 
                 correction = "isotropic", 
                 r = rsplit,
                 rmax = Rmax)
  Kcontrols <- Kest(contr, 
                    correction = "isotropic", 
                    r = rsplit,
                    rmax = Rmax)
  D <- Kcases$iso - Kcontrols$iso
  r <- Kcases$r
  plot(NULL, xlab = "Distance, Km", ylab = "Difference in K-functions",
       xlim = range(r), 
       ylim = c(-300, 250))
  polygon(c(r, rev(r)), c(qts[1, ], rev(qts[2, ])), 
          col = "lightgrey", border = NA)
  lines(r, D)
  abline(h = 0, col = "red", lty = "dashed")
}