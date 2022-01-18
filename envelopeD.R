envelopeD <- function(cases, contr, nsims) {
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
    
    Kcases <- Kest(simcases, correction = "isotropic")$iso 
    Kcontrols <- Kest(simcontr, correction = "isotropic")$iso
    dsim <- Kcases - Kcontrols 
    Dsim <- cbind(Dsim, dsim)
  }
  qts <- apply(Dsim, 1, quantile, probs = c(0.025, 0.975))
  
  Kcases <- Kest(cases, correction = "isotropic")
  Kcontrols <- Kest(contr, correction = "isotropic")
  D <- Kcases$iso - Kcontrols$iso
  r <- Kcases$r
  plot(NULL, xlab = "r", ylab = "Estimated D(r)",
       xlim = range(r), 
       ylim = c(min(c(D, qts[1, ])), max(D, qts[1, ])))
  polygon(c(r, rev(r)), c(qts[1, ], rev(qts[2, ])), 
          col = "lightgrey", border = NA)
  lines(r, D)
  abline(h = 0, col = "red", lty = "dashed")
}