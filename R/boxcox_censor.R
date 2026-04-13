######################################################################
########### Box-Cox with Right/Double censoring ############
######################################################################

boxcox_censor <- function(x, y, statusx, statusy) {
  
  xcen <- x
  ycen <- y
  eps <- 1e-12
  
  bc_transform <- function(data, lambda) {   # Box-Cox transformation function
    if (abs(lambda) < 1e-6) log(data)
    else (data^lambda - 1)/lambda
  }
  
  fx <- function(x, g) {
    z <- bc_transform(x, g[5])
    dnorm(z, g[1], g[2]) * (x^(g[5] - 1))
  }
  
  Sx <- function(x, g) {
    z <- bc_transform(x, g[5])
    pnorm(z, g[1], g[2], lower.tail = FALSE)
  }
  
  Fx <- function(x, g) {
    z <- bc_transform(x, g[5])
    pnorm(z, g[1], g[2])
  }
  
  fy <- function(y, g) {
    z <- bc_transform(y, g[5])
    dnorm(z, g[3], g[4]) * (y^(g[5] - 1))
  }
  
  Sy <- function(y, g) {
    z <- bc_transform(y, g[5])
    pnorm(z, g[3], g[4], lower.tail = FALSE)
  }
  
  Fy <- function(y, g) {
    z <- bc_transform(y, g[5])
    pnorm(z, g[3], g[4])
  }
  
  logL <- function(g) {
    val_x <- sum((statusx == 0)  * log(fx(xcen, g) + eps) +
                   (statusx == 1)  * log(Sx(xcen, g) + eps) +
                   (statusx == -1) * log(Fx(xcen, g) + eps))
    val_y <- sum((statusy == 0)  * log(fy(ycen, g) + eps) +
                   (statusy == 1)  * log(Sy(ycen, g) + eps) +
                   (statusy == -1) * log(Fy(ycen, g) + eps))
    val <- -(val_x + val_y)   # Minimizing  Negative Log-likelihood
    if (!is.finite(val)) return(1e10)   # If NegLogL value is Inf, then return 1e10 instead of crashing
    # And move-on from the given set of parameter values
    return(val)
  }
  
  initial <- c(mean(xcen), sd(xcen), mean(ycen), sd(ycen), 1)
  lower <- c(-Inf, 1e-4, -Inf, 1e-4, -Inf)
  upper <- c( Inf,  Inf,  Inf,  Inf,  Inf)
  
  fit <- optim(par = initial, fn = logL, method = "L-BFGS-B",
               lower = lower, upper = upper,
               control = list(maxit = 10000, pgtol = 1e-6))
  
  return(list(g = fit$par, gval = fit$value, exitflag = fit$convergence))
}
