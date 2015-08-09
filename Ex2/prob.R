library(spam);
library(MASS);

# Calculate the log of the eta probability (proportional)
etaProb <- function(eta, y, kv, ku, u, e) {
  eta = c(eta);
  y = c(y);
  n = length(eta);
  
  return(-0.5 * kv * t(eta) %*% eta + kv * t(eta) %*% u + t(eta) %*% y - t(exp(eta)) %*% e);
}

# Calculate the log of the x probability (proportional)
xProb <- function(x, y, kv, ku, e) {
  x = c(x);
  y = c(y);
  n = length(x) / 2;
  
  Q <- blockMatrix(kv, ku, R);
  
  return(-0.5 * t(x) %*% Q %*% x + t(x) %*% c(rep(0,n), y) - exp(x) %*% c(rep(0,n), e));
}

# x log proposal density (proportional)
xQ <- function(x, x0, y, kv, ku, e) {
  x = c(x);
  y = c(y);
  n = length(x) / 2;
  
  Q <- blockMatrix(kv, ku, R);
  
  eta0 = x0[(length(x0)/2+1):length(x0)];
  cv <- e * exp(eta0);
  b  <- y + e * exp(eta0) * (eta0 - 1);
  
  L <- c(rep(0, n), b);
  
  diagC <- diag.spam(c(rep(0, n), cv));
  Q <- Q + diagC;
  
  return(dmvnorm.canonical(x, L, Q, log=TRUE));
}

# eta log proposal density (proportional)
etaQ <- function(eta, eta0, y, kv, ku, u, e) {
  eta0 = c(eta0);
  eta = c(eta);
  u = c(u);
  y = c(y)
  
  cv <- e * exp(eta0);
  b <- y + e * exp(eta0) * (eta0 - 1);
  L <- b + kv * u;
  Q <- diag.spam(kv + c(cv));
  
  return(dmvnorm.canonical(eta, L, Q, log=TRUE));
}