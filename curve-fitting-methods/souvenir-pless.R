import::from(magrittr, `%>%`)
library(ggplot2)

source('https://mtrosset.pages.iu.edu/Courses/675/graph.r')
source('https://mtrosset.pages.iu.edu/Courses/675/stress.r')
source('~/dev/pabm-grdpg/functions.R')
source('~/dev/dissertation/curve-fitting-methods/functions.R')

weighted.mds <- function(D, d = 2, w = rep(1, ncol(D))) {
  n <- ncol(D)
  e <- rep(1, n)
  H <- diag(n) - tcrossprod(e, w)
  B <- -H %*% D %*% H / 2
  eigen.B <- eigen(B, symmetric = TRUE)
  Y <- sweep(eigen.B$vectors[, seq(d)], 2, sqrt(eigen.B$values[seq(d)]), `*`)
  Y <- Y / (sum(abs(eigen.B$values)) / sum(eigen.B$values[seq(d)]))
  return(Y)
}

compute.dist2.manifold <- function(Dg, De, w) {
  w.sum <- sum(w)
  (rowSums(sweep(Dg - De, 2, w, `*`)) / w.sum) ^ 2
}

update.weights <- function(D2, sigma2 = 1e-1) {
  W <- exp(-D2 / sigma2) %>% 
    sweep(., 1, rowSums(.), `/`) %>% 
    return()
}

f1 <- function(t) {
  x1 <- t ^ 2
  x2 <- 2 * t - 2 * t ^ 2
  return(cbind(x1, x2))
}

f2 <- function(t) {
  x1 <- 2 * t - 2 * t ^ 2
  x2 <- 1 - 2 * t + t ^ 2
  return(cbind(x1, x2))
}

n1 <- 2 ** 7
n2 <- n1
n <- n1 + n2
z <- c(rep(1, n1), rep(2, n2))

a <- 1
b <- 1
t1 <- rbeta(n1, a, b)
t2 <- rbeta(n2, a, b)
T <- rbind(t1, t2)

X1 <- f1(t1)
X2 <- f2(t2)
X <- rbind(X1, X2)
plot(X, col = z, xlim = c(0, 1), ylim = c(0, 1), asp = 1)

P <- X %*% t(X)
diag(P) <- 0

A <- draw.graph(P)
Xhat <- embedding(P, 2, 0)
plot(Xhat, col = z, asp = 1)

D <- Xhat %>% 
  dist() %>% 
  as.matrix() %>% 
  graph.eps(.2) %>% 
  graph.short()

W <- runif(n * 2) %>% 
  matrix(nrow = n, ncol = 2) %>% 
  sweep(., 1, rowSums(.), `/`)

Y1 <- weighted.mds(D, 2, W[, 1])
Y2 <- weighted.mds(D, 2, W[, 2])

d1 <- compute.dist2.manifold(D, as.matrix(dist(Y1)), W[, 1])
d2 <- compute.dist2.manifold(D, as.matrix(dist(Y2)), W[, 2])

W <- update.weights(cbind(d1, d2), 1e1)

zhat <- apply(W, 1, which.max)
if (mean(z == zhat) < .5) zhat <- 3 - zhat
table(z, zhat)
plot(Xhat, col = (z == zhat) + 1)

plot(Xhat, col = zhat, asp = 1)

eps <- 1e-3
W.prev <- W + 1 / eps
Dw <- 1 / eps
while (Dw > eps) {
  W.prev <- W
  Y1 <- weighted.mds(D, 2, W[, 1])
  Y2 <- weighted.mds(D, 2, W[, 2])
  d1 <- compute.dist2.manifold(D, as.matrix(dist(Y1)), W[, 1])
  d2 <- compute.dist2.manifold(D, as.matrix(dist(Y2)), W[, 2])
  W <- update.weights(cbind(d1, d2), 1e0)
  Dw <- norm(W - W.prev, 'F')
  print(Dw)
  print(table(z, apply(W, 1, which.max)))
}
