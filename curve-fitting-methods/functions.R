laplacian.eigenmap <- function(W, d = 2, normalized = TRUE) {
  n <- nrow(W)
  if (normalized) {
    L <- normalized.laplacian(W)
  } else {
    L <- graph.laplacian(W)
  }
  L.eigen <- eigen(L, symmetric = TRUE)
  sweep(L.eigen$vectors[, seq(n - 1, n - d)], 2, 
        sqrt(L.eigen$values[seq(n - 1, n - d)]), `/`)
}

graph.laplacian <- function(W) {
  return(diag(colSums(W)) - W)
}

spectral.clustering <- function(W, K = 2, d = K, normalized = TRUE) {
  Y <- laplacian.eigenmap(W, d, normalized)
  kmeans(Y, K)$cluster
}

param.curve.quadr.2 <- function(t., 
                                beta10, beta11, beta12, 
                                beta20, beta21, beta22) {
  x1 <- beta10 + beta11 * t. + beta12 * t. ^ 2
  x2 <- beta20 + beta21 * t. + beta22 * t. ^ 2
  return(cbind(x1, x2))
}

estimate.one.t.quadr.2 <- function(x1, x2, 
                                   beta10, beta11, beta12, 
                                   beta20, beta21, beta22) {
  gamma0 <- beta10 * beta11 + beta20 * beta21 - beta11 * x1 - beta21 * x2
  gamma1 <- 2 * beta10 * beta12 + 2 * beta20 * beta22 + 
    beta11 ^ 2 + beta21 ^ 2 - 
    2 * beta12 * x1 - 2 * beta22 * x2
  gamma2 <- 3 * beta11 * beta12 + 3 * beta21 * beta22
  gamma3 <- 2 * (beta12 ^ 2 + beta22 ^ 2)
  roots <- polyroot(c(gamma0, gamma1, gamma2, gamma3))
  roots <- roots[abs(roots) == abs(Re(roots))]
  roots <- Re(roots)
  # if (length(roots) > 1) print(roots)
  return(max(roots))
}

estimate.t.quadr.2 <- function(X, 
                               beta10, beta11, beta12, 
                               beta20, beta21, beta22) {
  apply(X, 1, function(x) {
    estimate.one.t.quadr.2(x[1], x[2], 
                           beta10, beta11, beta12, 
                           beta20, beta21, beta22)
  })
}

estimate.quadr.coefs.2 <- function(X, t.hat, intercept = TRUE) {
  sum.X <- colSums(X)
  sum.Xt <- colSums(sweep(X, 1, t.hat, `*`))
  sum.Xt2 <- colSums(sweep(X, 1, t.hat ^ 2, `*`))
  sum.t <- sum(t.hat)
  sum.t2 <- sum(t.hat ^ 2)
  sum.t3 <- sum(t.hat ^ 3)
  sum.t4 <- sum(t.hat ^ 4)
  
  if (intercept) {
    A.sub <- matrix(c(1, sum.t, sum.t2, 
                      sum.t, sum.t2, sum.t3, 
                      sum.t2, sum.t3, sum.t4), 
                    nrow = 3, ncol = 3)
    A <- Matrix::bdiag(A.sub, A.sub)
    
    b <- c(sum.X[1], sum.Xt[1], sum.Xt2[1], 
           sum.X[2], sum.Xt[2], sum.Xt2[2])
    return(solve(A, b))
  } else {
    A.sub <- matrix(c(sum.t2, sum.t3, sum.t3, sum.t4),
                    nrow = 2, ncol = 2)
    A <- Matrix::bdiag(A.sub, A.sub)
    b <- c(sum.Xt[1], sum.Xt2[1], sum.Xt[2], sum.Xt2[2])
    theta <- solve(A, b)
    return(c(0, theta[1:2], 0, theta[3:4]))
  }
}

estimate.quadr.curve.2 <- function(X, init.params, 
                                   intercept = TRUE,
                                   eps = 1e-6, maxit = 100) {
  if (missing(init.params)) {
    theta <- rnorm(6)
  } else {
    theta <- init.params
  }
  d.theta <- 1 / eps
  
  niter <- 0
  while (d.theta > eps) {
    theta.prev <- theta
    t.hat <- estimate.t.quadr.2(X,
                                theta[1], theta[2], theta[3],
                                theta[4], theta[5], theta[6])
    X.hat <- param.curve.quadr.2(t.hat, theta[1], theta[2], theta[3], 
                                 theta[4], theta[5], theta[6])
    theta <- estimate.quadr.coefs.2(X, t.hat, intercept)
    
    d.theta <- sum((theta - theta.prev) ^ 2) / sum(theta ^ 2)
    
    niter <- niter + 1
    if (niter >= maxit) {
      warning('failed to converge to a quadratic curve')
      break
    }
  }
  return(list(theta = theta, 
              X = X.hat,
              t = t.hat))
}

construct.bezier.model.matrix <- function(t.hat, r = 2) {
  sapply(seq(0, r), function(s) {
    choose(r, s) * (1 - t.hat) ^ (r - s) * t.hat ^ s
  })
}

estimate.bezier.curve.2 <- function(X, 
                                    init.params, 
                                    weights,
                                    eps = 1e-6, maxit = 100) {
  n <- nrow(X)
  r <- 2
  
  if (missing(init.params)) {
    theta <- rnorm(6)
  } else {
    theta <- init.params
  }
  if (missing(weights)) {
    W <- diag(n)
  } else {
    W <- diag(weights)
  }
  
  d.mse <- eps + 1
  
  
  niter <- 0
  while (d.mse > eps) {
    mse.prev <- mse
    t.hat <- estimate.t.quadr.2(X,
                                theta[1], theta[2], theta[3],
                                theta[4], theta[5], theta[6])
    # fit <- bezier::bezierCurveFit(X[order(t.hat, decreasing = FALSE), ],
    #                               max.control.points = 3)
    # theta <- extract.quad.params.from.bezier.fit(fit$p)
    T <- construct.bezier.model.matrix(t.hat, r)
    p <- solve(t(T) %*% W %*% T, t(T) %*% W %*% X)
    theta <- extract.quad.params.from.bezier.fit(lapply(seq(ncol(p)),
                                                        function(i) p[, i]))
    
    X.hat <- param.curve.quadr.2(t.hat, 
                                 theta[1], theta[2], theta[3], 
                                 theta[4], theta[5], theta[6])
    # plot(X, asp = 1)
    # lines(X.hat[order(t.hat), ])
    mse <- sum((X - X.hat) ^ 2)
    # print(mse)
    d.mse <- abs((mse - mse.prev) / mse)
    # d.mse <- (mse.prev - mse) / mse.prev
    # print(d.mse)
    # d.theta <- sum((theta - theta.prev) ^ 2) / (r + 1)
    
    niter <- niter + 1
    if (niter >= maxit) {
      warning('failed to converge to a bezier curve')
      break
    }
  }
  
  # t.hat[t.hat < 0] <- 0
  # t.hat[t.hat > 1] <- 1
  # X.hat <- bezier::bezier(t.hat, fit$p)
  # X.hat <- param.curve.quadr.2(t.hat, 
  #                              theta[1], theta[2], theta[3], 
  #                              theta[4], theta[5], theta[6])
  
  return(list(theta = theta, 
              X = X.hat,
              t = t.hat))
}

compute.distances.quadr.2 <- function(X, theta) {
  t. <- estimate.t.quadr.2(X, 
                           theta[1], theta[2], theta[3], 
                           theta[4], theta[5], theta[6])
  X.hat <- param.curve.quadr.2(t., 
                               theta[1], theta[2], theta[3], 
                               theta[4], theta[5], theta[6])
  apply(X - X.hat, 1, function(x) sum(x ^ 2))
}

extract.quad.params.from.bezier.fit <- function(p) {
  p1 <- p[[1]]
  p2 <- p[[2]]
  c(p1[1], 2 * (p1[2] - p1[1]), p1[1] - 2 * p1[2] + p1[3],
    p2[1], 2 * (p2[2] - p2[1]), p2[1] - 2 * p2[2] + p2[3])
}

manifold.clustering.quadr.2 <- function(A, K = 2, d = 2, 
                                        method = 'bezier',
                                        initialization = 'random',
                                        intercept = TRUE, 
                                        maxit = 200) {
  n <- nrow(A)
  
  if (length(initialization) == n) {
    z.hat <- initialization
  } else if (initialization == 'spectral') {
    z.hat <- spectral.clustering(A, K)
  } else if (initialization == 'random') {
    z.hat <- sample(seq(K), n, replace = TRUE)
  } else {
    stop('choose a valid initialization')
  }
  X <- embedding(A, d, 0)
  # if (mean(X[, 1]) < 0) X[, 1] <- -X[, 1]
  z.hat.prev <- sample(seq(K), n, replace = TRUE)
  
  niter <- 0
  
  while (!all(z.hat.prev == z.hat)) {
    z.hat.prev <- z.hat
    
    X1 <- X[z.hat == 1, ]
    X2 <- X[z.hat == 2, ]
    
    if (method == 'polynomial') {
      curve1 <- estimate.quadr.curve.2(X1, intercept = intercept)
      curve2 <- estimate.quadr.curve.2(X2, intercept = intercept)
      theta1 <- curve1$theta
      theta2 <- curve2$theta
      d1 <- compute.distances.quadr.2(X, theta1)
      d2 <- compute.distances.quadr.2(X, theta2)
      distances <- cbind(d1, d2)
    } else if (method == 'bezier') {
      if ('curve1' %in% ls()) {
        curve1 <- estimate.bezier.curve.2(X1)
        curve2 <- estimate.bezier.curve.2(X2)
      } else {
        curve1 <- estimate.bezier.curve.2(X1, curve1$theta)
        curve2 <- estimate.bezier.curve.2(X2, curve2$theta)
      }
      theta1 <- curve1$theta
      theta2 <- curve2$theta
      d1 <- compute.distances.quadr.2(X, theta1)
      d2 <- compute.distances.quadr.2(X, theta2)
      distances <- cbind(d1, d2)
    } else {
      stop('choose a valid method')
    }
    
    z.hat <- apply(distances, 1, which.min)
    
    niter <- niter + 1
    if (niter >= maxit) {
      warning('failed to converge to a clustering')
      break
    }
  }
  
  return(list(z = z.hat, 
              X = list(curve1$X, curve2$X),
              t = list(curve1$t, curve2$t),
              theta = list(theta1, theta2),
              niter = niter))
}
