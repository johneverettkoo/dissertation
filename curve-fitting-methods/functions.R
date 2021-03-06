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
  # roots <- polyroot(c(gamma0, gamma1, gamma2, gamma3))
  # roots <- roots[abs(roots) == abs(Re(roots))]
  # roots <- Re(roots)
  roots <- optimize(
    function(t) crossprod(c(x1 - beta10 - beta11 * t - beta12 * t^2, 
                            x2 - beta20 - beta21 * t - beta22 * t^2)), 
    c(-1, 1)
  )$minimum
  if (length(roots) > 1) print(roots)
  if (abs(roots) == 1) print(roots)
  return(max(roots))
}

estimate.one.t.bezier <- function(x, p, min.t = 0, max.t = 2, intercept = TRUE, tol = 1e-9) {
  polynomial <- function(t, p) {
    if (intercept) {
      r <- nrow(p) - 1
      T <- sapply(seq(0, r), function(s) {
        choose(r, s) * (1 - t) ^ (r - s) * t ^ s
      })
    } else {
      r <- nrow(p)
      T <- sapply(seq(r), function(s) {
        choose(r, s) * (1 - t) ^ (r - s) * t ^ s
      })
    }
    as.numeric(T %*% p)
  }
  l <- function(t) log(as.numeric(crossprod(polynomial(t, p) - x)))
  minima <- optimize(l, c(min.t, max.t), tol = tol)$minimum
  if (abs(minima) == max.t) warning(minima)
  if (length(minima) > 1) warning(minima)
  return(max(minima))
}

estimate.t.bezier <- function(X, p, t.prev, intercept = TRUE, tol = 1e-9) {
  if (missing(t.prev)) {
    t.hat <- apply(X, 1, function(x) {
      estimate.one.t.bezier(x, p, tol = tol, intercept = intercept)
    })
  } else {
    t.hat <- sapply(seq_along(t.prev), function(i) {
      estimate.one.t.bezier(X[i, ], p, t.prev[i] - 1, t.prev[i] + 1, 
                            tol = tol, intercept = intercept)
    })
    loss.new <- sapply(seq_along(t.hat),
                       function(i) loss.one.t(X[i, ], t.hat[i], p))
    loss.prev <- sapply(seq_along(t.prev),
                        function(i) loss.one.t(X[i, ], t.prev[i], p))
    prop.change <- (mean(loss.new < loss.prev))
    if (prop.change < 1) {
      warning(paste('Failed to find new t for', 1 - prop.change))
    }
    t.hat <- ifelse(loss.new <= loss.prev, t.hat, t.prev)
  }
  
  # if (mean(t.hat) < 0) t.hat <- -t.hat
  # t.hat[t.hat < 0] <- 0
  
  return(t.hat)
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
  
  n <- nrow(X)
  w <- n ^ 2 - 1
  
  niter <- 0
  while (d.theta > eps) {
    theta.prev <- theta
    t.hat <- estimate.t.quadr.2(X,
                                theta[1], theta[2], theta[3],
                                theta[4], theta[5], theta[6])
    # normalize to [0, 1] using umvue
    min.t <- min(t.hat)
    max.t <- max(t.hat)
    a <- (n * (n + 1) * min.t - (n + 1) * max.t) / w
    b <- (-(n + 1) * min.t + n * (n + 1) * max.t) / w
    t.hat <- (t.hat - a) / (b - a)
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

construct.bezier.model.matrix <- function(t.hat, r = 2, intercept = TRUE) {
  if (intercept) {
    sapply(seq(0, r), function(s) {
      choose(r, s) * (1 - t.hat) ^ (r - s) * t.hat ^ s
    })
  } else {
    sapply(seq(r), function(s) {
      choose(r, s) * (1 - t.hat) ^ (r - s) * t.hat ^ s
    })
  }
}

bezier.curve <- function(t, p, intercept = TRUE) {
  if (intercept) {
    r <- nrow(p) - 1
  } else {
    r <- nrow(p)
  }
  T <- construct.bezier.model.matrix(t, r, intercept = intercept)
  return(T %*% p)
}

normalize.umvue <- function(t.hat) {
  n <- length(t.hat)
  w <- n ^ 2 - 1
  min.t <- min(t.hat)
  max.t <- max(t.hat)
  a <- (n * (n + 1) * min.t - (n + 1) * max.t) / w
  b <- (-(n + 1) * min.t + n * (n + 1) * max.t) / w
  return((t.hat - a) / (b - a))
}

normalize.ecdf <- function(t.hat) {
  ecdf(t.hat)(t.hat)
}

estimate.bezier.curve.2 <- function(X, 
                                    init.params, 
                                    weights,
                                    eps = 1e-3, 
                                    maxit = 200,
                                    initialization = 'isomap', 
                                    normalize = TRUE, 
                                    normalize.method = 'umvue',
                                    eps.isomap = 2e-1, 
                                    intercept = TRUE) {
  n <- nrow(X)
  d <- ncol(X)
  r <- 2
  
  if (missing(weights)) {
    weights <- rep(1, n)
  }
  
  if (missing(init.params)) {
    if (initialization == 'random') {
      if (intercept) {
        p <- matrix(rnorm(d * (r + 1)), nrow = r + 1, ncol = d)
      } else {
        p <- matrix(rnorm(d * r), nrow = r, ncol = d)
      }
      t.hat <- estimate.t.bezier(X, p, intercept = intercept)
    } else if (initialization == 'isomap') {
      isomap.out <- estimate.bezier.curve.isomap(X, eps.isomap, weights, 
                                                 intercept = intercept)
      p <- isomap.out$p
      t.hat <- isomap.out$t
    } else {
      stop('initialization must be random or isomap')
    }
  } else {
    p <- init.params
    t.hat <- estimate.t.bezier(X, p, intercept = intercept)
  }
  
  if (normalize) {
    if (normalize.method == 'umvue') {
      t.hat <- normalize.umvue(t.hat)
    } else if (normalize.method == 'ecdf') {
      t.hat <- normalize.ecdf(t.hat)
    } else {
      stop('normalize.method must be umvue or ecdf')
    }
  }
  
  W <- diag(weights)
  
  mse <- bezier.mse(X, t.hat, p)
  d.mse <- 1 / eps
  
  mse.b <- c()
  update.what <- c()
  
  niter <- 0
  while (d.mse > eps) {
    mse.prev <- mse
    
    T <- construct.bezier.model.matrix(t.hat, r, intercept = intercept)
    p <- solve(t(T) %*% W %*% T, t(T) %*% W %*% X)
    # X.hat <- bezier.curve(t.hat, p)
    # print(bezier.mse(X, t.hat, p))
    mse.b <- c(mse.b, bezier.mse(X, t.hat, p))
    update.what <- c(update.what, 'p')
    # plot(X, asp = 1, main = 'update p')
    # lines(X.hat[order(t.hat), ], col = 2)
    
    t.hat <- estimate.t.bezier(X, p, t.hat, intercept = intercept)
    summary(t.hat)
    # print(bezier.mse(X, t.hat, p))
    mse.b <- c(mse.b, bezier.mse(X, t.hat, p))
    update.what <- c(update.what, 't')
    # X.hat <- bezier.curve(t.hat, p)
    # plot(X, asp = 1, main = 'update t')
    # lines(X.hat[order(t.hat), ], col = 2)
    # 
    mse <- bezier.mse(X, t.hat, p)
    
    if (normalize) {
      if (normalize.method == 'umvue') {
        t.hat <- normalize.umvue(t.hat)
      } else if (normalize.method == 'ecdf') {
        t.hat <- normalize.ecdf(t.hat)
      } else {
        stop('normalize.method must be umvue or ecdf')
      }
    }
    # print(mse)
    # d.mse <- abs((mse - mse.prev) / mse)
    if (mse < eps ^ 3) break
    d.mse <- (mse.prev - mse) / mse.prev
    # print(d.mse)
    # print(p)
    
    niter <- niter + 1
    if (niter >= maxit) {
      warning('failed to converge to a bezier curve')
      break
    }
  }
  
  X.hat <- bezier.curve(t.hat, p)
  theta <- extract.quad.params.from.bezier.fit(lapply(seq(ncol(p)),
                                                      function(i) p[, i]),
                                               intercept = intercept)
  
  return(list(theta = theta, 
              p = p, 
              X = X.hat,
              t = t.hat))
}

knn.graph <- function(X, k) {
  dist.matrix <- mds.edm1(X)
  knn.graph <- graph.knn(dist.matrix, k) * dist.matrix
  dist.matrix[dist.matrix == 0] <- Inf
  diag(dist.matrix) <- 0
  return(graph.short(dist.matrix))
}

estimate.bezier.curve.isomap <- function(X, k = lower(sqrt(nrow(X))), eps = .1, 
                                         weights, intercept = TRUE) {
  n <- nrow(X)
  d <- 2
  r <- 2
  
  w <- n ^ 2 - 1
  
  if (missing(weights)) {
    W <- diag(n)
  } else {
    W <- diag(weights)
  }
  
  # D <- X %>% 
  #   dist() %>% 
  #   as.matrix() %>% 
  #   graph.eps(eps) %>% 
  #   graph.short()
  D <- knn.graph(X, k)
  
  t.hat <- as.vector(cmdscale(sqrt(D), 1))
  t.hat <- ecdf(t.hat)(t.hat)
  # min.t <- min(t.hat)
  # max.t <- max(t.hat)
  # a <- (n * (n + 1) * min.t - (n + 1) * max.t) / w
  # b <- (-(n + 1) * min.t + n * (n + 1) * max.t) / w
  # t.hat <- (t.hat - a) / (b - a)
  
  T <- construct.bezier.model.matrix(t.hat, r, intercept = intercept)
  p <- solve(t(T) %*% W %*% T, t(T) %*% W %*% X)
  
  X.hat <- bezier.curve(t.hat, p, intercept = intercept)
  # plot(X, asp = 1)
  # lines(X.hat[order(t.hat), ])
  
  theta <- extract.quad.params.from.bezier.fit(lapply(seq(ncol(p)),
                                                      function(i) p[, i]),
                                               intercept = intercept)
  
  return(list(theta = theta, 
              p = p, 
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

extract.quad.params.from.bezier.fit <- function(p, intercept = TRUE) {
  p1 <- p[[1]]
  p2 <- p[[2]]
  if (!intercept) {
    p1 <- c(0, p1)
    p2 <- c(0, p2)
  }
  c(p1[1], 2 * (p1[2] - p1[1]), p1[1] - 2 * p1[2] + p1[3],
    p2[1], 2 * (p2[2] - p2[1]), p2[1] - 2 * p2[2] + p2[3])
}

manifold.clustering.quadr.2 <- function(A, K = 2, p = 2, q = 0, 
                                        initialization = 'random',
                                        intercept = TRUE, 
                                        maxit = 200,
                                        eps.isomap = 1e-2, 
                                        eps = 1e-3) {
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
  X <- embedding(A, p, q)
  # if (mean(X[, 1]) < 0) X[, 1] <- -X[, 1]
  z.hat.prev <- sample(seq(K), n, replace = TRUE)
  
  niter <- 0
  
  loss <- c()
  d.loss <- eps + 1
  new.loss <- eps + 1
  
  while (d.loss > eps) {
    z.hat.prev <- z.hat
    
    X1 <- X[z.hat == 1, ]
    X2 <- X[z.hat == 2, ]
    
    curve1 <- estimate.bezier.curve.2(X1, eps.isomap = eps.isomap, intercept = intercept)
    curve2 <- estimate.bezier.curve.2(X2, eps.isomap = eps.isomap, intercept = intercept)
    # if (!exists('curve1')) {
    #   curve1 <- estimate.bezier.curve.2(X1)
    #   curve2 <- estimate.bezier.curve.2(X2)
    # } else {
    #   curve1 <- estimate.bezier.curve.2(X1, curve1$p)
    #   curve2 <- estimate.bezier.curve.2(X2, curve2$p)
    # }
    # curve1 <- estimate.bezier.curve.isomap(X1, .25)
    # curve2 <- estimate.bezier.curve.isomap(X2, .25)
    theta1 <- curve1$theta
    theta2 <- curve2$theta
    d1 <- compute.distances.quadr.2(X, theta1)
    d2 <- compute.distances.quadr.2(X, theta2)
    distances <- cbind(d1, d2)
    
    prev.loss <- new.loss
    new.loss <- clustering.loss(X, z.hat, 
                                list(curve1$p, curve2$p),
                                list(curve1$t, curve2$t))
    d.loss <- abs(prev.loss - new.loss) / new.loss
    loss <- c(loss, new.loss)
    
    z.hat <- apply(distances, 1, which.min)
    # table(z.hat, z)
    # plot(X, col = z.hat)
    # loss
    
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
              niter = niter,
              loss = loss))
}

bezier.mse <- function(X, t, p, intercept = TRUE) {
  if (intercept) {
    r <- nrow(p) - 1L
  } else {
    r <- nrow(p)
  }
  T <- construct.bezier.model.matrix(t, r, intercept = intercept)
  # norm(X - T %*% p, 'F') ^ 2 / prod(dim(X))
  norm(X - T %*% p, 'F') ^ 2
}

loss.one.t <- function(x, t, p, intercept = TRUE) {
  if (intercept) {
    r <- nrow(p) - 1
    T <- sapply(seq(0, r), function(s) {
      choose(r, s) * (1 - t) ^ (r - s) * t ^ s
    })
  } else {
    r <- nrow(p)
    T <- sapply(seq(r), function(s) {
      choose(r, s) * (1 - t) ^ (r - s) * t ^ s
    })
  }
  g <- as.numeric(T %*% p)
  as.numeric(crossprod(x - g))
}

clustering.loss <- function(X, z.hat, p.list, t.list, intercept = TRUE) {
  K <- max(z.hat)
  sapply(seq(K), function(k) {
    X.k <- X[z.hat == k, ]
    p.k <- p.list[[k]]
    t.k <- t.list[[k]]
    bezier.mse(X.k, t.k, p.k, intercept = intercept)
  }) %>% 
    sum()
}
