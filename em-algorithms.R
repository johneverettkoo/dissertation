em.sbm <- function(A, pi.start, theta.start, maxit = 100, eps = 1e-6) {
  pi <- pi.start
  theta <- theta.start
  
  n <- ncol(A)
  K <- ncol(theta)
  
  iter <- 0
  
  while(TRUE) {
    pi.old <- pi
    theta.old <- theta
    # e-step
    pi <- update.pi.sbm(pi.old, A, theta)
    
    # m-step
    theta <- update.theta.sbm(pi, A)
    
    if (Matrix::norm(pi.old - pi, 'F') < eps) break
    
    iter <- iter + 1
    if (iter >= maxit) {
      warning('failed to converge')
      break
    }
  }
  
  return(list(pi = pi,
              theta = theta, 
              iter = iter))
}

update.pi.sbm.ik <- function(pi, A, theta, i, k) {
  log.pi.ik <- 0
  n <- ncol(A)
  K <- ncol(theta)
  for (j in seq(n)) {
    for (l in seq(K)) {
      if (j != i) {
        log.pi.ik <- log.pi.ik + 
          pi[j, l] * 
          (A[i, j] * log(theta[k, l]) + 
             (1 - A[i, j]) * log(1 - theta[k, l]))
      }
    }
  }
  return(log.pi.ik)
}

update.pi.sbm <- function(pi, A, theta) {
  pi.new <- pi
  n <- nrow(pi)
  K <- ncol(pi)
  for (i in seq(n)) {
    for (k in seq(K)) {
      pi.new[i, k] <- update.pi.sbm.ik(pi, A, theta, i, k)
    }
  }
  pi.new <- pi.new - mean(pi.new)
  pi.new <- exp(pi.new)
  pi.new <- t(apply(pi.new, 1, function(x) x / sum(x)))
  return(pi.new)
}

update.theta.sbm.kl <- function(pi, A, k, l) {
  n <- nrow(pi)
  numerator <- 0
  denominator <- 0
  for (i in seq(2, n)) {
    for (j in seq(i - 1)) {
      numerator <- numerator + A[i, j] * pi[i, k] * pi[j, l]
      denominator <- denominator + pi[i, k] * pi[j, l]
    }
  }
  return(numerator / denominator)
}

update.theta.sbm <- function(pi, A) {
  K <- ncol(pi)
  theta.new <- matrix(0, nrow = K, ncol = K)
  for (k in seq(K)) {
    for (l in seq(K)) {
      theta.new[k, l] <- update.theta.sbm.kl(pi, A, k, l)
    }
  }
  return(theta.new)
}

em.dcbm <- function(A, pi.start, theta.start, omega.start, 
                    maxit = 100, eps = 1e-6) {
  pi <- pi.start
  theta <- theta.start
  omega <- omega.start
  
  n <- nrow(pi)
  K <- ncol(pi)
  
  iter <- 0
  
  while(TRUE) {
    pi.old <- pi
    theta.old <- theta
    omega.old <- omega
    
    if (max(theta) > 1) {
      theta <- theta / max(theta)
    }
    if (max(omega) > 1) {
      omega <- omega / max(omega)
    }
    pi <- update.pi.dcbm(pi.old, A, theta, omega)
    theta <- update.theta.dcbm(pi, A, omega)
    omega <- update.omega.dcbm(pi, A, theta, omega)
    
    if (Matrix::norm(pi.old - pi, 'F') < eps) break
    
    iter <- iter + 1
    if (iter >= maxit) {
      warning('failed to converge')
      break
    }
  }
  return(list(pi = pi,
              theta = theta,
              omega = omega,
              iter = iter))
}

update.pi.dcbm <- function(pi, A, theta, omega) {
  pi.new <- pi
  n <- nrow(pi)
  K <- ncol(pi)
  for (i in seq(n)) {
    for (k in seq(K)) {
      pi.new[i, k] <- update.pi.dcbm.ik(pi, A, theta, omega, i, k)
    }
  }
  pi.new <- pi.new - mean(pi.new)
  pi.new <- exp(pi.new)
  pi.new <- t(apply(pi.new, 1, function(x) x / sum(x)))
  return(pi.new)
}

update.pi.dcbm.ik <- function(pi, A, theta, omega, i, k) {
  log.pi.ik <- 0
  n <- ncol(A)
  K <- ncol(theta)
  for (j in seq(n)) {
    for (l in seq(K)) {
      if (j != i) {
        log.pi.ik <- log.pi.ik + 
          pi[j, l] * 
          (A[i, j] * log(theta[k, l] * omega[i] * omega[j]) + 
             (1 - A[i, j]) * log(1 - theta[k, l] * omega[i] * omega[j]))
      }
    }
  }
  return(log.pi.ik)
}

update.theta.dcbm <- function(pi, A, omega) {
  K <- ncol(pi)
  theta.new <- matrix(0, nrow = K, ncol = K)
  for (k in seq(K)) {
    for (l in seq(K)) {
      theta.new[k, l] <- update.theta.dcbm.kl(pi, A, omega, k, l)
    }
  }
  return(theta.new)
}

update.theta.dcbm.kl <- function(pi, A, omega, k, l) {
  n <- nrow(pi)
  numerator <- 0
  denominator <- 0
  for (i in seq(2, n)) {
    for (j in seq(i - 1)) {
      numerator <- numerator + A[i, j] * pi[i, k] * pi[j, l]
      denominator <- denominator + pi[i, k] * pi[j, l] * omega[i] * omega[j]
    }
  }
  return(numerator / denominator)
}

update.omega.dcbm <- function(pi, A, theta, omega) {
  n <- nrow(pi)
  K <- ncol(pi)
  omega.new <- omega
  for (i in seq(n)) {
    omega.new[i] <- update.omega.dcbm.i(pi, A, theta, omega, i)
  }
  return(omega.new)
}

update.omega.dcbm.i <- function(pi, A, theta, omega, i) {
  n <- nrow(pi)
  K <- ncol(pi)
  numerator <- 0
  denominator <- 0
  for (j in seq(n)) {
    if (j != i) {
      for (k in seq(K)) {
        for (l in seq(K)) {
          numerator <- numerator + pi[i, k] * pi[j, l] * A[i, j]
          denominator <- denominator + pi[i, k] * pi[j, l] * theta[k, l] * omega[j]
        }
      }
    }
  }
  return(numerator / denominator)
}

em.pabm <- function(A, pi.start, lambda.start, 
                    maxit = 100, eps = 1e-9) {
  pi <- pi.start
  lambda <- lambda.start
  
  n <- nrow(pi)
  K <- ncol(pi)
  
  iter <- 0
  
  while(TRUE) {
    pi.old <- pi
    lambda.old <- lambda
    
    if (max(lambda) >= 1) {
      lambda <- lambda / max(lambda + 1e-6)
    }
    pi <- update.pi.pabm(pi.old, A, lambda)
    lambda <- update.lambda.pabm(pi, A, lambda)
    
    if (Matrix::norm(pi.old - pi, 'F') < eps) break
    
    iter <- iter + 1
    if (iter >= maxit) {
      warning('failed to converge')
      break
    }
  }
  return(list(pi = pi,
              lambda = lambda, 
              iter = iter))
}

update.pi.pabm <- function(pi, A, lambda) {
  n <- nrow(pi)
  K <- ncol(pi)
  pi.new <- pi
  
  for (i in seq(n)) {
    for (k in seq(K)) {
      pi.new[i, k] <- update.pi.pabm.ik(pi, A, lambda, i, k)
    }
  }
  
  pi.new <- pi.new - mean(pi.new)
  pi.new <- exp(pi.new)
  pi.new <- t(apply(pi.new, 1, function(x) x / sum(x)))
  return(pi.new)
}

update.pi.pabm.ik <- function(pi, A, lambda, i, k) {
  log.pi.ik <- 0
  n <- nrow(pi)
  K <- ncol(pi)
  for (j in seq(n)) {
    for (l in seq(K)) {
      if (j != i) {
        log.pi.ik <- log.pi.ik + 
          pi[j, l] * 
          (A[i, j] * log(lambda[i, k] * lambda[j, l]) + 
             (1 - A[i, j]) * log(1 - lambda[i, k] * lambda[j, l]))
      }
    }
  }
  return(log.pi.ik)
}

update.lambda.pabm <- function(pi, A, lambda) {
  n <- nrow(pi)
  K <- ncol(pi)
  lambda.new <- matrix(NA, nrow = n, ncol = K)
  
  for (i in seq(n)) {
    for (k in seq(K)) {
      lambda.new[i, k] <- update.lambda.pabm.ik(pi, A, lambda, i, k)
    }
  }
  
  return(lambda.new)
}

update.lambda.pabm.ik <- function(pi, A, lambda, i, k) {
  n <- nrow(pi)
  K <- ncol(pi)
  
  numerator <- 0
  denominator <- 0
  
  for (j in seq(n)) {
    if (j != i) {
      for (l in seq(K)) {
        if (l != k) {
          numerator <- numerator + 
            pi[i, k] * pi[j, l] * A[i, j]
          denominator <- denominator + 
            pi[i, k] * pi[j, l] * lambda[j, l]
        }
      }
    }
  }
  return(numerator / denominator)
}