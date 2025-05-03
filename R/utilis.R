#' Estimate mean outcome at balance.target via residual balancing
#'
#' @param XW the input features for the sub-population of interest
#' @param YW the observed responses for the sub-population of interest
#' @param balance.target the desired center of the dataset
#' @param allow.negative.weights whether negative gammas are allowed for balancing
#' @param zeta tuning parameter for selecting approximately balancing weights
#' @param fit.method the method used to fit mu(x) = E[YW | XW = x]
#' @param alpha tuning paramter for glmnet
#' @param optimizer which optimizer to use for approximate balancing
#' @param bound.gamma whether upper bound on gamma should be imposed
#' @param verbose whether the optimizer should print progress information
#' @param K hyperparameter for optimization
#'
#' @return Estimate for E[YW | XW = balance.target], along with variance estimate
#'
#' @export residualBalance.mean
residualBalance.mean = function(XW, YW,
                                balance.target,
                                allow.negative.weights = FALSE,
                                zeta,
                                K = 2.75,
                                fit.method = c("elnet", "none"),
                                alpha,
                                optimizer = c("mosek", "pogs", "pogs.dual", "quadprog", "pogs.new"),
                                bound.gamma = TRUE,
                                verbose = FALSE) {
  
  fit.method = match.arg(fit.method)
  optimizer = match.arg(optimizer)
  
  gamma = approx.balance(XW, balance.target, zeta = zeta, K, allow.negative.weights = allow.negative.weights, optimizer = optimizer, bound.gamma=bound.gamma, verbose=verbose)
  
  if (fit.method == "elnet") {
    
    lasso.fit = glmnet::cv.glmnet(XW, YW, alpha = alpha, nfold = 3)
    mu.lasso = predict(lasso.fit, newx = matrix(balance.target, 1, length(balance.target)))
    
    residuals = YW - predict(lasso.fit, newx = XW)
    mu.residual = sum(gamma * residuals)
    
    var.hat = sum(gamma^2 * residuals^2) *
      # degrees of freedom correction
      length(gamma) / max(1, length(gamma) - sum(coef(lasso.fit) != 0))
    
  } else if (fit.method == "none") {
    
    mu.lasso = 0
    mu.residual = sum(gamma * YW)
    
    var.hat = NA
    
  } else {
    
    stop("Invalid choice of fitting method.")
    
  }
  
  mu.hat = mu.lasso + mu.residual
  c(mu.hat, var.hat)
}

#' Compute approximately balancing weights
#'
#' Returns the minimizer of:
#'   (1 - zeta) ||gamma||^2 + zeta ||M'gamma - balance.target||_infty^2 (*)
#'
#' @param M the feature matrix, see (*)
#' @param balance.target the target solution, see (*)
#' @param zeta tuning parameter, see (*)
#' @param allow.negative.weights are the gammas allowed to be negative?
#' @param optimizer Which optimizer to use? Mosek is a commercial solver, but free
#'                  academic licenses are available. Needs to be installed separately.
#'                  Pogs runs ADMM and may be useful for large problems, and
#'                  must be installed separately. Quadprog is the default
#'                  R solver.
#' @param bound.gamma whether upper bound on gamma should be imposed
#' @param gamma.max specific upper bound for gamma (ignored if bound.gamma = FALSE)
#' @param verbose whether the optimizer should print progress information
#'
#' @return gamma, the minimizer of (*)
#'
#' @export approx.balance
approx.balance = function(M,
                          balance.target,
                          zeta = 0.5,
                          K = 1,
                          cv.K = FALSE,  # NEW: 是否对 K 进行交叉验证
                          allow.negative.weights = FALSE,
                          optimizer = c("mosek", "pogs", "pogs.dual", "quadprog", "pogs.new"),
                          bound.gamma = FALSE,
                          gamma.max = 1/nrow(M)^(2/3),
                          verbose = FALSE) {
  
  if (zeta <= 0 || zeta >= 1) {
    stop("approx.balance: zeta must be between 0 and 1")
  }
  
  optimizer = match.arg(optimizer)
  
  # 当 optimizer = "pogs.new" 且 cv.K = TRUE 时，进行搜索选择最优 K
  if (cv.K && optimizer == "pogs.new") {
    K_values <- seq(1, 3, length.out = 15)  # 生成 K 取值范围
    errors <- numeric(length(K_values))  # 存储 balance error
    
    for (i in seq_along(K_values)) {
      gamma_cv <- approx.balance.pogs.new(M, balance.target, K = K_values[i],
                                          allow.negative.weights = allow.negative.weights, 
                                          bound.gamma = bound.gamma, gamma.max = gamma.max)
      
      balance_error <- max(abs(t(M) %*% gamma_cv - balance.target))  # 计算 balance error
      errors[i] <- balance_error
    }
    
    # 选择 balance error 最小的 K
    K <- K_values[which.min(errors)]
    if (verbose) cat("Selected optimal K:", K, "\n")
  }
  
  if (optimizer == "mosek") {
    if (suppressWarnings(require("Rmosek", quietly = TRUE))) {
      gamma = approx.balance.mosek.dual(M, balance.target, zeta, allow.negative.weights, bound.gamma, gamma.max, verbose)
    } else {
      if (suppressWarnings(require("pogs", quietly = TRUE))) {
        warning("The mosek optimizer is not installed. Using pogs instead.")
        optimizer = "pogs"
      } else {
        warning("Neither mosek nor pogs optimizers are installed. Using quadprog instead.")
        optimizer = "quadprog"
      }
    }
  }
  
  if (optimizer %in% c("pogs", "pogs.dual", "pogs.new")) {
    if (suppressWarnings(require("pogs", quietly = TRUE))) {
      if (optimizer == "pogs") {
        gamma = approx.balance.pogs(M, balance.target, zeta, allow.negative.weights, bound.gamma, gamma.max, verbose)
      } else if (optimizer == "pogs.new"){
        gamma = approx.balance.pogs.new(M, balance.target, K, allow.negative.weights, bound.gamma, gamma.max)
      } else {
        if (bound.gamma) {warning("bound.gamma = TRUE not implemented for this optimizer")}
        gamma = approx.balance.pogs.dual(M, balance.target, zeta, allow.negative.weights, verbose)
      }
    } else {
      warning("The POGS optimizer is not installed. Using quadprog instead.")
      optimizer = "quadprog"
    }
  }
  
  if (optimizer == "quadprog") {
    if (bound.gamma) {warning("bound.gamma = TRUE not implemented for this optimizer")}
    gamma = approx.balance.quadprog(M, balance.target, zeta, allow.negative.weights)
  }
  
  gamma
}

#################################################################################################################
approx.balance.pogs.new = function(M,                # Design matrix X_c
                                   balance.target,   # Target vector 
                                   K,            # Trade-off parameter for balance constraints
                                   allow.negative.weights = FALSE,  # Whether gamma can be negative
                                   bound.gamma = FALSE,  # Whether to impose max gamma constraint
                                   verbose = FALSE) {
  
  # Ensure balance.target is a column vector
  balance.target = matrix(balance.target, ncol = 1) # xi
  
  # Compute the bound for the ||ξ - X_c^T γ||_∞ constraint
  nn = nrow(M)
  balance.bound = K * sqrt(log(ncol(M)) / (nn))
  
  # Compute the maximum allowed γ_i values, ensuring it's not too small
  gamma.max = max(nn^(-2/3))  # the gamma constraint
  
  # Define the quadratic loss function ||γ||_2^2
  g = list(h = kSquare())
  
  # Define constraints:
  # (1) Infinity norm constraint: ||ξ - X_c^T γ||_∞ ≤ balance.bound
  # (2) Sum constraint: sum(gamma_i) = 1
  f = list(h = c(kIndLe0(2 * ncol(M)), kIndEq0(2)),  # (1) 2p inequality constraints, (2) sum constraint
           b = c(rep(balance.bound, 2 * ncol(M)), 1, 1))  # Right-hand side values for the constraints
  
  # Construct the constraint matrix A
  # The first two blocks enforce ||M^T gamma - balance.target||_∞ ≤ balance.bound
  # The third block enforces sum(gamma) = 1
  A = rbind(
    cbind(t(M), -balance.target, -balance.bound),  # M^T γ - balance.target ≤ balance.bound
    cbind(-t(M), balance.target, -balance.bound),  # -M^T γ + balance.target ≤ balance.bound
    c(rep(1, nrow(M)), 0, 0),  # sum(gamma) = 1
    c(rep(0, nrow(M)), 1, 0)   # Auxiliary variable Z constraint
  )
  
  # If we do not allow negative weights, we enforce gamma ≥ 0
  if (!allow.negative.weights) {
    f$h = c(f$h, kIndLe0(nrow(M)))  # Add inequality constraint gamma ≥ 0
    f$b = c(f$b, rep(0, nrow(M)))  # Right-hand side for the new constraints
    A = rbind(A, cbind(diag(-1, nrow(M)), 0, 0))  # Add negative identity to enforce non-negativity
  }
  
  # If we bound gamma, we enforce gamma_i ≤ gamma.max
  if (bound.gamma) {
    f$h = c(f$h, kIndLe0(nrow(M)))  # Add upper bound constraint
    f$b = c(f$b, rep(0, nrow(M)))  # Right-hand side for the constraints
    A = rbind(A, cbind(diag(1, nrow(M)), -gamma.max, 0))  # Enforce gamma ≤ gamma.max
  }
  
  # Solve using POGS optimizer
  pogs.solution = tryCatch({
    pogs(A, f, g, params = list(rel_tol=1e-4, abs_tol=1e-5, verbose=2*as.numeric(verbose)))
  }, error = function(e) {
    cat("Warning: POGS failed, returning zero vector.\n")
    return(list(x = rep(0, nrow(M))))
  })
  
  # Extract the gamma solution and normalize it to sum to 1
  gamma = pogs.solution$x[1:nrow(M)]
  gamma / sum(gamma)
}
#################################################################################################################
# Find approximately balancing weights using quadprog
approx.balance.quadprog = function(M,
                                   balance.target,
                                   zeta = 0.5,
                                   allow.negative.weights = FALSE) {
  # The system is effectively
  # minimize zeta * delta^2 + (1 - zeta) * ||gamma||^2
  # subject to
  #   sum gamma = 1
  #   delta + (M'gamma)_j >= balance.target_j
  #   delta - (M'gamma)_j >= -balance.target_j
  #
  # The last two constraints mean that
  # delta = ||M'gamma - balance.target||_infty
  
  Dmat = diag(c(zeta, rep(1 - zeta, nrow(M))))
  dvec = rep(0, 1 + nrow(M))
  Amat = cbind(
    c(0, rep(1, nrow(M))),
    rbind(rep(1, ncol(M)), M ),
    rbind(rep(1, ncol(M)), -M))
  bvec = c(1, balance.target, -balance.target)
  
  if (!allow.negative.weights) {
    LB = 1/nrow(M)/10000
    Amat = cbind(Amat, rbind(rep(0, nrow(M)), diag(rep(1, nrow(M)))))
    bvec = c(bvec, rep(LB, nrow(M)))
  }
  
  balance.soln = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  gamma = balance.soln$solution[-1]
  gamma
}

# Find approximately balancing weights using mosek, using dual of QP
approx.balance.mosek.dual = function(M,
                                     balance.target,
                                     zeta = 0.5,
                                     allow.negative.weights = FALSE,
                                     bound.gamma = FALSE,
                                     gamma.max = 1/nrow(M)^(2/3),
                                     verbose = FALSE) {
  
  # The primal problem is:
  # Minimize 1/2 x' diag(qvec) x
  # subject to Ax <= b,
  # where the constraints indexted by "equality" are required to be equalities
  #
  # Here, the vector x is interpreted as (max.imbalance, gamma)
  qvec.primal = 2 * c(zeta, rep(1 - zeta, nrow(M)))
  
  tM = Matrix::t(M)
  nvar = length(qvec.primal)
  
  A.primal.list = list(
    matrix(c(0, rep(1, nrow(M))), 1, nvar),
    cbind(rep(-1, ncol(M)), tM),
    cbind(rep(-1, ncol(M)), -tM),
    if(!allow.negative.weights) {
      cbind(0, Matrix::diag(-1, nrow(M), nrow(M)))
    } else {
      numeric()
    },
    if(bound.gamma) {
      cbind(0, Matrix::diag(1, nrow(M), nrow(M)))
    } else {
      numeric()
    }
  )
  A.primal = Reduce(rbind, A.primal.list)
  
  bvec = c(1,
           balance.target,
           -balance.target,
           if(!allow.negative.weights) {
             rep(0, nrow(M))
           } else {
             numeric()
           },
           if (bound.gamma) {
             rep(gamma.max, nrow(M))
           } else {
             numeric()
           })
  
  equality.primal = c(TRUE, rep(FALSE, 2*ncol(M) + (as.numeric(!allow.negative.weights) + as.numeric(bound.gamma)) * nrow(M)))
  
  # The dual problem is then
  # Minimize 1/2 lambda A diag(1/qvec) A' lambda + b * lambda
  # Subject to lambda >= 0 for the lambdas corresponding to inequality constraints
  
  #A.dual = diag(1/sqrt(qvec.primal)) %*% t(A.primal)
  A.dual = 1/sqrt(qvec.primal) * Matrix::t(A.primal) #this is the same thing, but faster
  
  # Next we turn this into a conic program:
  # Minimize t
  # Subject to bvec * lambda + q - t = 0
  #	A.dual * lambda - mu = 0
  #	lambda >= 0 (for the inequality constraints)
  #	2 * q >= ||mu||_2^2
  # Where we note that the last constraint is a rotated cone.
  #
  # Below, the solution vector is (lambda, mu, q, t, ONE), where ONE
  # is just a variable constrained to be 1.
  A.conic = rbind(c(bvec, rep(0, nvar), 1, -1, 0),
                  cbind(A.dual, Matrix::diag(-1, nvar, nvar), 0, 0, 0),
                  c(rep(0, length(bvec) + nvar + 2), 1))
  rhs.conic = c(rep(0, 1 + nvar), 1)
  
  blx.conic = rep(-Inf, ncol(A.conic))
  blx.conic[which(!equality.primal)] = 0
  bux.conic = rep(Inf, ncol(A.conic))
  
  obj.conic = c(rep(0, length(bvec) + nvar + 1), 1, 0)
  
  mosek.problem <- list()
  mosek.problem$sense <- "min"
  mosek.problem$c <- obj.conic
  mosek.problem$bx <- rbind(blx = blx.conic, bux = bux.conic)
  mosek.problem$A <- as(A.conic, "CsparseMatrix")
  mosek.problem$bc <- rbind(blc = rhs.conic, buc = rhs.conic)
  mosek.problem$cones <- cbind(list("RQUAD", c(length(bvec) + nvar + 1, length(bvec) + nvar + 3, length(bvec) + 1:nvar)))
  
  if (verbose) {
    mosek.out = Rmosek::mosek(mosek.problem)
  } else {
    mosek.out = Rmosek::mosek(mosek.problem, opts=list(verbose=0))
  }
  
  primal = -1/qvec.primal * (t(A.primal) %*% mosek.out$sol$itr$xx[1:nrow(A.primal)])
  delta = primal[1]
  gamma = primal[1 + 1:nrow(M)]
  gamma/sum(gamma)
}


# Find approximately balancing weights using pogs
approx.balance.pogs = function(M,
                               balance.target,
                               zeta = 0.5,
                               allow.negative.weights = FALSE,
                               bound.gamma = FALSE,
                               gamma.max = 1/nrow(M)^(2/3),
                               verbose = FALSE) {
  
  # Our original problem is the following:
  #
  # Minimize ||gamma||_2^2 + delta^2 / lambda^2, subject to
  # ||M'gamma - balance.target||_infty <= delta and sum gamma_i = 1.
  #
  # Equivalently, in notation recognized by POGS, we can write
  #
  # Minimize ||gamma||_2^2 + Z^2
  # subject to I[1:2p] <= 0, J[2p+1, 2p+2] = 1, where
  #
  #     (M'  -v -lambda)   (gamma)
  # I = (-M' v  -lambda) * (ONE  )
  #     (1'  0  0      )   (Z    )
  #     (0   1  0      )
  #
  # and v denotes balance.target.
  
  lambda = sqrt((1 - zeta) / zeta)
  g = list(h = kSquare())
  
  f = list(h = c(kIndLe0(2 * ncol(M)), kIndEq0(2)),
           b = c(rep(0, 2 * ncol(M)), 1, 1))
  A = rbind(cbind(t(M), -balance.target, -lambda),
            cbind(-t(M), balance.target, -lambda),
            c(rep(1, nrow(M)), 0, 0),
            c(rep(0, nrow(M)), 1, 0))
  
  if (!allow.negative.weights) {
    f$h = c(f$h, kIndLe0(nrow(M)))
    f$b = c(f$b, rep(0, nrow(M)))
    A = rbind(A, cbind(diag(-1, nrow(M)), 0, 0))
  }
  
  if (bound.gamma) {
    f$h = c(f$h, kIndLe0(nrow(M)))
    f$b = c(f$b, rep(0, nrow(M)))
    A = rbind(A, cbind(diag(1, nrow(M)), -gamma.max, 0))
  }
  
  pogs.solution = pogs(A, f, g, params = list(rel_tol=1e-4, abs_tol=1e-5, verbose=2*as.numeric(verbose)))
  
  gamma = pogs.solution$x[1:nrow(M)]
  gamma/sum(gamma)
}

# Find approximately balancing weights using pogs
approx.balance.pogs.dual = function(M,
                                    balance.target,
                                    zeta = 0.5,
                                    allow.negative.weights = FALSE,
                                    verbose = FALSE) {
  
  # The primal problem is:
  # Minimize 1/2 x' diag(qvec) x
  # subject to Ax <= b,
  # where the constraints indexted by "equality" are required to be equalities
  #
  # Here, the vector x is interpreted as (max.imbalance, gamma)
  qvec.primal = 2 * c(zeta, rep(1 - zeta, nrow(M)))
  A.primal = rbind(
    c(0, rep(1, nrow(M))),
    cbind(rep(-1, ncol(M)), t(M) ),
    cbind(rep(-1, ncol(M)), -t(M)))
  b.primal = c(1, balance.target, -balance.target)
  equality.primal = c(TRUE, rep(FALSE, 2*ncol(M)))
  
  if (!allow.negative.weights) {
    A.primal = rbind(A.primal, cbind(0, diag(-1, nrow(M), nrow(M))))
    b.primal = c(b.primal, rep(0, nrow(M)))
    equality.primal = c(equality.primal, rep(FALSE, nrow(M)))
  }
  
  # The dual problem is then
  # Minimize 1/2 lambda A diag(1/qvec) A' lambda + b * lambda
  # Subject to lambda >= 0 for the lambdas corresponding to inequality constraint
  #
  # This is equivalent to
  # Minimize 1/2 ||c||_2^2 + d
  # Subject to c = diag(1/sqrt(qvec)) A' lambda
  #  d = b * lambda
  #  lambda >= 0
  
  A.pogs = rbind(diag(1/sqrt(qvec.primal)) %*% t(A.primal), b.primal)
  f.pogs = list(h = c(kSquare(ncol(A.primal)), kIdentity(1)))
  g.pogs = list(h = kIndGe0(length(equality.primal)))
  g.pogs$h[equality.primal] = kZero(sum(equality.primal))
  
  
  pogs.solution = pogs(A.pogs, f.pogs, g.pogs, params = list(rel_tol=1e-4, abs_tol=1e-5, verbose=2*as.numeric(verbose)))
  
  primal = -diag(1/sqrt(qvec.primal)) %*% t(A.primal) %*% pogs.solution$x
  delta = primal[1]
  gamma = primal[1 + 1:nrow(M)]
  gamma / sum(gamma)
}