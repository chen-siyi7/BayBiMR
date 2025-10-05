#' Bidirectional BayBiMR sampler with latent updates in C++
#'
#' @param bX,bY Numeric vectors of marginal SNP associations.
#' @param sX,sY Numeric vectors of standard errors (positive).
#' @param n_iter,burnin,thin Sampler controls.
#' @param tau_theta Prior SD for causal effects.
#' @param api Beta prior params for inclusion probs (length 2).
#' @param nu4_0,S4_0 IW prior for 4x4 slab covariance.
#' @param init List of initial values.
#' @param update_cov_every Refresh period for Sigma4.
#' @param return_samples If TRUE, stores latent draws.
#' @param seed RNG seed.
#' @return A data.frame with posterior summaries; PIPs, draws as attributes.
#' @export
BayBiMR <- function(
  bX, sX, bY, sY,
  n_iter = 6000, burnin = 2000, thin = 1L,
  tau_theta = 0.3,
  api = c(1, 19),
  nu4_0 = 6, S4_0 = diag(4)*0.05,
  init = list(
    theta_XY = 0.05, theta_YX = 0.00,
    pi_g = 0.2, pi_d = 0.2,
    Sigma4 = diag(c(0.08, 0.06, 0.08, 0.05))
  ),
  update_cov_every = 5L,
  return_samples = FALSE,
  seed = NULL
){
  if (!is.null(seed)) set.seed(seed)
  p <- length(bX)
  stopifnot(length(bY)==p, length(sX)==p, length(sY)==p)
  if (burnin >= n_iter) stop("'burnin' must be < n_iter")
  if (thin < 1L) thin <- 1L
  if (update_cov_every < 1L) update_cov_every <- 1L

  eps <- 1e-12
  sX2 <- pmax(sX^2, eps); sY2 <- pmax(sY^2, eps)
  uX  <- 1 / sX2;         uY  <- 1 / sY2
  tau_theta2_inv <- 1 / (tau_theta^2)
  log2pi <- log(2*pi)

  logpdf2_vec <- function(v1, v2, c11, c12, c22) {
    det <- pmax(c11*c22 - c12*c12, 1e-20)
    a <-  c22 / det
    b <- -c12 / det
    c <-  c11 / det
    q <- v1*(a*v1 + b*v2) + v2*(b*v1 + c*v2)
    -(log2pi + 0.5*log(det) + 0.5*q)
  }

  rIW <- function(df, S) {
    Wi <- stats::rWishart(1, df = df, Sigma = solve(S))[,,1]
    solve(Wi)
  }

  theta_XY <- init$theta_XY; theta_YX <- init$theta_YX
  pi_g <- min(max(init$pi_g, eps), 1 - eps)
  pi_d <- min(max(init$pi_d, eps), 1 - eps)

  Sigma4 <- init$Sigma4
  Sigma4 <- (Sigma4 + t(Sigma4)) / 2
  diag(Sigma4) <- pmax(diag(Sigma4), 1e-8)
  J4 <- solve(Sigma4)

  z <- rbinom(p, 1, pi_g)
  w <- rbinom(p, 1, pi_d)
  x <- bX; y <- bY
  gamma <- numeric(p); delta <- numeric(p)

  keep <- (burnin + 1):n_iter
  keep <- keep[seq(1, length(keep), by = thin)]
  n_keep <- length(keep)

  draws_theta <- matrix(NA_real_, nrow = n_keep, ncol = 2,
                        dimnames = list(NULL, c("theta_XY","theta_YX")))
  draws_pi <- matrix(NA_real_, nrow = n_keep, ncol = 2,
                     dimnames = list(NULL, c("pi_g","pi_d")))
  pip_z_sum <- numeric(p); pip_w_sum <- numeric(p)
  Sigma4_draws <- array(NA_real_, dim = c(n_keep, 4, 4))
  if (return_samples) {
    draws_x <- matrix(NA_real_, nrow = n_keep, ncol = p)
    draws_y <- matrix(NA_real_, nrow = n_keep, ncol = p)
    draws_g <- matrix(NA_real_, nrow = n_keep, ncol = p)
    draws_d <- matrix(NA_real_, nrow = n_keep, ncol = p)
  }

  k <- 0L
  for (it in 1:n_iter) {

    zw <- cbind(as.integer(z), as.integer(w))
    lat <- baymix_update_latents_cpp(bX, bY, uX, uY, theta_XY, theta_YX, zw, J4)
    x     <- lat$x
    gamma <- lat$gamma
    y     <- lat$y
    delta <- lat$delta

    den_XY <- sum(uY * x^2) + tau_theta2_inv
    vth <- 1 / den_XY
    mu  <- vth * sum(uY * x * (bY - gamma))
    theta_XY <- rnorm(1, mu, sqrt(vth))

    den_YX <- sum(uX * y^2) + tau_theta2_inv
    vth <- 1 / den_YX
    mu  <- vth * sum(uX * y * (bX - delta))
    theta_YX <- rnorm(1, mu, sqrt(vth))

    Sxg <- Sigma4[c(1,2), c(1,2)]
    Syd <- Sigma4[c(3,4), c(3,4)]
    S11 <- Sxg[1,1]; S12 <- Sxg[1,2]; S22 <- Sxg[2,2]
    C1_11 <- S11 + sX2
    C1_12 <- theta_XY * S11 + S12
    C1_22 <- theta_XY^2 * S11 + 2*theta_XY*S12 + S22 + sY2
    C0_11 <- S11 + sX2
    C0_12 <- theta_XY * S11
    C0_22 <- theta_XY^2 * S11 + sY2
    lp1 <- log(pmax(pi_g, eps))     + logpdf2_vec(bX, bY, C1_11, C1_12, C1_22)
    lp0 <- log(pmax(1 - pi_g, eps)) + logpdf2_vec(bX, bY, C0_11, C0_12, C0_22)
    z <- rbinom(p, 1, plogis(lp1 - lp0))

    S11 <- Syd[1,1]; S12 <- Syd[1,2]; S22 <- Syd[2,2]
    D1_11 <- theta_YX^2 * S11 + 2*theta_YX*S12 + S22 + sX2
    D1_12 <- theta_YX * S11 + S12
    D1_22 <- S11 + sY2
    E0_11 <- theta_YX^2 * S11 + sX2
    E0_12 <- theta_YX * S11
    E0_22 <- S11 + sY2
    lp1 <- log(pmax(pi_d, eps))     + logpdf2_vec(bX, bY, D1_11, D1_12, D1_22)
    lp0 <- log(pmax(1 - pi_d, eps)) + logpdf2_vec(bX, bY, E0_11, E0_12, E0_22)
    w <- rbinom(p, 1, plogis(lp1 - lp0))

    if (it %% update_cov_every == 0L) {
      active <- (z == 1L) | (w == 1L)
      if (any(active)) {
        Z <- cbind(x[active], gamma[active], y[active], delta[active])
        Sigma4 <- rIW(nu4_0 + nrow(Z), S4_0 + crossprod(Z))
      } else {
        Sigma4 <- rIW(nu4_0, S4_0)
      }
      Sigma4 <- (Sigma4 + t(Sigma4)) / 2
      diag(Sigma4) <- pmax(diag(Sigma4), 1e-10)
      J4 <- solve(Sigma4)
    }

    pi_g <- min(max(stats::rbeta(1, api[1] + sum(z), api[2] + p - sum(z)), eps), 1 - eps)
    pi_d <- min(max(stats::rbeta(1, api[1] + sum(w), api[2] + p - sum(w)), eps), 1 - eps)

    if (it %in% keep) {
      k <- k + 1L
      draws_theta[k,] <- c(theta_XY, theta_YX)
      draws_pi[k,]    <- c(pi_g, pi_d)
      pip_z_sum <- pip_z_sum + z
      pip_w_sum <- pip_w_sum + w
      Sigma4_draws[k,,] <- Sigma4
      if (return_samples) {
        draws_x[k,] <- x; draws_y[k,] <- y
        draws_g[k,] <- gamma; draws_d[k,] <- delta
      }
    }
  }

  post <- as.data.frame(draws_theta)
  est <- colMeans(post); sdv <- apply(post, 2, sd)
  ci  <- t(apply(post, 2, stats::quantile, probs = c(0.025, 0.975)))
  colnames(ci) <- c("ci_lo","ci_hi")

  res <- data.frame(
    direction   = c("X_to_Y","Y_to_X"),
    estimate    = as.numeric(est),
    sd          = as.numeric(sdv),
    ci_lo       = as.numeric(ci[,1]),
    ci_hi       = as.numeric(ci[,2]),
    post_P_gt0  = colMeans(post > 0),
    post_P_lt0  = colMeans(post < 0),
    mode        = "latent_corr4_cpp"
  )
  attr(res, "pi_draws")       <- draws_pi
  attr(res, "pip_z")          <- pip_z_sum / n_keep
  attr(res, "pip_w")          <- pip_w_sum / n_keep
  attr(res, "Sigma4_draws")   <- Sigma4_draws
  if (return_samples) {
    attr(res, "x_draws")     <- draws_x
    attr(res, "y_draws")     <- draws_y
    attr(res, "gamma_draws") <- draws_g
    attr(res, "delta_draws") <- draws_d
  }
  res
}
