#' Parametric bootstrap for BayBiMR 
#'
#' @inheritParams BayBiMR
#' @param B Number of bootstrap replicates.
#' @param progress Show a text progress bar.
#' @return A list with bootstrap `summary`, replicate `draws`, any `errors`, and `meta`.
#' @export
BayBiMR_DP <- function(
  bX, sX, bY, sY,
  B = 200,
  n_iter = 6000, burnin = 2000, thin = 2,
  update_cov_every = 5,
  seed = 123,
  progress = TRUE
) {
  stopifnot(is.numeric(bX), is.numeric(sX), is.numeric(bY), is.numeric(sY))
  p <- length(bX)
  if (!all(lengths(list(sX, bY, sY)) == p)) stop("bX, sX, bY, sY must have same length.")
  if (any(sX <= 0) || any(sY <= 0)) stop("All sX and sY must be positive.")
  set.seed(seed)

  extract_two_estimates <- function(df) {
    stopifnot(is.data.frame(df))
    has_dir <- "direction" %in% names(df)
    if (has_dir) {
      ord <- match(c("X_to_Y","Y_to_X"), df$direction)
    } else {
      ord <- c(1L, 2L)
    }
    bad_names <- c("sd","se","ci_lo","ci_hi","post_P_gt0","post_P_lt0",
                   "postPgt0","postPlt0","T_used")
    num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    num_cols <- setdiff(num_cols, bad_names)
    cand_cols <- union(intersect("estimate", num_cols), num_cols)

    val <- c(NA_real_, NA_real_)
    for (col in cand_cols) {
      tmp <- suppressWarnings(as.numeric(df[[col]]))
      if (all(is.finite(tmp[ord]), na.rm = TRUE)) {
        val <- as.numeric(tmp[ord][1:2])
        break
      }
    }
    if (any(!is.finite(val))) {
      for (col in num_cols) {
        tmp <- suppressWarnings(as.numeric(df[[col]]))
        if (sum(is.finite(tmp)) >= 2) {
          idx <- which(is.finite(tmp))[1:2]
          val <- as.numeric(tmp[idx])
          break
        }
      }
    }
    if (any(!is.finite(val))) stop("Could not locate numeric estimate column.")
    names(val) <- c("theta_XY","theta_YX")
    val
  }

  one_rep <- function(b_idx) {
    bX_star <- bX + stats::rnorm(p, 0, sX)
    bY_star <- bY + stats::rnorm(p, 0, sY)
    BayBiMR(
      bX = bX_star, sX = sX, bY = bY_star, sY = sY,
      n_iter = n_iter, burnin = burnin, thin = thin,
      update_cov_every = update_cov_every, seed = seed + b_idx
    )
  }

  if (progress) pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
  on.exit(if (progress) close(pb), add = TRUE)

  est_XY <- rep(NA_real_, B); est_YX <- rep(NA_real_, B)
  ok <- logical(B)
  err_msg <- character(B)
  first_good_str <- NULL

  for (b in seq_len(B)) {
    res <- tryCatch(one_rep(b), error = identity)
    if (inherits(res, "error")) {
      err_msg[b] <- conditionMessage(res); ok[b] <- FALSE
    } else {
      if (is.null(first_good_str)) first_good_str <- capture.output(str(res))
      got <- tryCatch(extract_two_estimates(res), error = identity)
      if (inherits(got, "error")) {
        err_msg[b] <- paste("Parse error:", conditionMessage(got),
                            "Returned cols=", paste(names(res), collapse=", "))
        ok[b] <- FALSE
      } else {
        est_XY[b] <- got["theta_XY"]; est_YX[b] <- got["theta_YX"]
        ok[b] <- TRUE
      }
    }
    if (progress) utils::setTxtProgressBar(pb, b)
  }

  T_used <- sum(ok)
  if (T_used == 0) {
    top_err <- sort(table(err_msg[nchar(err_msg) > 0]), decreasing = TRUE)
    msg <- if (length(top_err)) names(top_err)[1] else "Unknown error."
    stop("All bootstrap replicates failed. Example error: ", msg,
         if (!is.null(first_good_str)) "\nExample successful result structure:\n" else "",
         if (!is.null(first_good_str)) paste(first_good_str, collapse = "\n") else "")
  }

  boot_mean_XY <- mean(est_XY[ok])
  boot_mean_YX <- mean(est_YX[ok])
  boot_se_XY   <- stats::sd(est_XY[ok])
  boot_se_YX   <- stats::sd(est_YX[ok])
  ci_XY        <- stats::quantile(est_XY[ok], c(0.025, 0.975), names = FALSE)
  ci_YX        <- stats::quantile(est_YX[ok], c(0.025, 0.975), names = FALSE)

  summary <- data.frame(
    direction = c("theta_XY","theta_YX"),
    boot_mean = c(boot_mean_XY, boot_mean_YX),
    boot_se   = c(boot_se_XY,   boot_se_YX),
    ci_lo     = c(ci_XY[1],     ci_YX[1]),
    ci_hi     = c(ci_XY[2],     ci_YX[2]),
    T_used    = T_used,
    stringsAsFactors = FALSE
  )

  out <- list(
    summary = summary,
    draws   = data.frame(replicate = which(ok),
                         theta_XY = est_XY[ok],
                         theta_YX = est_YX[ok]),
    errors  = data.frame(replicate = which(!ok), message = err_msg[!ok],
                         row.names = NULL),
    meta    = list(B = B, T_used = T_used, seed = seed,
                   n_iter = n_iter, burnin = burnin, thin = thin,
                   update_cov_every = update_cov_every)
  )
  class(out) <- c("BayBiMR_DP", class(out))
  out
}
