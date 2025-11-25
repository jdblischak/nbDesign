#' Simulate Recurrent Events with Fixed Follow-up
#'
#' Simulates recurrent events for a clinical trial with piecewise constant enrollment,
#' exponential failure rates (Poisson process), and piecewise exponential dropout.
#'
#' @param enroll_rate A data frame with columns \code{rate} and \code{duration} defining
#'   the piecewise constant enrollment rates.
#' @param fail_rate A data frame with columns \code{treatment} and \code{rate} defining
#'   the exponential failure rate for each treatment group.
#' @param dropout_rate A data frame with columns \code{treatment}, \code{rate}, and \code{duration}
#'   defining the piecewise constant dropout rates.
#' @param max_followup Numeric. Maximum duration of follow-up for each individual
#'   (relative to their randomization time).
#' @param n Total sample size. If NULL, it is estimated from \code{enroll_rate}.
#'   If provided, enrollment stops when \code{n} subjects are recruited.
#'
#' @return A data frame (tibble) with columns:
#'   \describe{
#'     \item{id}{Subject identifier}
#'     \item{treatment}{Treatment group}
#'     \item{enroll_time}{Time of enrollment relative to trial start}
#'     \item{tte}{Time to event or censoring relative to randomization}
#'     \item{calendar_time}{Calendar time of event or censoring (enroll_time + tte)}
#'     \item{event}{Binary indicator: 1 for event, 0 for censoring}
#'   }
#'   Multiple rows per subject are returned (one for each event, plus one for the final censoring time).
#'
#' @export
#'
#' @importFrom stats rexp runif
nb_sim <- function(enroll_rate, fail_rate, dropout_rate = NULL, max_followup = NULL, n = NULL) {
  # 1. Generate Enrollment
  # Simplified implementation of piecewise constant enrollment
  # If n is provided, we simulate until n. If not, we assume enroll_rate defines the full period.

  # Validate inputs
  if (is.null(enroll_rate) || !is.data.frame(enroll_rate)) stop("enroll_rate must be a data frame")
  if (is.null(fail_rate) || !is.data.frame(fail_rate)) stop("fail_rate must be a data frame")
  if (is.null(max_followup)) stop("max_followup must be provided")

  # Calculate total duration and expected N if n not provided
  if (is.null(n)) {
    n <- sum(enroll_rate$rate * enroll_rate$duration)
    n <- round(n) # Integer n
  }

  # Generate enrollment times
  # We need to sample n times from the piecewise density defined by enroll_rate
  total_enroll_duration <- sum(enroll_rate$duration)

  # Create a cumulative distribution function for enrollment time
  # Rates: r_1 for d_1, r_2 for d_2...
  # Mass in period i: m_i = r_i * d_i
  # Total mass M = sum(m_i). Note: This M might not equal n if n was forced.
  # We treat the piecewise constant function as a PDF (scaled).

  mass <- enroll_rate$rate * enroll_rate$duration
  cum_mass <- c(0, cumsum(mass))
  total_mass <- tail(cum_mass, 1)

  # Assign random uniform values scaled to total mass
  u_enroll <- runif(n, 0, total_mass)

  # Find which interval each u falls into
  # This can be vectorized with cut or findInterval
  interval_idx <- findInterval(u_enroll, cum_mass, rightmost.closed = TRUE)

  # For a point u in interval i (starting at cum_mass[i]), the time within interval is:
  # (u - cum_mass[i]) / rate[i]
  # Total time = sum(duration[1..i-1]) + time_within

  cum_duration <- c(0, cumsum(enroll_rate$duration))

  enroll_times <- numeric(n)

  # Vectorized calculation might be tricky with differing rates, so we iterate or use vectorized math
  rates_vec <- enroll_rate$rate[interval_idx]
  start_mass_vec <- cum_mass[interval_idx]
  start_time_vec <- cum_duration[interval_idx]

  enroll_times <- start_time_vec + (u_enroll - start_mass_vec) / rates_vec
  enroll_times <- sort(enroll_times) # Usually ids are ordered by enrollment

  # Assign treatments
  # Assuming simple randomization based on fail_rate rows?
  # sim_pw_surv usually takes n as argument but treatment allocation is often separate or implied.
  # Here we'll assume equal allocation or derived from fail_rate rows if multiple treatments.
  treatments <- unique(fail_rate$treatment)
  n_arms <- length(treatments)

  # Simple block randomization or just random assignment?
  # User prompt didn't specify, but "multiple treatment groups" implies we need to assign them.
  # We'll use simple random assignment for now.
  assigned_trt <- sample(treatments, n, replace = TRUE)

  # 2. Generate Event and Censoring Times

  out_list <- vector("list", n)

  for (i in 1:n) {
    trt <- assigned_trt[i]
    ent <- enroll_times[i]

    # Get failure rate for this treatment
    lambda <- fail_rate$rate[fail_rate$treatment == trt]
    if (length(lambda) == 0) stop(paste("No failure rate found for treatment", trt))
    if (length(lambda) > 1) lambda <- lambda[1] # Take first if multiple (shouldn't happen if properly formatted)

    # Dropout Time
    # Piecewise exponential dropout
    # If dropout_rate is NULL or rate 0, dropout_time is Inf
    dropout_time <- Inf
    if (!is.null(dropout_rate)) {
      # Filter dropout rates for this treatment (if specific) or use general
      dr <- dropout_rate
      if ("treatment" %in% names(dr)) {
        dr <- dr[dr$treatment == trt, ]
      }

      # Generate piecewise exponential time
      # Algorithm:
      # 1. Start at t=0.
      # 2. For each interval, generate exponential.
      # 3. If time < duration, event occurs. Else, subtract duration and move to next interval.

      t_curr <- 0
      generated <- FALSE
      for (j in 1:nrow(dr)) {
        rate_j <- dr$rate[j]
        dur_j <- dr$duration[j]

        if (rate_j > 0) {
          e <- rexp(1, rate_j)
          if (e < dur_j) {
            dropout_time <- t_curr + e
            generated <- TRUE
            break
          }
        }
        # Survived this interval
        t_curr <- t_curr + dur_j
      }

      # If survived all defined intervals, check if there is a final infinite interval implied?
      # simtrial usually implies last interval extends to infinity if not handled.
      # If not generated yet, and we ran out of intervals, does it continue with last rate?
      # Usually sim_pw_surv assumes inputs cover the duration.
      # If not generated, leave as Inf? Or assume last rate continues?
      # Let's assume last rate continues if not generated, akin to simtrial behavior often.
      if (!generated) {
        last_rate <- tail(dr$rate, 1)
        if (last_rate > 0) {
          dropout_time <- t_curr + rexp(1, last_rate)
        }
      }
    }

    # End of Follow-up for this subject
    # min(dropout, max_followup)
    # Note: Trial cutoff is NOT applied here per instructions.

    end_time <- min(dropout_time, max_followup)

    # Generate Recurrent Events (Poisson Process)
    # Inter-arrival times ~ Exp(lambda)

    event_times <- numeric(0)
    cum_t <- 0

    while (TRUE) {
      gap <- rexp(1, lambda)
      cum_t <- cum_t + gap

      if (cum_t <= end_time) {
        event_times <- c(event_times, cum_t)
      } else {
        break
      }
    }

    # Construct records
    # One row per event + one row for censoring (end_time)

    n_events <- length(event_times)

    ids <- rep(i, n_events + 1)
    trts <- rep(trt, n_events + 1)
    ents <- rep(ent, n_events + 1)
    ttes <- c(event_times, end_time)
    evts <- c(rep(1, n_events), 0) # 1 for event, 0 for censor

    # Calculate calendar time
    # calendar_time = enroll_time + tte
    cal_times <- ents + ttes

    out_list[[i]] <- data.frame(
      id = ids,
      treatment = trts,
      enroll_time = ents,
      tte = ttes,
      calendar_time = cal_times,
      event = evts
    )
  }

  # Combine all
  do.call(rbind, out_list)
}
