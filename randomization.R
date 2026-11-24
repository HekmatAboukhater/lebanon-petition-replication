# Finite-distribution tests of the sharp null (no signing response changes under disclosure).
# Polynomial coefficient j+1 is the probability of j public-arm signers.
# Add nonnegative terms directly to preserve the small tail probabilities.
add_counts <- function(p, q) {
  result <- numeric(length(p) + length(q) - 1L)
  for (j in which(q > 0)) {
    idx <- seq_along(p) + j - 1L
    result[idx] <- result[idx] + p * q[j]
  }
  result
}
individual_counts <- function(rows) {
  n <- nrow(rows)
  signed <- sum(rows$signed)
  public <- sum(rows$public)
  dhyper(0:public, signed, n - signed, public)
}
cluster_counts <- function(rows) {
  units <- split(rows, rows$psu)
  stopifnot(all(vapply(units, function(u) length(unique(u$public)) == 1L, logical(1))),
            length(unique(vapply(units, nrow, integer(1)))) == 1L)
  totals <- vapply(units, function(u) sum(u$signed), numeric(1))
  treated <- sum(vapply(units, function(u) u$public[1], numeric(1)))
  stopifnot(treated > 0, treated < length(units))
  possible <- combn(totals, treated, FUN = sum)
  tabulate(possible + 1L, nbins = max(possible) + 1L) / length(possible)
}
randomization_tests <- function(obs) {
  stopifnot(all(tapply(obs$public, obs$block, mean) == 0.5))
  small <- startsWith(obs$block, "small_")
  by_psu <- split(which(small), obs$psu[small])
  mixed_small <- unlist(by_psu[vapply(by_psu, function(i) length(unique(obs$public[i])) > 1L,
                                      logical(1))], use.names = FALSE)
  irregular <- obs$block %in% unique(obs$block[mixed_small])
  stopifnot(sum(irregular) == 16, length(unique(obs$block[irregular])) == 1,
            all(obs$female %in% 0:1))
  one <- function(label, method, omit_irregular = FALSE, only_clusters = FALSE) {
    use <- rep(TRUE, nrow(obs))
    if (omit_irregular) use <- !irregular
    if (only_clusters) use <- small & !irregular
    x <- obs[use, ]
    freeze <- irregular[use]
    small_x <- small[use]
    cl <- x[small_x & !freeze, ]
    indiv <- x[!small_x, ]
    components <- lapply(split(cl, cl$block), cluster_counts)
    if (nrow(indiv)) {
      cells <- if (method == "sex") interaction(indiv$psu, indiv$female, drop = TRUE) else indiv$psu
      components <- c(components, lapply(split(indiv, cells), individual_counts))
    }
    dist <- Reduce(add_counts, components, init = 1)
    stopifnot(abs(sum(dist) - 1) < 1e-10)
    # For every allowed assignment each regression block stays exactly half public.
    # The block-adjusted public-minus-private coefficient is (2*S - total) / (N/2).
    frozen_public_signers <- sum(x$public[freeze] * x$signed[freeze])
    possible <- (seq_along(dist) - 1L) + frozen_public_signers
    observed_signers <- sum(x$public * x$signed)
    total <- sum(x$signed)
    observed_distance <- abs(2 * observed_signers - total)
    two_sided <- sum(dist[abs(2 * possible - total) >= observed_distance])
    data.frame(design = label, n = nrow(x), frozen_records = sum(freeze),
               individually_reassigned = nrow(indiv), clustered_records = nrow(cl),
               reassigned_clusters = length(unique(cl$psu)),
               estimate_pp = 100 * (2 * observed_signers - total) / (nrow(x)/2),
               sharp_null_p = two_sided)
  }
  rbind(one("All records, realized sex margins, ambiguous block fixed", "sex"),
        one("Omit ambiguous block, realized sex margins", "sex", omit_irregular = TRUE),
        one("All records, PSU margins only, ambiguous block fixed", "psu"),
        one("Unambiguous cluster-assigned blocks only", "sex", only_clusters = TRUE))
}
