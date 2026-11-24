# Run from the repository root: Rscript replicate.R [archived.dta] [output_directory]
if (!requireNamespace("readstata13", quietly = TRUE)) stop("Install the R package readstata13 first.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Expected at most an input file and an output directory.")
out <- if (length(args) > 1L) args[2] else "."
dir.create(out, recursive = TRUE, showWarnings = FALSE)
input <- if (length(args)) args[1] else file.path(tempdir(), "lebsurv_petition_pf_REP.dta")
if (!length(args) && !file.exists(input)) {
  options(timeout = max(180, getOption("timeout")))
  download.file("https://dataverse.harvard.edu/api/access/datafile/3124574", input,
                mode = "wb", method = "libcurl", headers = c("User-Agent" = "Mozilla/5.0"))
}
stopifnot(unname(tools::md5sum(input)) == "352c4160ae1fd6a95536e746bdb8dcb7")
items <- c("q42_b_sanc_zaim", "q42_d_sanc_fam", "q42_e_sanc_sect")
fields <- c("m_impute", "_mi_id", "id_psu", "strata_exp", "T_pet_pub", "pet_sign_final",
            "female", paste0(items, "_REC1b"), paste0(items, "_REC2"))
d <- readstata13::read.dta13(input, select.cols = fields, convert.factors = FALSE)
d <- as.data.frame(d[, fields])
names(d) <- c("imputation", "id", "psu", "block", "public", "signed", "female",
              "elite_item", "family_item", "community_item", "elite", "family", "community")
d <- d[order(d$imputation, d$id), ]
samples <- split(d, d$imputation)
observed <- samples[["0"]]
completed <- samples[as.character(1:10)]
fixed <- c("id", "psu", "block", "public", "signed", "female")
stopifnot(length(samples) == 11L, !anyDuplicated(observed$id),
          nrow(observed) == 2496L, all(observed$signed %in% 0:1), all(observed$public %in% 0:1))
for (s in completed) {
  stopifnot(!anyNA(s), all(vapply(fixed, function(k) identical(s[[k]], observed[[k]]), logical(1))))
}

# Equal-weight version of the authors' Anderson index: solve the correlation
# system for the weights of standardized items, then split at the within-imputation median.
fear_index <- function(s) {
  z <- scale(as.matrix(s[c("elite_item", "family_item", "community_item")]))
  w <- solve(cor(z), rep(1, ncol(z)))
  v <- drop(z %*% (w / sum(w)))
  as.numeric(v > median(v))
}

# Absorb block effects as in Stata areg. Center block effects to have mean zero
# over respondents, recovering the reported grand-mean-normalized intercept.
# Include the absorbed block degrees of freedom in Stata's CR1 correction.
cluster_fit <- function(y, x, block, psu) {
  block <- factor(block)
  demean <- function(v) v - ave(v, block, FUN = mean) + mean(v)
  X <- cbind(constant = 1, apply(x, 2, demean))
  colnames(X) <- c("constant", colnames(x))
  fit <- lm.fit(X, demean(y))
  bread <- chol2inv(chol(crossprod(X)))
  scores <- rowsum(X * as.numeric(fit$residuals), psu, reorder = FALSE)
  n <- nrow(X)
  g <- nrow(scores)
  k <- ncol(X) + nlevels(block) - 1L
  correction <- (g / (g - 1)) * ((n - 1) / (n - k))
  vc <- correction * bread %*% crossprod(scores) %*% bread
  list(b = fit$coefficients, variance = diag(vc))
}
fit_column <- function(s, model) {
  x <- cbind(public = s$public)
  if (model != "main") {
    interaction <- s$public * s[[model]]
    x <- cbind(x, interaction = interaction)
    controls <- if (model == "index") "index" else c("elite", "family", "community")
    x <- cbind(x, as.matrix(s[controls]))
  }
  cluster_fit(s$signed, x, s$block, s$psu)
}
completed <- lapply(completed, function(s) { s$index <- fear_index(s); s })
models <- c("main", "index", "elite", "family", "community")
# Published estimates and SEs at the paper's two-decimal precision. Within each
# column the order is constant, public, interaction (if any), then controls.
published_b <- list(c(.70,-.20), c(.70,-.16,-.08,.01), c(.70,-.17,-.10,.05,-.03,-.01),
                    c(.69,-.14,-.14,.00,.04,-.01), c(.70,-.16,-.09,.00,-.03,.03))
published_se <- list(c(.01,.02), c(.02,.03,.04,.03), c(.02,.03,.04,.03,.02,.03),
                     c(.02,.03,.04,.02,.03,.03), c(.02,.03,.04,.02,.02,.03))
rows <- lapply(seq_along(models), function(i) {
  fits <- lapply(completed, fit_column, model = models[i])
  b <- do.call(cbind, lapply(fits, "[[", "b"))
  within <- rowMeans(do.call(cbind, lapply(fits, "[[", "variance")))
  between <- apply(b, 1, var)
  se <- sqrt(within + (1 + 1/ncol(b)) * between)
  data.frame(column = i, model = models[i], term = rownames(b),
             estimate = unname(rowMeans(b)), se = unname(se), published = published_b[[i]],
             published_se = published_se[[i]], n = nrow(observed), imputations = ncol(b))
})
table1 <- do.call(rbind, rows)
table1$match_estimate <- round(table1$estimate, 2) == table1$published
table1$match_se <- round(table1$se, 2) == table1$published_se
write.csv(table1, file.path(out, "table-1.csv"), row.names = FALSE)
arms <- aggregate(cbind(signed = observed$signed, assigned = rep(1, nrow(observed))),
                  list(public = observed$public), sum)
arms$proportion <- arms$signed / arms$assigned
write.csv(arms, file.path(out, "signing-rates.csv"), row.names = FALSE)
source("randomization.R", local = TRUE)
tests <- randomization_tests(observed)
write.csv(tests, file.path(out, "randomization-tests.csv"), row.names = FALSE)
main <- subset(table1, column == 1 & term == "public")
stopifnot(isTRUE(all.equal(main$estimate, diff(arms$proportion), check.attributes = FALSE)),
          isTRUE(all.equal(main$estimate, tests$estimate_pp[1] / 100, tolerance = 1e-10)))

# Show the observed action separately from the tests of heterogeneous effects.
# Intervals use a normal approximation to the clustered, imputation-pooled SE.
png(file.path(out, "petition-results.png"), width = 1530, height = 830, res = 145)
par(mar = c(4.3, 12.4, 5.4, 4.6), bg = "#faf8f3", fg = "#263744", family = "sans",
    col.axis = "#485962", col.lab = "#263744", lend = "round")
rr <- rbind(main, subset(table1, term == "interaction"))
lo <- 100*(rr$estimate - 1.96*rr$se)
hi <- 100*(rr$estimate + 1.96*rr$se)
yy <- c(5.6,3.9,2.9,1.9,0.9)
plot(0, type = "n", xlim = c(-27,4), ylim = c(0.1,6.1), xaxs = "i",
     axes = FALSE, xlab = "", ylab = "")
abline(v = c(-25,-20,-15,-10,-5,0), col = "#dde1df", lwd = 0.8)
abline(v = 0, col = "#77878a", lwd = 1.3)
abline(h = 4.85, col = "#becac8")
axis(1, at = seq(-25,0,5), lwd = 0, lwd.ticks = 1, tck = -0.018, mgp = c(2,0.5,0))
axis(2, at = yy, labels = c("All respondents", "Overall fear index", "Political or sectarian elites",
                           "Family, friends, neighbors", "Confessional community"),
     tick = FALSE, las = 1, line = 0.2, cex.axis = 0.87)
cols <- c("#a15538", rep("#24566b", 4))
segments(lo, yy, hi, yy, lwd = 2.2, col = cols)
points(100*rr$estimate, yy, pch = c(18,19,19,19,19), cex = c(1.5,1,1,1,1), col = cols)
text(3, yy, sprintf("%+.1f", 100*rr$estimate), adj = 0, xpd = NA, cex = 0.85, col = cols)
text(-27, 4.52, "Additional effect among those reporting higher fear", pos = 4, offset = 0,
     col = "#536368", cex = 0.77)
mtext("Difference in petition signing (percentage points)", 1, line = 2.25, cex = 0.9)
mtext("What changes when a petition requires a name?", 3, line = 3.0, adj = 0, cex = 1.35, font = 2)
mtext("Lebanon  |  2,496 participants  |  Table 1 estimates and approximate 95% intervals",
      3, line = 1.25, adj = 0, cex = 0.79, col = "#485962")
invisible(dev.off())
print(table1, digits = 4, row.names = FALSE)
print(tests, digits = 5, row.names = FALSE)
