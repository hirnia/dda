# Internal helper: summarise a numeric vector with the chosen aggregation
# statistic. Shared by dda.bagging, reaggregate_bagging and summary_ols.
#' @noRd
dda_agg <- function(x, agg_stat, trim_prob = 0.10, win_prob = 0.10) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)

  switch(agg_stat,
         "mean"       = mean(x),
         "median"     = median(x),
         "trimmed"    = mean(x, trim = trim_prob),
         "winsorized" = {
           q_low  <- quantile(x, probs = win_prob,     na.rm = TRUE, names = FALSE)
           q_high <- quantile(x, probs = 1 - win_prob, na.rm = TRUE, names = FALSE)
           x[x < q_low]  <- q_low
           x[x > q_high] <- q_high
           mean(x)
         },
         "midhinge" = {
           q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
           mean(q)
         },
         "tukey" = {
           q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
           (q[1] + 2 * q[2] + q[3]) / 4
         }
  )
}
