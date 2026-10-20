# =============================================================================
# Constrained portfolio optimization: minimum-variance, index tracking and
# efficient frontier, with an out-of-sample evaluation
# =============================================================================
#
# Expected data (Excel file, simple MONTHLY returns, in decimals):
#   - Date  : month
#   - Bond  : monthly change in the US 10-year yield, in percentage points
#   - ETF   : return of a Nasdaq-100 ETF, used as the benchmark
#   - INDEX : ignored (unusable series, see README)
#   - all other columns: returns of the stocks in the universe
#
# Approach:
#   1. TRAINING period (first IS_MONTHS months): every estimate (covariance,
#      weights, frontier, solver comparison) is computed here.
#   2. TEST period (remaining months): portfolios are evaluated on data that
#      was never used for any estimation.
#
# Run from the project folder: Rscript portfolio_optimization.R
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Parameters
# -----------------------------------------------------------------------------

DATA_PATH   <- "returns.xlsx"        # path to the data file
DATA_SHEET  <- 1                     # sheet name or number
BENCH_COL   <- "ETF"                 # benchmark: the Nasdaq-100 ETF (the INDEX column
                                     # of the file is unusable, see README)
NON_ASSET   <- c("Date", "Bond", "INDEX", "ETF")  # columns that are not stocks
RF_UNIT     <- "yield_change"        # "yield_change" | "zero" | "monthly_decimal" | "annual_decimal" | "annual_percent"
# "yield_change": the Bond column holds the monthly CHANGE in the US 10-year yield, in
# percentage points (e.g. -0.414 in January 2020: 1.92% -> 1.51%). The yield level is
# rebuilt from its end-of-December-2019 value, then converted into a monthly rate.
RF_START_LEVEL <- 1.92               # US 10-year yield on 31 Dec 2019, in %
IS_MONTHS   <- 48                    # length of the training period (months)
LOWER_BOUND <- 0.0001                # minimum weight per stock (long-only)
UPPER_BOUND <- 0.20                  # maximum weight per stock
N_TARGETS   <- 15                    # number of points on the efficient frontier
PERIODS     <- 12                    # periods per year (monthly data)
OUTPUT_DIR  <- "outputs"

required <- c("readxl", "quadprog", "nloptr", "ggplot2")
missing  <- required[!required %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing)

library(readxl)
library(quadprog)
library(nloptr)
library(ggplot2)

dir.create(file.path(OUTPUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 1. Data loading and cleaning
# -----------------------------------------------------------------------------

to_numeric <- function(x) as.numeric(gsub(",", ".", as.character(x)))  # handles decimal commas

raw <- as.data.frame(read_excel(DATA_PATH, sheet = DATA_SHEET))
stopifnot(all(c("Date", "Bond", BENCH_COL) %in% names(raw)))

asset_cols <- setdiff(names(raw), NON_ASSET)
R_raw      <- suppressWarnings(sapply(raw[asset_cols], to_numeric))
# empty or non-numeric columns (e.g. "...104" columns created by Excel): dropped
empty_cols <- colSums(!is.na(R_raw)) == 0
if (any(empty_cols)) {
  message("Ignored columns (empty or non-numeric): ", paste(asset_cols[empty_cols], collapse = ", "))
  R_raw      <- R_raw[, !empty_cols, drop = FALSE]
  asset_cols <- asset_cols[!empty_cols]
}
# stocks with missing months: dropped, rather than deleting those months for the whole universe
incomplete <- colSums(is.na(R_raw)) > 0
if (any(incomplete)) {
  message("Stocks dropped (missing data): ", paste(asset_cols[incomplete], collapse = ", "))
  R_raw      <- R_raw[, !incomplete, drop = FALSE]
  asset_cols <- asset_cols[!incomplete]
}
rf_raw     <- to_numeric(raw$Bond)
bench_raw  <- to_numeric(raw[[BENCH_COL]])
dates_raw  <- as.Date(raw$Date)

# a single valid-row filter for every series (guaranteed alignment)
valid <- !is.na(dates_raw) & !is.na(rf_raw) & !is.na(bench_raw) & complete.cases(R_raw)
dates <- dates_raw[valid]
R     <- R_raw[valid, , drop = FALSE]
bench <- bench_raw[valid]
rf    <- switch(RF_UNIT,
  zero            = rep(0, sum(valid)),
  yield_change    = (1 + (RF_START_LEVEL + cumsum(rf_raw))[valid] / 100)^(1 / PERIODS) - 1,
  monthly_decimal = rf_raw[valid],
  annual_decimal  = (1 + rf_raw[valid])^(1 / PERIODS) - 1,
  annual_percent  = (1 + rf_raw[valid] / 100)^(1 / PERIODS) - 1,
  stop("Unknown RF_UNIT")
)

ord   <- order(dates)                  # chronological order
dates <- dates[ord]; R <- R[ord, , drop = FALSE]; bench <- bench[ord]; rf <- rf[ord]

n_obs    <- nrow(R)
n_assets <- ncol(R)

if (median(abs(R)) > 0.5) stop("Asset columns look like prices, not returns.")
if (IS_MONTHS >= n_obs - 1) stop("IS_MONTHS leaves fewer than 2 months for the test period.")
if (LOWER_BOUND * n_assets > 1 || UPPER_BOUND * n_assets < 1) stop("Infeasible weight bounds.")

is_idx  <- 1:IS_MONTHS
oos_idx <- (IS_MONTHS + 1):n_obs
lb <- rep(LOWER_BOUND, n_assets)
ub <- rep(UPPER_BOUND, n_assets)

cat(sprintf("%d stocks, %d months (%s -> %s)\n", n_assets, n_obs, min(dates), max(dates)))

# --- Data check: a benchmark consistent with the universe should be strongly
#     correlated with the equal-weighted portfolio of the stocks (typically > 0.8)
ew <- rowMeans(R)
cat("\n=== Data check ===\n")
for (col in intersect(c("INDEX", "ETF", "Bond"), names(raw))) {
  x <- to_numeric(raw[[col]])[valid][ord]
  cat(sprintf("%-6s | monthly mean %8.4f | std dev %7.4f | min %8.4f | max %8.4f | corr. with 1/N %6.2f\n",
              col, mean(x, na.rm = TRUE), sd(x, na.rm = TRUE), min(x, na.rm = TRUE), max(x, na.rm = TRUE),
              cor(x, ew, use = "complete.obs")))
}
cat(sprintf("1/N    | monthly mean %8.4f | std dev %7.4f\n", mean(ew), sd(ew)))
cat("First dates:", format(head(dates, 3)), "| last:", format(tail(dates, 3)), "\n")
if (cor(bench, ew) < 0.5) warning("Column ", BENCH_COL, " is weakly correlated with the stocks: check the data.")

cat(sprintf("Training: %d months | Test: %d months (from %s)\n\n",
            length(is_idx), length(oos_idx), dates[oos_idx[1]]))


# -----------------------------------------------------------------------------
# 2. Tools: covariance, solvers, metrics
# -----------------------------------------------------------------------------

# Ledoit-Wolf (2004) covariance, target = scaled identity matrix.
# With more stocks than months (e.g. 99 stocks, 48 months), the sample covariance
# is singular: the optimizer would exploit "risk-free" directions that are only
# estimation noise. Shrinkage fixes this.
ledoit_wolf <- function(X) {
  X <- scale(X, center = TRUE, scale = FALSE)
  n <- nrow(X); p <- ncol(X)
  S <- crossprod(X) / n
  target <- diag(sum(diag(S)) / p, p)
  d2 <- sum((S - target)^2)
  b2 <- sum(apply(X, 1, function(x) sum((tcrossprod(x) - S)^2))) / n^2
  delta <- min(b2, d2) / d2                      # shrinkage intensity, between 0 and 1
  Sigma <- delta * target + (1 - delta) * S
  dimnames(Sigma) <- list(colnames(X), colnames(X))
  attr(Sigma, "shrinkage") <- delta
  Sigma
}

valid_weights <- function(w, mu = NULL, target = NULL, tol = 1e-6) {
  ok <- abs(sum(w) - 1) < tol && all(w >= lb - tol) && all(w <= ub + tol)
  if (!is.null(target)) ok <- ok && abs(sum(w * mu) - target) < tol
  ok
}

# min w' Sigma w  s.t.  sum(w) = 1, lb <= w <= ub  [and mu'w = target if provided]
# Solver 1: exact quadratic programming (quadprog)
solve_qp <- function(Sigma, mu = NULL, target = NULL) {
  n    <- ncol(Sigma)
  Amat <- cbind(rep(1, n), diag(n), -diag(n))
  bvec <- c(1, lb, -ub)
  meq  <- 1
  if (!is.null(target)) {
    Amat <- cbind(rep(1, n), mu, diag(n), -diag(n))
    bvec <- c(1, target, lb, -ub)
    meq  <- 2
  }
  res <- tryCatch(
    solve.QP(Dmat = 2 * Sigma + diag(1e-10, n), dvec = rep(0, n),
             Amat = Amat, bvec = bvec, meq = meq),
    error = function(e) NULL
  )
  if (is.null(res) || !valid_weights(res$solution, mu, target)) return(rep(NA_real_, n))
  setNames(res$solution, colnames(Sigma))
}

# Solver 2: generic nonlinear optimization (nloptr, SLSQP), for comparison
solve_nlopt <- function(Sigma, mu = NULL, target = NULL) {
  n  <- ncol(Sigma)
  eq <- function(w) {
    if (is.null(target)) {
      list(constraints = sum(w) - 1, jacobian = matrix(1, 1, n))
    } else {
      list(constraints = c(sum(w) - 1, sum(w * mu) - target), jacobian = rbind(rep(1, n), mu))
    }
  }
  res <- nloptr(
    x0 = rep(1 / n, n),
    eval_f = function(w) list(objective = as.numeric(t(w) %*% Sigma %*% w),
                              gradient  = as.numeric(2 * Sigma %*% w)),
    lb = lb, ub = ub,
    eval_g_eq = eq,
    opts = list(algorithm = "NLOPT_LD_SLSQP", xtol_rel = 1e-10, ftol_rel = 1e-12, maxeval = 5000)
  )
  if (res$status < 0 || !valid_weights(res$solution, mu, target, tol = 1e-4)) return(rep(NA_real_, n))
  setNames(res$solution, colnames(Sigma))
}

# Tracking portfolio: min Var(R w - benchmark). Since sum(w) = 1,
# R w - benchmark = (R - benchmark) w: a minimum-variance on active returns.
tracking_weights <- function(R_win, bench_win) {
  solve_qp(ledoit_wolf(sweep(R_win, 1, bench_win, "-")))
}

# Highest (or lowest) return attainable under the weight bounds (greedy algorithm)
feasible_return <- function(mu, maximize = TRUE) {
  w <- lb
  remaining <- 1 - sum(w)
  for (i in order(mu, decreasing = maximize)) {
    add <- min(ub[i] - lb[i], remaining)
    w[i] <- w[i] + add
    remaining <- remaining - add
    if (remaining <= 1e-12) break
  }
  sum(w * mu)
}

perf_metrics <- function(r, rf, bench, periods = PERIODS) {
  ex       <- r - rf
  ex_b     <- bench - rf
  active   <- r - bench
  wealth   <- cumprod(c(1, 1 + r))
  max_dd   <- min(wealth / cummax(wealth) - 1)
  downside <- sqrt(mean(pmin(ex, 0)^2))
  beta     <- cov(ex, ex_b) / var(ex_b)
  q05      <- quantile(r, 0.05, names = FALSE)
  te       <- sd(active) * sqrt(periods)
  c(
    "Annualized return"          = mean(r) * periods,
    "Annualized volatility"      = sd(r) * sqrt(periods),
    "Sharpe ratio"               = mean(ex) / sd(ex) * sqrt(periods),
    "Sortino ratio"              = if (downside > 0) mean(ex) / downside * sqrt(periods) else NA,
    "Omega ratio (threshold 0)"  = sum(pmax(r, 0)) / sum(pmax(-r, 0)),
    "Beta"                       = beta,
    "Treynor ratio"              = mean(ex) * periods / beta,
    "Jensen's alpha (annual)"    = (mean(ex) - beta * mean(ex_b)) * periods,
    "Max drawdown"               = max_dd,
    "Calmar ratio"               = mean(r) * periods / abs(max_dd),
    "VaR 95% (monthly)"          = -q05,
    "ES 95% (monthly)"           = -mean(r[r <= q05]),
    "Tracking error"             = te,
    "Information ratio"          = if (te > 1e-12) mean(active) * periods / te else NA,
    "Correlation with benchmark" = cor(r, bench)
  )
}

metrics_table <- function(strategies, idx) {
  sapply(strategies, perf_metrics, rf = rf[idx], bench = bench[idx])
}


# -----------------------------------------------------------------------------
# 3. Training period: portfolios and efficient frontier
# -----------------------------------------------------------------------------

R_is     <- R[is_idx, , drop = FALSE]
Sigma_is <- ledoit_wolf(R_is)
mu_is    <- colMeans(R_is)
cat(sprintf("Ledoit-Wolf shrinkage intensity: %.2f\n", attr(Sigma_is, "shrinkage")))

w_gmv   <- solve_qp(Sigma_is)                 # global minimum-variance
w_track <- tracking_weights(R_is, bench[is_idx])
stopifnot(!anyNA(w_gmv), !anyNA(w_track))

strategies_is <- list(
  "Benchmark"           = bench[is_idx],
  "Equal weight (1/N)"  = rowMeans(R_is),
  "Minimum variance"    = as.numeric(R_is %*% w_gmv),
  "Index tracking"      = as.numeric(R_is %*% w_track)
)
table_is <- metrics_table(strategies_is, is_idx)
cat("\n=== Training period (in-sample: optimistic by construction) ===\n")
print(round(table_is, 4))
write.csv(round(table_is, 6), file.path(OUTPUT_DIR, "metrics_in_sample.csv"))

# --- Efficient frontier: from the minimum-variance portfolio to the highest attainable return
target_returns <- seq(sum(w_gmv * mu_is), feasible_return(mu_is, maximize = TRUE),
                      length.out = N_TARGETS)
# the last point (maximum return) has a single solution, often poorly resolved numerically
target_returns <- head(target_returns, -1)

W_qp    <- t(sapply(target_returns, function(m) solve_qp(Sigma_is, mu_is, m)))
W_nlopt <- t(sapply(target_returns, function(m) solve_nlopt(Sigma_is, mu_is, m)))
colnames(W_qp) <- colnames(W_nlopt) <- asset_cols

port_vol <- function(W) apply(W, 1, function(w) sqrt(as.numeric(t(w) %*% Sigma_is %*% w) * PERIODS))

solver_comparison <- data.frame(
  target_return_annual = target_returns * PERIODS,
  vol_quadprog         = port_vol(W_qp),
  vol_nloptr           = port_vol(W_nlopt),
  max_weight_gap       = apply(abs(W_qp - W_nlopt), 1, max)
)
cat("\n=== Solver comparison along the frontier ===\n")
print(round(solver_comparison, 5))
write.csv(solver_comparison, file.path(OUTPUT_DIR, "solver_comparison.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# 4. Test period: out-of-sample evaluation
# -----------------------------------------------------------------------------

R_oos <- R[oos_idx, , drop = FALSE]

# Rolling minimum-variance: each month t, weights estimated on all months < t
rolling_returns <- sapply(oos_idx, function(t) {
  w <- solve_qp(ledoit_wolf(R[1:(t - 1), , drop = FALSE]))
  sum(R[t, ] * w)
})

strategies_oos <- list(
  "Benchmark"                        = bench[oos_idx],
  "Equal weight (1/N)"               = rowMeans(R_oos),
  "Minimum variance (fixed weights)" = as.numeric(R_oos %*% w_gmv),
  "Minimum variance (rolling)"       = rolling_returns,
  "Index tracking"                   = as.numeric(R_oos %*% w_track)
)
table_oos <- metrics_table(strategies_oos, oos_idx)
cat(sprintf("\n=== Test period (out-of-sample, %d months) ===\n", length(oos_idx)))
print(round(table_oos, 4))
write.csv(round(table_oos, 6), file.path(OUTPUT_DIR, "metrics_out_of_sample.csv"))


# -----------------------------------------------------------------------------
# 5. Charts
# -----------------------------------------------------------------------------

theme_set(theme_minimal(base_size = 11) +
            theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank()))

# 5.1 Efficient frontier (training period)
assets_df <- data.frame(vol = sqrt(diag(Sigma_is) * PERIODS), ret = mu_is * PERIODS)
frontier_df <- rbind(
  data.frame(vol = port_vol(W_qp),    ret = target_returns * PERIODS, solver = "quadprog"),
  data.frame(vol = port_vol(W_nlopt), ret = target_returns * PERIODS, solver = "nloptr")
)
p_frontier <- ggplot() +
  geom_point(data = assets_df, aes(vol, ret), colour = "grey70", size = 1.2) +
  geom_path(data = subset(frontier_df, solver == "quadprog"), aes(vol, ret), colour = "#1f3d5c", linewidth = 1) +
  geom_point(data = subset(frontier_df, solver == "nloptr"), aes(vol, ret), colour = "#e07a5f", size = 2, shape = 1) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Efficient frontier (training period)",
       subtitle = "Line: quadprog · circles: nloptr · grey dots: individual stocks",
       x = "Annualized volatility", y = "Annualized expected return")
ggsave(file.path(OUTPUT_DIR, "figures", "efficient_frontier.png"), p_frontier, width = 8, height = 5, dpi = 150)
print(p_frontier)   # also shows the chart in RStudio's Plots tab

# 5.2 Allocation along the frontier (top 5 stocks + others)
weights_long <- function(W, solver) {
  data.frame(
    target = rep(sprintf("%.1f%%", target_returns * PERIODS * 100), times = ncol(W)),
    stock  = rep(colnames(W), each = nrow(W)),
    weight = as.vector(W),
    solver = solver
  )
}
top5 <- names(sort(colMeans(W_qp, na.rm = TRUE), decreasing = TRUE))[1:5]
alloc <- rbind(weights_long(W_qp, "quadprog"), weights_long(W_nlopt, "nloptr"))
alloc$group <- ifelse(alloc$stock %in% top5, alloc$stock, "Others")
alloc <- aggregate(weight ~ target + solver + group, data = alloc, FUN = sum)
alloc$group  <- factor(alloc$group, levels = c(top5, "Others"))
alloc$target <- factor(alloc$target, levels = unique(sprintf("%.1f%%", target_returns * PERIODS * 100)))
p_alloc <- ggplot(alloc, aes(target, weight, fill = group)) +
  geom_col(width = 0.8) +
  facet_wrap(~solver, ncol = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Allocation along the efficient frontier", x = "Annualized target return", y = "Weight", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "figures", "frontier_allocation.png"), p_alloc, width = 9, height = 6, dpi = 150)
print(p_alloc)

# 5.3 Minimum-variance portfolio weights (top 15)
gmv_df <- data.frame(stock = names(w_gmv), weight = w_gmv)
gmv_df <- head(gmv_df[order(-gmv_df$weight), ], 15)
p_gmv <- ggplot(gmv_df, aes(reorder(stock, weight), weight)) +
  geom_col(fill = "#1f3d5c") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Minimum-variance portfolio: top 15 holdings", x = NULL, y = "Weight")
ggsave(file.path(OUTPUT_DIR, "figures", "gmv_weights.png"), p_gmv, width = 7, height = 5, dpi = 150)
print(p_gmv)

# 5.4 Out-of-sample cumulative growth
growth <- do.call(rbind, lapply(names(strategies_oos), function(s) {
  data.frame(date = c(dates[IS_MONTHS], dates[oos_idx]),
             value = cumprod(c(1, 1 + strategies_oos[[s]])),
             strategy = s)
}))
p_oos <- ggplot(growth, aes(date, value, colour = strategy)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3) +
  labs(title = "Growth of $1 over the test period",
       x = NULL, y = "Portfolio value", colour = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "figures", "out_of_sample_growth.png"), p_oos, width = 9, height = 5, dpi = 150)
print(p_oos)

# 5.5 Pie chart of the minimum-variance weights
# Weights above PIE_THRESHOLD get their own slice; the rest are grouped into "Others"
# so the labels stay readable.
PIE_THRESHOLD <- 0.02
pie_w <- sort(w_gmv[w_gmv >= PIE_THRESHOLD], decreasing = TRUE)
if (1 - sum(pie_w) > 1e-6) pie_w <- c(pie_w, Others = 1 - sum(pie_w))
pie_df <- data.frame(stock = factor(names(pie_w), levels = names(pie_w)), weight = as.numeric(pie_w))
pie_df$label <- paste0(pie_df$stock, "\n", round(pie_df$weight * 100, 1), "%")
pie_cols <- setNames(hcl.colors(nrow(pie_df), "Dark 3"), pie_df$stock)
if ("Others" %in% names(pie_cols)) pie_cols["Others"] <- "grey80"
p_pie <- ggplot(pie_df, aes(x = 1, y = weight, fill = stock)) +
  geom_col(width = 1, colour = "white") +
  geom_text(aes(x = 1.75, label = label), position = position_stack(vjust = 0.5), size = 2.8) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = pie_cols) +
  xlim(0.5, 2) +
  labs(title = "Minimum-variance portfolio weights",
       subtitle = sprintf("Holdings of %.0f%% or more; the rest grouped into Others", PIE_THRESHOLD * 100)) +
  theme_void() +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")
ggsave(file.path(OUTPUT_DIR, "figures", "gmv_weights_pie.png"), p_pie, width = 7, height = 7, dpi = 150)
print(p_pie)

# 5.6 Distribution of portfolio weights along the frontier, by solver (boxplot)
WEIGHT_THRESHOLD <- 0.001   # near-zero weights are left out so the boxes stay readable
box_df <- rbind(
  data.frame(weight = as.vector(W_qp),    solver = "quadprog"),
  data.frame(weight = as.vector(W_nlopt), solver = "nloptr")
)
box_df <- subset(box_df, !is.na(weight) & weight > WEIGHT_THRESHOLD)
p_box <- ggplot(box_df, aes(solver, weight, fill = solver)) +
  geom_boxplot(outlier.size = 1.2) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Distribution of portfolio weights along the frontier",
       subtitle = sprintf("Weights above %.1f%%", WEIGHT_THRESHOLD * 100),
       x = "Solver", y = "Weight") +
  theme(legend.position = "none")
ggsave(file.path(OUTPUT_DIR, "figures", "weights_distribution_boxplot.png"), p_box, width = 7, height = 5, dpi = 150)
print(p_box)

# 5.7 Allocation of the highest-return frontier portfolio, by solver (grouped bars, top 12)
last_qp    <- max(which(rowSums(!is.na(W_qp)) > 0))
last_nlopt <- max(which(rowSums(!is.na(W_nlopt)) > 0))
solver_alloc <- rbind(
  data.frame(stock = asset_cols, weight = W_qp[last_qp, ],       solver = "quadprog"),
  data.frame(stock = asset_cols, weight = W_nlopt[last_nlopt, ], solver = "nloptr")
)
solver_alloc <- subset(solver_alloc, !is.na(weight))
avg_w  <- aggregate(weight ~ stock, data = solver_alloc, FUN = mean)
top12  <- head(avg_w$stock[order(-avg_w$weight)], 12)
solver_alloc <- subset(solver_alloc, stock %in% top12)
solver_alloc$stock <- factor(solver_alloc$stock, levels = top12)
p_solver <- ggplot(solver_alloc, aes(stock, weight, fill = solver)) +
  geom_col(position = "dodge", width = 0.75) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Highest-return frontier portfolio: allocation by solver",
       subtitle = sprintf("Target return %.1f%% per year · top 12 holdings",
                          target_returns[last_qp] * PERIODS * 100),
       x = NULL, y = "Weight", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
ggsave(file.path(OUTPUT_DIR, "figures", "solver_allocation.png"), p_solver, width = 9, height = 5, dpi = 150)
print(p_solver)

cat(sprintf("\nResults and charts saved in %s/\n", OUTPUT_DIR))
