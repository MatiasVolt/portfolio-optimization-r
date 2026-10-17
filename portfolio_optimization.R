# =============================================================================
# Optimisation de portefeuille : minimum-variance, réplication d'indice et
# frontière efficiente, avec une évaluation hors échantillon
# =============================================================================
#
# Données attendues (fichier Excel, rendements MENSUELS simples, en décimal) :
#   - Date  : date de fin de mois
#   - Bond  : taux sans risque (unité à préciser dans RF_UNIT ci-dessous)
#   - INDEX : rendement de l'indice de référence
#   - toutes les autres colonnes : rendements des actions de l'univers
#
# Démarche :
#   1. Période d'APPRENTISSAGE (IS_MONTHS premiers mois) : toutes les estimations
#      (covariance, poids, frontière, comparaison de solveurs) se font ici.
#   2. Période de TEST (mois restants) : les portefeuilles sont évalués sur des
#      données qui n'ont servi à aucune estimation.
#
# Exécuter depuis le dossier du projet : Rscript portfolio_optimization.R
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Paramètres
# -----------------------------------------------------------------------------

DATA_PATH   <- "data/returns.xlsx"   # chemin du fichier de données
DATA_SHEET  <- 1                     # nom ou numéro de la feuille
RF_UNIT     <- "monthly_decimal"     # "monthly_decimal" | "annual_decimal" | "annual_percent"
IS_MONTHS   <- 60                    # longueur de la période d'apprentissage (mois)
LOWER_BOUND <- 0.0001                # poids minimal par action (long-only)
UPPER_BOUND <- 0.20                  # poids maximal par action
N_TARGETS   <- 15                    # nombre de points sur la frontière efficiente
PERIODS     <- 12                    # périodes par an (données mensuelles)
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
# 1. Chargement et nettoyage des données
# -----------------------------------------------------------------------------

to_numeric <- function(x) as.numeric(gsub(",", ".", as.character(x)))  # gère la virgule décimale

raw <- as.data.frame(read_excel(DATA_PATH, sheet = DATA_SHEET))
stopifnot(all(c("Date", "Bond", "INDEX") %in% names(raw)))

asset_cols <- setdiff(names(raw), c("Date", "Bond", "INDEX"))
R_raw      <- suppressWarnings(sapply(raw[asset_cols], to_numeric))
# colonnes vides ou non numériques (ex. colonnes "...104" créées par Excel) : retirées
empty_cols <- colSums(!is.na(R_raw)) == 0
if (any(empty_cols)) {
  message("Colonnes ignorées (vides ou non numériques) : ", paste(asset_cols[empty_cols], collapse = ", "))
  R_raw      <- R_raw[, !empty_cols, drop = FALSE]
  asset_cols <- asset_cols[!empty_cols]
}
rf_raw     <- to_numeric(raw$Bond)
bench_raw  <- to_numeric(raw$INDEX)
dates_raw  <- as.Date(raw$Date)

# un seul filtre de lignes valides pour toutes les séries (alignement garanti)
valid <- !is.na(dates_raw) & !is.na(rf_raw) & !is.na(bench_raw) & complete.cases(R_raw)
dates <- dates_raw[valid]
R     <- R_raw[valid, , drop = FALSE]
bench <- bench_raw[valid]
rf    <- switch(RF_UNIT,
  monthly_decimal = rf_raw[valid],
  annual_decimal  = (1 + rf_raw[valid])^(1 / PERIODS) - 1,
  annual_percent  = (1 + rf_raw[valid] / 100)^(1 / PERIODS) - 1,
  stop("RF_UNIT inconnu")
)

ord   <- order(dates)                  # ordre chronologique
dates <- dates[ord]; R <- R[ord, , drop = FALSE]; bench <- bench[ord]; rf <- rf[ord]

n_obs    <- nrow(R)
n_assets <- ncol(R)

if (median(abs(R)) > 0.5) stop("Les colonnes d'actifs ressemblent à des prix, pas à des rendements.")
if (IS_MONTHS >= n_obs - 1) stop("IS_MONTHS laisse moins de 2 mois pour la période de test.")
if (LOWER_BOUND * n_assets > 1 || UPPER_BOUND * n_assets < 1) stop("Bornes de poids infaisables.")

is_idx  <- 1:IS_MONTHS
oos_idx <- (IS_MONTHS + 1):n_obs
lb <- rep(LOWER_BOUND, n_assets)
ub <- rep(UPPER_BOUND, n_assets)

cat(sprintf("%d actions, %d mois (%s -> %s)\n", n_assets, n_obs, min(dates), max(dates)))
cat(sprintf("Apprentissage : %d mois | Test : %d mois (à partir de %s)\n\n",
            length(is_idx), length(oos_idx), dates[oos_idx[1]]))


# -----------------------------------------------------------------------------
# 2. Outils : covariance, solveurs, métriques
# -----------------------------------------------------------------------------

# Covariance de Ledoit-Wolf (2004), cible = matrice identité mise à l'échelle.
# Avec plus d'actions que de mois (ex. 100 actions, 60 mois), la covariance
# empirique est singulière : l'optimiseur exploiterait alors des directions
# « sans risque » qui ne sont que du bruit d'estimation. Le shrinkage corrige cela.
ledoit_wolf <- function(X) {
  X <- scale(X, center = TRUE, scale = FALSE)
  n <- nrow(X); p <- ncol(X)
  S <- crossprod(X) / n
  target <- diag(sum(diag(S)) / p, p)
  d2 <- sum((S - target)^2)
  b2 <- sum(apply(X, 1, function(x) sum((tcrossprod(x) - S)^2))) / n^2
  delta <- min(b2, d2) / d2                      # intensité de shrinkage, entre 0 et 1
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

# min w' Sigma w  s.c.  sum(w) = 1, lb <= w <= ub  [et mu'w = target si fourni]
# Solveur 1 : programmation quadratique exacte (quadprog)
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

# Solveur 2 : optimisation non linéaire générique (nloptr, SLSQP), pour comparaison
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

# Portefeuille de réplication : min Var(R w - indice). Comme sum(w) = 1,
# R w - indice = (R - indice) w : c'est un min-variance sur les rendements actifs.
tracking_weights <- function(R_win, bench_win) {
  solve_qp(ledoit_wolf(sweep(R_win, 1, bench_win, "-")))
}

# Bornes de rendement atteignables sous les contraintes de poids (algorithme glouton)
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
    "Rendement annualisé"       = mean(r) * periods,
    "Volatilité annualisée"     = sd(r) * sqrt(periods),
    "Sharpe"                    = mean(ex) / sd(ex) * sqrt(periods),
    "Sortino"                   = if (downside > 0) mean(ex) / downside * sqrt(periods) else NA,
    "Omega (seuil 0)"           = sum(pmax(r, 0)) / sum(pmax(-r, 0)),
    "Bêta"                      = beta,
    "Treynor"                   = mean(ex) * periods / beta,
    "Alpha de Jensen (annuel)"  = (mean(ex) - beta * mean(ex_b)) * periods,
    "Max drawdown"              = max_dd,
    "Calmar"                    = mean(r) * periods / abs(max_dd),
    "VaR 95 % (mensuelle)"      = -q05,
    "ES 95 % (mensuelle)"       = -mean(r[r <= q05]),
    "Tracking error"            = te,
    "Information ratio"         = if (te > 1e-12) mean(active) * periods / te else NA,
    "Corrélation avec l'indice" = cor(r, bench)
  )
}

metrics_table <- function(strategies, idx) {
  sapply(strategies, perf_metrics, rf = rf[idx], bench = bench[idx])
}


# -----------------------------------------------------------------------------
# 3. Période d'apprentissage : portefeuilles et frontière efficiente
# -----------------------------------------------------------------------------

R_is     <- R[is_idx, , drop = FALSE]
Sigma_is <- ledoit_wolf(R_is)
mu_is    <- colMeans(R_is)
cat(sprintf("Intensité de shrinkage Ledoit-Wolf : %.2f\n", attr(Sigma_is, "shrinkage")))

w_gmv   <- solve_qp(Sigma_is)                 # minimum-variance global
w_track <- tracking_weights(R_is, bench[is_idx])
stopifnot(!anyNA(w_gmv), !anyNA(w_track))

strategies_is <- list(
  "Indice"                = bench[is_idx],
  "Équipondéré (1/N)"     = rowMeans(R_is),
  "Minimum-variance"      = as.numeric(R_is %*% w_gmv),
  "Réplication d'indice"  = as.numeric(R_is %*% w_track)
)
table_is <- metrics_table(strategies_is, is_idx)
cat("\n=== Apprentissage (in-sample : performances optimistes par construction) ===\n")
print(round(table_is, 4))
write.csv(round(table_is, 6), file.path(OUTPUT_DIR, "metrics_in_sample.csv"))

# --- Frontière efficiente : du portefeuille min-variance au rendement max atteignable
target_returns <- seq(sum(w_gmv * mu_is), feasible_return(mu_is, maximize = TRUE),
                      length.out = N_TARGETS)
# le dernier point (rendement max) n'a qu'une solution, souvent mal résolue numériquement
target_returns <- head(target_returns, -1)

W_qp    <- t(sapply(target_returns, function(m) solve_qp(Sigma_is, mu_is, m)))
W_nlopt <- t(sapply(target_returns, function(m) solve_nlopt(Sigma_is, mu_is, m)))
colnames(W_qp) <- colnames(W_nlopt) <- asset_cols

port_vol <- function(W) apply(W, 1, function(w) sqrt(as.numeric(t(w) %*% Sigma_is %*% w) * PERIODS))

solver_comparison <- data.frame(
  rendement_cible_annuel = target_returns * PERIODS,
  vol_quadprog           = port_vol(W_qp),
  vol_nloptr             = port_vol(W_nlopt),
  ecart_max_poids        = apply(abs(W_qp - W_nlopt), 1, max)
)
cat("\n=== Comparaison des solveurs le long de la frontière ===\n")
print(round(solver_comparison, 5))
write.csv(solver_comparison, file.path(OUTPUT_DIR, "solver_comparison.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# 4. Période de test : évaluation hors échantillon
# -----------------------------------------------------------------------------

R_oos <- R[oos_idx, , drop = FALSE]

# Min-variance glissant : chaque mois t, poids estimés sur tous les mois < t
rolling_returns <- sapply(oos_idx, function(t) {
  w <- solve_qp(ledoit_wolf(R[1:(t - 1), , drop = FALSE]))
  sum(R[t, ] * w)
})

strategies_oos <- list(
  "Indice"                         = bench[oos_idx],
  "Équipondéré (1/N)"              = rowMeans(R_oos),
  "Minimum-variance (poids fixes)" = as.numeric(R_oos %*% w_gmv),
  "Minimum-variance (glissant)"    = rolling_returns,
  "Réplication d'indice"           = as.numeric(R_oos %*% w_track)
)
table_oos <- metrics_table(strategies_oos, oos_idx)
cat(sprintf("\n=== Test (out-of-sample, %d mois) ===\n", length(oos_idx)))
print(round(table_oos, 4))
write.csv(round(table_oos, 6), file.path(OUTPUT_DIR, "metrics_out_of_sample.csv"))


# -----------------------------------------------------------------------------
# 5. Graphiques
# -----------------------------------------------------------------------------

theme_set(theme_minimal(base_size = 11) +
            theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank()))

# 5.1 Frontière efficiente (apprentissage)
assets_df <- data.frame(vol = sqrt(diag(Sigma_is) * PERIODS), ret = mu_is * PERIODS)
frontier_df <- rbind(
  data.frame(vol = port_vol(W_qp),    ret = target_returns * PERIODS, solveur = "quadprog"),
  data.frame(vol = port_vol(W_nlopt), ret = target_returns * PERIODS, solveur = "nloptr")
)
p_frontier <- ggplot() +
  geom_point(data = assets_df, aes(vol, ret), colour = "grey70", size = 1.2) +
  geom_path(data = subset(frontier_df, solveur == "quadprog"), aes(vol, ret), colour = "#1f3d5c", linewidth = 1) +
  geom_point(data = subset(frontier_df, solveur == "nloptr"), aes(vol, ret), colour = "#e07a5f", size = 2, shape = 1) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Frontière efficiente (période d'apprentissage)",
       subtitle = "Ligne : quadprog · cercles : nloptr · points gris : actions individuelles",
       x = "Volatilité annualisée", y = "Rendement espéré annualisé")
ggsave(file.path(OUTPUT_DIR, "figures", "efficient_frontier.png"), p_frontier, width = 8, height = 5, dpi = 150)

# 5.2 Composition le long de la frontière (5 actions principales + autres)
weights_long <- function(W, solver) {
  data.frame(
    cible  = rep(sprintf("%.1f %%", target_returns * PERIODS * 100), times = ncol(W)),
    action = rep(colnames(W), each = nrow(W)),
    poids  = as.vector(W),
    solveur = solver
  )
}
top5 <- names(sort(colMeans(W_qp, na.rm = TRUE), decreasing = TRUE))[1:5]
alloc <- rbind(weights_long(W_qp, "quadprog"), weights_long(W_nlopt, "nloptr"))
alloc$groupe <- ifelse(alloc$action %in% top5, alloc$action, "Autres")
alloc <- aggregate(poids ~ cible + solveur + groupe, data = alloc, FUN = sum)
alloc$groupe <- factor(alloc$groupe, levels = c(top5, "Autres"))
alloc$cible  <- factor(alloc$cible, levels = unique(sprintf("%.1f %%", target_returns * PERIODS * 100)))
p_alloc <- ggplot(alloc, aes(cible, poids, fill = groupe)) +
  geom_col(width = 0.8) +
  facet_wrap(~solveur, ncol = 1) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Allocation le long de la frontière", x = "Rendement cible annualisé", y = "Poids", fill = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "figures", "frontier_allocation.png"), p_alloc, width = 9, height = 6, dpi = 150)

# 5.3 Poids du portefeuille minimum-variance (15 principaux)
gmv_df <- data.frame(action = names(w_gmv), poids = w_gmv)
gmv_df <- head(gmv_df[order(-gmv_df$poids), ], 15)
p_gmv <- ggplot(gmv_df, aes(reorder(action, poids), poids)) +
  geom_col(fill = "#1f3d5c") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Portefeuille minimum-variance : 15 premières lignes", x = NULL, y = "Poids")
ggsave(file.path(OUTPUT_DIR, "figures", "gmv_weights.png"), p_gmv, width = 7, height = 5, dpi = 150)

# 5.4 Croissance cumulée hors échantillon
growth <- do.call(rbind, lapply(names(strategies_oos), function(s) {
  data.frame(date = c(dates[IS_MONTHS], dates[oos_idx]),
             valeur = cumprod(c(1, 1 + strategies_oos[[s]])),
             strategie = s)
}))
p_oos <- ggplot(growth, aes(date, valeur, colour = strategie)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = 0.3) +
  labs(title = "Valeur de 1 € investi sur la période de test",
       x = NULL, y = "Valeur du portefeuille", colour = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "figures", "out_of_sample_growth.png"), p_oos, width = 9, height = 5, dpi = 150)

cat(sprintf("\nRésultats et figures enregistrés dans %s/\n", OUTPUT_DIR))
