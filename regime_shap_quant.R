# =============================================================================
#  RegimeShap-Quant
#  A Hybrid Rule-Based + Gradient-Boosting Regime Framework for
#  Volatility-Targeted Allocation with Shapley Explainability
# -----------------------------------------------------------------------------
#  Pipeline
#    1. Market data   
#    2. Features      
#    3. Rule layer    
#    4. ML layer      
#    5. Hybrid arbiter
#    6. Allocator      
#    7. Back-test      
#    8. Statistics  
#    9. Explainability
#   10. Robustness
# =============================================================================

suppressPackageStartupMessages({
  library(gbm); library(ggplot2); library(zoo); library(scales); library(parallel)
})

# ------------------------------ 0. CONFIGURATION -----------------------------
CFG <- list(
  data_source  = "simulated",            # "simulated"  or  "yahoo"
  tickers      = c(EQ = "SPY", BD = "TLT", GD = "GLD", CM = "DBC"),
  start        = "2006-01-01", end = "2025-12-31",   # used for yahoo only
  seed         = 2026,
  n_days       = 5040,                   # ~20 years (simulated)
  target_vol   = 0.10,                   # annualised volatility target
  max_weight   = 1.00,                   # no leverage
  tc_bps       = 5,                      # one-way cost, basis points
  band         = 0.05,                   # no-trade band on equity weight
  train_min    = 1260,                   # first training window (~5y)
  retrain_every= 126,                    # re-fit every ~6 months
  horizon      = 21,                     # forecast horizon (trading days)
  gbm_trees    = 600, gbm_depth = 2, gbm_shrink = 0.02,   # max trees; early-stopped per window
  shap_explain = 150, shap_perm = 50, shap_bg = 100,
  n_boot       = 1000, block = 10,
  n_mc         = 20,                     # independent simulated markets
  cores        = max(1, min(8, detectCores() - 1)),
  quick        = identical(Sys.getenv("RSQ_QUICK"), "1"),   # smoke-test switch
  out_dir      = "outputs"
)
if (CFG$quick) { CFG$n_mc <- 0; CFG$shap_explain <- 20; CFG$shap_perm <- 10; CFG$n_boot <- 100 }
ANN <- sqrt(252)
dir.create(file.path(CFG$out_dir, "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CFG$out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
REG <- c("Calm", "Stress", "Crisis")
save_tab <- function(x, name) write.csv(x, file.path(CFG$out_dir, "tables", name), row.names = FALSE)
save_fig <- function(p, name, w = 3.5, h = 2.6)
  ggsave(file.path(CFG$out_dir, "figures", name), p, width = w, height = h, dpi = 300, bg = "white")
theme_ieee <- function() theme_bw(base_size = 7, base_family = "serif") +
  theme(panel.grid.minor = element_blank(), legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 6), legend.text = element_text(size = 6),
        plot.margin = margin(3, 4, 3, 3))

# ------------------------------ 1. DATA --------------------------------------
# Ground-truth 3-state Markov chain with fat tails and stochastic volatility.
simulate_market <- function(n, seed) {
  set.seed(seed)
  P <- matrix(c(.985, .013, .002,
                .035, .960, .005,
                .010, .060, .930), 3, byrow = TRUE)
  mu <- rbind(c(.12, .03, .04, .05), c(0, .04, .07, -.02), c(-.45, .10, .10, -.30))
  sg <- rbind(c(.11, .06, .13, .14), c(.22, .09, .16, .20), c(.50, .12, .22, .35))
  mk <- function(a, b, c, d, e, f) { m <- diag(4); m[upper.tri(m)] <- c(a, b, c, d, e, f)
                                      m[lower.tri(m)] <- t(m)[lower.tri(m)]; m }
  # upper-tri order: (1,2) (1,3) (2,3) (1,4) (2,4) (3,4)
  C <- list(mk(-.20, .05, .20, .25, -.05, .30),
            mk(-.35, .10, .25, .45, -.10, .35),
            mk(-.45, .15, .30, .65, -.15, .45))
  state <- integer(n); state[1] <- 1L
  for (t in 2:n) state[t] <- sample.int(3, 1, prob = P[state[t - 1], ])
  h <- numeric(n); for (t in 2:n) h[t] <- .97 * h[t - 1] + .06 * rnorm(1)   # vol clustering
  dfree <- 5
  Z <- matrix(rt(n * 4, dfree) * sqrt((dfree - 2) / dfree), n, 4)           # fat tails
  for (k in 1:3) { i <- which(state == k); if (length(i)) Z[i, ] <- Z[i, , drop = FALSE] %*% chol(C[[k]]) }
  R <- mu[state, ] / 252 + sg[state, ] * exp(h) / ANN * Z
  R <- pmax(pmin(R, 0.25), -0.25)
  d <- seq(as.Date("2005-01-03"), by = "day", length.out = ceiling(n * 1.5))
  d <- d[!weekdays(d) %in% c("Saturday", "Sunday")][1:n]
  out <- data.frame(date = d, R); names(out)[2:5] <- c("EQ", "BD", "GD", "CM")
  out$regime <- state
  attr(out, "sg_eq") <- sg[, 1]
  out
}

load_yahoo <- function() {
  if (!requireNamespace("quantmod", quietly = TRUE)) stop("Install 'quantmod' for real data.")
  px <- lapply(CFG$tickers, function(tk)
    quantmod::Ad(quantmod::getSymbols(tk, from = CFG$start, to = CFG$end, auto.assign = FALSE)))
  px <- na.omit(do.call(merge, px)); colnames(px) <- names(CFG$tickers)
  p <- zoo::coredata(px); r <- p[-1, ] / p[-nrow(p), ] - 1
  out <- data.frame(date = zoo::index(px)[-1], r); out$regime <- NA_integer_
  attr(out, "sg_eq") <- c(.11, .22, .50)
  out
}

get_data <- function(seed = CFG$seed)
  if (CFG$data_source == "yahoo") load_yahoo() else simulate_market(CFG$n_days, seed)

# ------------------------------ 2. FEATURES ----------------------------------
FEATS <- c("ret_1", "ret_5", "ret_21", "mom_126", "rv_10", "rv_21", "rv_63",
           "vol_ratio", "vol_of_vol", "down_vol", "skew_63", "dd_252",
           "corr_eq_bd", "bd_mom_63", "gd_mom_63")

build_features <- function(d, horizon = CFG$horizon) {
  r <- d$EQ; n <- length(r)
  rsd  <- function(x, k) rollapplyr(x, k, sd, fill = NA)
  rsum <- function(x, k) rollsumr(x, k, fill = NA)
  f <- data.frame(date = d$date)
  f$ret_1 <- r; f$ret_5 <- rsum(r, 5); f$ret_21 <- rsum(r, 21); f$mom_126 <- rsum(r, 126)
  f$rv_10 <- rsd(r, 10) * ANN; f$rv_21 <- rsd(r, 21) * ANN; f$rv_63 <- rsd(r, 63) * ANN
  f$vol_ratio  <- f$rv_10 / f$rv_63
  f$vol_of_vol <- rsd(f$rv_21, 63)
  f$down_vol   <- sqrt(rollmeanr(pmin(r, 0)^2, 21, fill = NA)) * ANN
  f$skew_63    <- rollapplyr(r, 63, function(v) mean(((v - mean(v)) / sd(v))^3), fill = NA)
  px <- cumprod(1 + r)
  f$dd_252     <- px / rollmaxr(px, 252, fill = NA) - 1      # drawdown vs 1-year high
  f$corr_eq_bd <- rollapplyr(cbind(r, d$BD), 63, function(m) cor(m[, 1], m[, 2]),
                             by.column = FALSE, fill = NA)
  f$bd_mom_63  <- rsum(d$BD, 63); f$gd_mom_63 <- rsum(d$GD, 63)
  # target: log of next-21-day realised volatility (annualised)
  m <- rollmeanr(r^2, horizon, fill = NA)
  f$fwd_vol <- c(sqrt(m[(horizon + 1):n]) * ANN, rep(NA, horizon))
  f$y <- log(f$fwd_vol)
  f
}

# ------------------------------ 3. RULE LAYER --------------------------------
# Domain-knowledge "signature" rules with Schmitt-trigger hysteresis for Crisis.
# Thresholds are calibrated ONCE on independent development paths (seeds 101-103)
# by maximising macro-F1 (see calibrate_rules) and then frozen -- the evaluation
# paths never influence them.
RULES <- list(c_rv10 = 0.40, c_rv21 = 0.35, c_dd = -0.25, c_exit = 0.26,
              s_rv21 = 0.15, s_dd = -0.12)   # output of calibrate_rules() on seeds 101-103

rule_regime <- function(f, R = RULES) {
  n <- nrow(f); st <- rep(1L, n); in_crisis <- FALSE
  for (t in seq_len(n)) {
    if (anyNA(c(f$rv_10[t], f$rv_21[t], f$dd_252[t]))) { st[t] <- 1L; next }
    enter <- (f$rv_10[t] > R$c_rv10 & f$rv_21[t] > R$c_rv21) | (f$dd_252[t] < R$c_dd & f$rv_21[t] > R$c_rv21)
    if (!in_crisis && enter) in_crisis <- TRUE
    else if (in_crisis && f$rv_21[t] < R$c_exit) in_crisis <- FALSE
    st[t] <- if (in_crisis) 3L else if (f$rv_21[t] > R$s_rv21 || f$dd_252[t] < R$s_dd) 2L else 1L
  }
  st
}

calibrate_rules <- function(seeds = 101:103) {
  grid <- expand.grid(c_rv10 = c(.30, .35, .40, .45), c_rv21 = c(.25, .30, .35), c_dd = c(-.15, -.20, -.25),
                      c_exit = c(.18, .22, .26), s_rv21 = c(.15, .17, .19, .21), s_dd = c(-.08, -.10, -.12))
  dev <- lapply(seeds, function(sd) { d <- simulate_market(CFG$n_days, sd); f <- build_features(d)
                                       k <- complete.cases(f[, FEATS]); list(f = f[k, ], tru = d$regime[k]) })
  grid$macroF1 <- sapply(seq_len(nrow(grid)), function(i) {
    R <- as.list(grid[i, 1:6])
    mean(sapply(dev, function(z) { pr <- rule_regime(z$f, R)
      cm <- table(factor(z$tru, 1:3), factor(pr, 1:3)); p <- diag(cm) / colSums(cm); r <- diag(cm) / rowSums(cm)
      mean(2 * p * r / (p + r), na.rm = TRUE) })) })
  grid[order(-grid$macroF1), ][1:5, ]
}

# ------------------------------ 4. ML LAYER (walk-forward) -------------------
walk_forward <- function(f, cfg = CFG, keep_last = FALSE) {
  n <- nrow(f); ok <- which(complete.cases(f[, FEATS]))
  s0 <- ok[1] + cfg$train_min - 1
  refits <- seq(s0, n - 1, by = cfg$retrain_every)
  p_ml <- p_har <- rep(NA_real_, n); last <- NULL
  for (s in refits) {
    tr <- ok[ok <= s - cfg$horizon]; tr <- tr[!is.na(f$y[tr])]   # purge: labels known at s
    tr_df <- f[tr, c("y", FEATS)]
    fit <- gbm(y ~ ., data = tr_df, distribution = "gaussian", n.trees = cfg$gbm_trees,
               interaction.depth = cfg$gbm_depth, shrinkage = cfg$gbm_shrink, bag.fraction = 0.5,
               train.fraction = 0.8, n.minobsinnode = 30, verbose = FALSE)  # chronological 80/20 split
    best <- max(50L, gbm.perf(fit, method = "test", plot.it = FALSE))        # early stopping
    te <- s:min(s + cfg$retrain_every - 1, n - 1)
    p_ml[te] <- predict(fit, f[te, FEATS], n.trees = best)
    har_df <- data.frame(y = f$y[tr], l10 = log(f$rv_10[tr]), l21 = log(f$rv_21[tr]), l63 = log(f$rv_63[tr]))
    har <- lm(y ~ l10 + l21 + l63, data = har_df)                # HAR-type benchmark
    p_har[te] <- predict(har, data.frame(l10 = log(f$rv_10[te]), l21 = log(f$rv_21[te]), l63 = log(f$rv_63[te])))
    if (keep_last) last <- list(fit = fit, train = tr_df, s = s, best = best)
  }
  list(p_ml = p_ml, p_har = p_har, oos = s0:(n - 1), last = last)
}

# ------------------------------ 5-7. ALLOCATION & BACK-TEST ------------------
clip <- function(x, lo = 0, hi = CFG$max_weight) pmin(pmax(x, lo), hi)
apply_band <- function(w, band = CFG$band) {
  out <- w; for (i in 2:length(w)) out[i] <- if (abs(w[i] - out[i - 1]) > band) w[i] else out[i - 1]; out
}

build_weights <- function(f, d, wf, rule, tv = CFG$target_vol) {
  i <- wf$oos
  W <- list(
    "Equity B&H"     = rep(1, length(i)),
    "60/40"          = rep(0.6, length(i)),
    "Naive Vol-Tgt"  = clip(tv / f$rv_21[i]),
    "HAR Vol-Tgt"    = clip(tv / exp(wf$p_har[i])),
    "Rule-Only"      = c(1, .5, 0)[rule[i]],
    "ML-Only"        = clip(tv / exp(wf$p_ml[i])),
    "Hybrid-Override"= ifelse(rule[i] == 3L, 0, clip(tv / exp(wf$p_ml[i]))),
    # Conservative arbiter: ML weight capped by a regime ceiling (1 / 0.75 / 0.25)
    "Hybrid-Cap"     = pmin(clip(tv / exp(wf$p_ml[i])), c(1, .75, .25)[rule[i]])
  )
  if (!all(is.na(d$regime))) W[["Oracle (true regime)"]] <- clip(tv / attr(d, "sg_eq")[d$regime[i]])
  for (k in setdiff(names(W), c("Equity B&H", "60/40"))) W[[k]] <- apply_band(W[[k]])
  W
}

run_backtest <- function(d, wf, W, tc_bps = CFG$tc_bps) {
  i <- wf$oos; r_eq <- d$EQ[i + 1]; r_def <- 0.5 * d$BD[i + 1] + 0.5 * d$GD[i + 1]
  R <- lapply(names(W), function(k) {
    w <- W[[k]]; dw <- abs(w - c(w[1], head(w, -1)))
    cost <- if (k %in% c("Equity B&H", "60/40")) 0 else 2 * tc_bps / 1e4 * dw
    w * r_eq + (1 - w) * r_def - cost
  }); names(R) <- names(W)
  list(ret = as.data.frame(R, check.names = FALSE), dates = d$date[i + 1])
}

# ------------------------------ 8. METRICS & STATISTICS ----------------------
sharpe <- function(x) mean(x) / sd(x) * ANN
metrics <- function(r, w, tv = CFG$target_vol) {
  n <- length(r); eq <- cumprod(1 + r); dd <- eq / cummax(eq) - 1
  q <- quantile(r, 0.05)
  data.frame(AnnReturn = prod(1 + r)^(252 / n) - 1, AnnVol = sd(r) * ANN, Sharpe = sharpe(r),
             Sortino = mean(r) / sqrt(mean(pmin(r, 0)^2)) * ANN, MaxDD = min(dd),
             Calmar = (prod(1 + r)^(252 / n) - 1) / abs(min(dd)), CVaR5 = mean(r[r <= q]),
             Turnover = sum(abs(diff(w))) * 252 / n)
}
perf_table <- function(bt, W) {
  do.call(rbind, lapply(names(W), function(k) cbind(Strategy = k, metrics(bt$ret[[k]], W[[k]]))))
}
boot_sharpe_diff <- function(a, b, B = CFG$n_boot, blk = CFG$block) {
  n <- length(a); nb <- ceiling(n / blk); d_hat <- sharpe(a) - sharpe(b)
  bd <- replicate(B, {
    st <- sample.int(n, nb, TRUE); ix <- (as.vector(outer(0:(blk - 1), st, "+")) - 1) %% n + 1
    sharpe(a[ix[1:n]]) - sharpe(b[ix[1:n]])
  })
  c(diff = d_hat, lo = unname(quantile(bd, .025)), hi = unname(quantile(bd, .975)),
    p = mean(abs(bd - d_hat) >= abs(d_hat)))
}
qlike <- function(rv_real, vol_fc) mean(log(vol_fc^2) + rv_real^2 / vol_fc^2)

# ------------------------------ 9. SHAPLEY EXPLAINABILITY --------------------
# Monte-Carlo permutation estimator of interventional Shapley values
# (Strumbelj & Kononenko 2014); same quantity SHAP approximates (Lundberg & Lee 2017).
shapley_mc <- function(fit, X, Xbg, n_perm, ntrees) {
  p <- ncol(X); phi <- matrix(0, nrow(X), p, dimnames = list(NULL, colnames(X)))
  base <- numeric(nrow(X)); fx <- predict(fit, X, n.trees = ntrees)
  for (i in seq_len(nrow(X))) {
    x <- as.numeric(X[i, ]); M <- matrix(0, n_perm * (p + 1), p); ord <- matrix(0L, n_perm, p)
    for (m in seq_len(n_perm)) {
      o <- sample.int(p); z <- as.numeric(Xbg[sample.int(nrow(Xbg), 1), ]); ord[m, ] <- o
      rows <- (m - 1) * (p + 1) + 1:(p + 1); M[rows[1], ] <- z
      for (j in 1:p) { z[o[j]] <- x[o[j]]; M[rows[j + 1], ] <- z }
    }
    colnames(M) <- colnames(X)
    pr <- matrix(predict(fit, as.data.frame(M), n.trees = ntrees), p + 1)
    base[i] <- mean(pr[1, ])
    for (m in seq_len(n_perm)) for (j in 1:p) phi[i, ord[m, j]] <- phi[i, ord[m, j]] + pr[j + 1, m] - pr[j, m]
  }
  list(phi = phi / n_perm, base = base, fx = fx)
}

# ------------------------------ PLOT HELPERS ---------------------------------
plot_architecture <- function() {
  b <- data.frame(
    x = c(.5, .5, .24, .76, .5, .5, .5, .5), y = c(.95, .81, .64, .64, .47, .33, .19, .05),
    w = c(.78, .78, .42, .42, .78, .78, .78, .78), h = .085,
    lab = c("Market data (returns of EQ, BD, GD, CM)", "Leakage-free feature engineering (15 features)",
            "RULE LAYER\nvol / drawdown signatures", "ML LAYER\ngradient-boosting vol forecast",
            "Hybrid arbiter: crisis flag overrides, else ML", "Vol-targeted allocator + no-trade band",
            "Purged walk-forward back-test + bootstrap tests", "Shapley explainer + analytics tables"),
    fill = c("#d9d9d9", "#d9d9d9", "#fdd9b5", "#b7d7f0", "#d5e8c8", "#d5e8c8", "#d9d9d9", "#f0d9f0"))
  seg <- data.frame(x = c(.5, .5, .5, .24, .76, .5, .5, .5), y = c(.907, .767, .767, .597, .597, .427, .287, .147),
                    xe = c(.5, .24, .76, .5, .5, .5, .5, .5), ye = c(.853, .683, .683, .513, .513, .373, .233, .093))
  ggplot() +
    geom_rect(data = b, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2 - .01, ymax = y + h / 2 + .01, fill = fill),
              colour = "black", linewidth = .25) +
    geom_text(data = b, aes(x, y, label = lab), size = 2.1, family = "serif", lineheight = .85) +
    geom_segment(data = seg, aes(x = x, y = y, xend = xe, yend = ye), linewidth = .25,
                 arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
    scale_fill_identity() + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_void()
}

plot_waterfall <- function(phi_row, base, fx, title, top = 7) {
  o <- order(abs(phi_row), decreasing = TRUE); keep <- o[1:top]
  df <- data.frame(feat = c(names(phi_row)[keep], "other"), v = c(phi_row[keep], sum(phi_row[-keep])))
  df <- df[order(abs(df$v)), ]; df$feat <- factor(df$feat, levels = df$feat)
  ggplot(df, aes(v, feat, fill = v > 0)) + geom_col(width = .7) +
    scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"), guide = "none") +
    geom_vline(xintercept = 0, linewidth = .3) +
    labs(x = "Shapley contribution to log-vol forecast", y = NULL,
         title = sprintf("%s: base %.2f -> forecast %.2f (vol %.1f%%)", title, base, fx, 100 * exp(fx))) +
    theme_ieee() + theme(plot.title = element_text(size = 6))
}

# ============================== MAIN RUN =====================================
cat("== RegimeShap-Quant: main run (seed", CFG$seed, ") ==\n")
t0 <- Sys.time()
d  <- get_data(); f <- build_features(d); rule <- rule_regime(f)
wf <- walk_forward(f, keep_last = TRUE)
W  <- build_weights(f, d, wf, rule); bt <- run_backtest(d, wf, W)
tab_perf <- perf_table(bt, W)
cat("Back-test finished in", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

oos <- wf$oos
# ---- Statistical tests vs. benchmarks
bench <- c("Equity B&H", "60/40", "Naive Vol-Tgt", "HAR Vol-Tgt", "ML-Only", "Rule-Only", "Hybrid-Override")
set.seed(CFG$seed)
tab_boot <- do.call(rbind, lapply(bench, function(b) {
  z <- boot_sharpe_diff(bt$ret[["Hybrid-Cap"]], bt$ret[[b]])
  data.frame(Benchmark = b, SharpeDiff = z["diff"], CI_lo = z["lo"], CI_hi = z["hi"], p_value = z["p"], row.names = NULL)
}))

# ---- Forecast accuracy (out-of-sample)
y <- f$y[oos]; ok <- !is.na(y); rvr <- f$fwd_vol[oos][ok]
fc <- list("GBM (ML layer)" = exp(wf$p_ml[oos])[ok], "HAR-OLS" = exp(wf$p_har[oos])[ok],
           "Trailing RV21" = f$rv_21[oos][ok])
tab_fc <- do.call(rbind, lapply(names(fc), function(k) data.frame(
  Forecaster = k, RMSE_logvol = sqrt(mean((log(fc[[k]]) - y[ok])^2)),
  MAE_logvol = mean(abs(log(fc[[k]]) - y[ok])), QLIKE = qlike(rvr, fc[[k]]),
  Corr = cor(log(fc[[k]]), y[ok]))))

# ---- Regime diagnostics (only meaningful with ground truth)
tab_reg <- NULL; tab_conf <- NULL
if (!all(is.na(d$regime))) {
  tru <- d$regime[oos]; pr <- rule[oos]
  cm <- table(True = factor(REG[tru], REG), Rule = factor(REG[pr], REG))
  tab_conf <- as.data.frame.matrix(cm); tab_conf <- cbind(True = rownames(tab_conf), tab_conf)
  prec <- diag(cm) / colSums(cm); rec <- diag(cm) / rowSums(cm); f1 <- 2 * prec * rec / (prec + rec)
  ep <- rle(tru == 3L); ends <- cumsum(ep$lengths); starts <- ends - ep$lengths + 1
  ev <- which(ep$values)
  delay <- sapply(ev, function(e) { s <- starts[e]; z <- which(pr[s:ends[e]] == 3L)[1]; if (is.na(z)) NA else z - 1 })
  tab_reg <- data.frame(Accuracy = sum(diag(cm)) / sum(cm), MacroF1 = mean(f1, na.rm = TRUE),
                        CrisisRecall = unname(rec[3]), CrisisPrecision = unname(prec[3]),
                        CrisisEpisodes = length(ev), Detected = sum(!is.na(delay)),
                        MeanDelayDays = mean(delay, na.rm = TRUE))
  save_tab(tab_conf, "rule_confusion.csv"); save_tab(tab_reg, "regime_diagnostics.csv")
}
# crisis-period protection: strategy returns on days that are truly Crisis
if (!all(is.na(d$regime))) {
  cr <- d$regime[oos + 1] == 3L
  tab_crisis <- data.frame(Strategy = names(W),
                           MeanDailyRet_Crisis_bps = sapply(names(W), function(k) 1e4 * mean(bt$ret[[k]][cr])),
                           AvgEquityWeight_Crisis = sapply(names(W), function(k) mean(W[[k]][cr])))
  save_tab(tab_crisis, "crisis_period_protection.csv")
}

# ---- Sensitivity (forecasts fixed, only allocation / cost vary)
sens <- do.call(rbind, lapply(c(.06, .08, .10, .12, .15), function(tv) {
  Wt <- build_weights(f, d, wf, rule, tv); bt2 <- run_backtest(d, wf, Wt)
  data.frame(TargetVol = tv, Sharpe = sharpe(bt2$ret[["Hybrid-Cap"]]), MaxDD = metrics(bt2$ret[["Hybrid-Cap"]], Wt[["Hybrid-Cap"]])$MaxDD,
             RealisedVol = sd(bt2$ret[["Hybrid-Cap"]]) * ANN)
}))
sens_tc <- do.call(rbind, lapply(c(0, 5, 10, 25, 50), function(tc) {
  bt2 <- run_backtest(d, wf, W, tc); data.frame(Cost_bps = tc, Sharpe = sharpe(bt2$ret[["Hybrid-Cap"]]),
                                               AnnReturn = metrics(bt2$ret[["Hybrid-Cap"]], W[["Hybrid-Cap"]])$AnnReturn)
}))

# ---- Shapley analysis on the last walk-forward model
cat("Computing Shapley values ...\n")
set.seed(CFG$seed)
lf <- wf$last; Xall <- f[oos, FEATS]
reg_all <- if (!all(is.na(d$regime))) d$regime[oos] else rule[oos]
cand <- which(complete.cases(Xall))                      # stratified by regime so Crisis is represented
pick <- sort(unlist(lapply(1:3, function(k) { ii <- cand[reg_all[cand] == k]
  if (length(ii)) sample(ii, min(round(CFG$shap_explain / 3), length(ii))) })))
Xe <- Xall[pick, ]; Xbg <- lf$train[sample(nrow(lf$train), CFG$shap_bg), FEATS]
sh <- shapley_mc(lf$fit, Xe, Xbg, CFG$shap_perm, lf$best)
eff_err <- abs(rowSums(sh$phi) - (sh$fx - sh$base))
imp <- sort(colMeans(abs(sh$phi)), decreasing = TRUE)
tab_shap <- data.frame(Feature = names(imp), MeanAbsShapley = as.numeric(imp))
reg_pick <- if (!all(is.na(d$regime))) d$regime[oos][pick] else rule[oos][pick]
tab_shap_reg <- do.call(rbind, lapply(1:3, function(k) {
  if (sum(reg_pick == k) < 2) return(NULL)
  data.frame(Regime = REG[k], Feature = colnames(sh$phi), MeanAbsShapley = colMeans(abs(sh$phi[reg_pick == k, , drop = FALSE])),
             N = sum(reg_pick == k), row.names = NULL)
}))
tab_shap_meta <- data.frame(N_explained = nrow(Xe), Perms = CFG$shap_perm, Background = CFG$shap_bg,
                            Mean_efficiency_error = mean(eff_err), Max_efficiency_error = max(eff_err))

# ---- Save tables
save_tab(tab_perf, "performance.csv"); save_tab(tab_boot, "bootstrap_tests.csv")
save_tab(tab_fc, "forecast_accuracy.csv"); save_tab(sens, "sensitivity_target_vol.csv")
save_tab(sens_tc, "sensitivity_costs.csv"); save_tab(tab_shap, "shap_global.csv")
save_tab(tab_shap_reg, "shap_by_regime.csv"); save_tab(tab_shap_meta, "shap_meta.csv")

# ---- Figures
save_fig(plot_architecture(), "fig1_architecture.png", 3.4, 4.2)

dd <- data.frame(date = d$date, px = cumprod(1 + d$EQ), reg = if (all(is.na(d$regime))) NA else REG[d$regime])
if (!all(is.na(d$regime))) {
  rr <- rle(d$regime); e <- cumsum(rr$lengths); s <- e - rr$lengths + 1
  rect <- data.frame(xmin = d$date[s], xmax = d$date[e], reg = factor(REG[rr$values], REG))
  p2 <- ggplot() +
    geom_rect(data = rect, aes(xmin = xmin, xmax = xmax, ymin = 1e-2, ymax = Inf, fill = reg), alpha = .35) +
    geom_line(data = dd, aes(date, px), linewidth = .3) + scale_y_log10(labels = label_number(accuracy = 1)) +
    scale_fill_manual(values = c(Calm = "#ffffff00", Stress = "#fdb863", Crisis = "#d7301f"), name = "True regime") +
    coord_cartesian(ylim = range(dd$px) * c(.9, 1.1)) +
    labs(x = NULL, y = "Equity index (log)") + theme_ieee() + theme(legend.position = "bottom")
  save_fig(p2, "fig2_simulated_market.png", 3.5, 2.3)
}

eqc <- data.frame(date = rep(bt$dates, ncol(bt$ret)), strat = rep(names(bt$ret), each = nrow(bt$ret)),
                  eq = as.vector(apply(bt$ret, 2, function(x) cumprod(1 + x))))
eqc$dd <- ave(eqc$eq, eqc$strat, FUN = function(x) x / cummax(x) - 1)
show <- c("Equity B&H", "Naive Vol-Tgt", "ML-Only", "Hybrid-Override", "Hybrid-Cap")
pal  <- c("Equity B&H" = "grey40", "Naive Vol-Tgt" = "#4daf4a", "ML-Only" = "#377eb8", "Hybrid-Override" = "#ff7f00", "Hybrid-Cap" = "#e41a1c")
p3a <- ggplot(subset(eqc, strat %in% show), aes(date, eq, colour = strat)) + geom_line(linewidth = .35) +
  scale_y_log10() + scale_colour_manual(values = pal, name = NULL) + labs(x = NULL, y = "Growth of 1 (log)") +
  theme_ieee() + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
p3b <- ggplot(subset(eqc, strat %in% show), aes(date, dd, colour = strat)) + geom_line(linewidth = .3) +
  scale_colour_manual(values = pal, guide = "none") + scale_y_continuous(labels = percent) +
  labs(x = NULL, y = "Drawdown") + theme_ieee()
save_fig(p3a, "fig3_equity_curves.png", 3.5, 2.5); save_fig(p3b, "fig4_drawdowns.png", 3.5, 1.9)

fd <- data.frame(date = d$date[oos], real = f$fwd_vol[oos], ml = exp(wf$p_ml[oos]), har = exp(wf$p_har[oos]))
fd <- na.omit(fd)
p5 <- ggplot(fd, aes(date)) + geom_line(aes(y = real, colour = "Realised (21d fwd)"), linewidth = .3) +
  geom_line(aes(y = ml, colour = "GBM forecast"), linewidth = .3) +
  scale_colour_manual(values = c("Realised (21d fwd)" = "grey50", "GBM forecast" = "#e41a1c"), name = NULL) +
  scale_y_continuous(labels = percent) + labs(x = NULL, y = "Annualised volatility") +
  theme_ieee() + theme(legend.position = "bottom")
save_fig(p5, "fig5_vol_forecast.png", 3.5, 2.2)

# Beeswarm-style global Shapley
lo <- do.call(rbind, lapply(names(imp), function(nm) {
  v <- Xe[[nm]]; data.frame(feat = nm, phi = sh$phi[, nm], val = rank(v) / length(v))
}))
lo$feat <- factor(lo$feat, levels = rev(names(imp)))
p6 <- ggplot(lo[lo$feat %in% names(imp)[1:10], ], aes(phi, feat, colour = val)) +
  geom_vline(xintercept = 0, linewidth = .3) + geom_jitter(height = .18, size = .6, alpha = .8) +
  scale_colour_gradient(low = "#2166ac", high = "#b2182b", name = "Feature\nvalue", breaks = c(0, 1), labels = c("Low", "High")) +
  labs(x = "Shapley value (log-vol units)", y = NULL) + theme_ieee()
save_fig(p6, "fig6_shap_beeswarm.png", 3.5, 2.8)

top8 <- names(imp)[1:8]
p7 <- ggplot(subset(tab_shap_reg, Feature %in% top8), aes(MeanAbsShapley, factor(Feature, levels = rev(top8)), fill = factor(Regime, REG))) +
  geom_col(position = position_dodge(.8), width = .75) +
  scale_fill_manual(values = c(Calm = "#92c5de", Stress = "#fdb863", Crisis = "#d7301f"), name = "Regime") +
  labs(x = "Mean |Shapley value|", y = NULL) + theme_ieee()
save_fig(p7, "fig7_shap_by_regime.png", 3.5, 2.6)

# Local explanations: most extreme high-vol and low-vol forecast instances
ix_hi <- which.max(sh$fx); ix_lo <- which.min(sh$fx)
save_fig(plot_waterfall(sh$phi[ix_hi, ], sh$base[ix_hi], sh$fx[ix_hi], "High-risk day"), "fig8_local_high.png", 3.5, 2.2)
save_fig(plot_waterfall(sh$phi[ix_lo, ], sh$base[ix_lo], sh$fx[ix_lo], "Calm day"), "fig9_local_low.png", 3.5, 2.2)
local_info <- data.frame(case = c("High-risk", "Calm"), date = as.character(d$date[oos][pick][c(ix_hi, ix_lo)]),
                         true_regime = if (all(is.na(d$regime))) NA else REG[d$regime[oos][pick][c(ix_hi, ix_lo)]],
                         forecast_vol = exp(sh$fx[c(ix_hi, ix_lo)]))
save_tab(local_info, "shap_local_cases.csv")

# ============================== MONTE CARLO (many markets) ===================
mc_one <- function(seed) {
  d <- simulate_market(CFG$n_days, seed); f <- build_features(d); rule <- rule_regime(f)
  wf <- walk_forward(f); W <- build_weights(f, d, wf, rule); bt <- run_backtest(d, wf, W)
  tp <- perf_table(bt, W); tp$seed <- seed
  hy <- bt$ret[["Hybrid-Cap"]]; tp$HybridCapBeatsSharpe <- sapply(tp$Strategy, function(k) sharpe(hy) > sharpe(bt$ret[[k]]))
  tp
}
if (CFG$data_source == "simulated" && CFG$n_mc > 0) {
  cat("Monte Carlo over", CFG$n_mc, "independent markets on", CFG$cores, "cores ...\n")
  mc <- do.call(rbind, mclapply(1000 + seq_len(CFG$n_mc), function(sd) {
    fn <- file.path(CFG$out_dir, "tables", sprintf("mc_seed_%d.rds", sd))   # checkpoint per path
    if (file.exists(fn)) return(readRDS(fn)); r <- mc_one(sd); saveRDS(r, fn); r }, mc.cores = CFG$cores))
  save_tab(mc, "monte_carlo_raw.csv")
  agg <- do.call(rbind, lapply(split(mc, mc$Strategy), function(g) data.frame(
    Strategy = g$Strategy[1], Sharpe_mean = mean(g$Sharpe), Sharpe_sd = sd(g$Sharpe),
    MaxDD_mean = mean(g$MaxDD), MaxDD_sd = sd(g$MaxDD), Ret_mean = mean(g$AnnReturn),
    Vol_mean = mean(g$AnnVol), Calmar_mean = mean(g$Calmar),
    HybridCapWinRate_Sharpe = mean(g$HybridCapBeatsSharpe), N = nrow(g))))
  save_tab(agg, "monte_carlo_summary.csv")
  print(agg, digits = 3)
}

# ============================== CONSOLE REPORT ===============================
options(width = 160)
cat("\n--- Out-of-sample performance (seed", CFG$seed, ") ---\n"); print(tab_perf, digits = 3, row.names = FALSE)
cat("\n--- Bootstrap Sharpe tests (Hybrid-Cap minus benchmark) ---\n"); print(tab_boot, digits = 3, row.names = FALSE)
cat("\n--- Volatility forecast accuracy ---\n"); print(tab_fc, digits = 4, row.names = FALSE)
if (!is.null(tab_reg)) { cat("\n--- Rule-layer regime diagnostics ---\n"); print(tab_reg, digits = 3); print(tab_conf) }
cat("\n--- Shapley global importance ---\n"); print(head(tab_shap, 10), digits = 3, row.names = FALSE)
print(tab_shap_meta, digits = 3)
cat("\nTotal runtime:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min.  Outputs in ./", CFG$out_dir, "\n")

