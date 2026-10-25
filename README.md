# RegimeShap-Quant

**A Hybrid Rule-Based and Gradient-Boosting Framework with Shapley Explainability for Regime-Aware Volatility-Targeted Allocation**

## Overview

RegimeShap-Quant is an R-based quantitative finance research framework designed to investigate whether combining market-regime rules, machine-learning volatility forecasts, and explainable AI improves volatility-targeted portfolio allocation.

The framework combines a rule-based regime detection layer with a Gradient Boosting Machine (GBM) volatility forecaster. A conservative hybrid decision layer uses these signals to determine portfolio exposure, while Monte Carlo Shapley values help explain the model's volatility forecasts.

The project focuses not only on portfolio risk management but also on understanding when a model succeeds and when it fails.

## Key Features

- **Market regime detection:** Identifies Calm, Stress, and Crisis conditions using volatility and drawdown indicators.
- **Machine-learning forecasting:** Uses gradient boosting to forecast forward-looking volatility.
- **Hybrid decision framework:** Combines rule-based signals and machine-learning forecasts.
- **Volatility targeting:** Adjusts portfolio exposure according to estimated risk.
- **Explainable AI:** Uses Monte Carlo Shapley values to investigate feature contributions.
- **Walk-forward evaluation:** Evaluates forecasts using out-of-sample procedures.
- **Transaction-cost modelling:** Incorporates transaction costs and a no-trade band.
- **Statistical evaluation:** Includes block-bootstrap comparisons, forecast error metrics, and robustness experiments.
- **Reproducible experiments:** Uses configurable parameters and random seeds.

## Research Objectives

1. Investigate whether machine-learning volatility forecasts improve on conventional trailing-volatility targeting.
2. Evaluate whether regime-based rules provide additional value during stressed market conditions.
3. Compare hybrid allocation with alternative volatility-targeting strategies.
4. Use Shapley explanations to identify influential features and diagnose forecast limitations.
5. Examine the effect of transaction costs and model parameters on portfolio performance.

## Methodology

The framework follows this general pipeline:

1. Generate simulated market data or retrieve historical market data.
2. Calculate financial features from asset returns.
3. Identify market regimes using rule-based indicators.
4. Forecast future volatility using gradient boosting.
5. Combine forecasts and regime signals through the allocation decision layer.
6. Construct volatility-targeted portfolio weights.
7. Run the backtest with transaction costs and trading constraints.
8. Evaluate risk-adjusted returns, drawdowns, and forecast accuracy.
9. Analyse feature contributions using Monte Carlo Shapley values.
10. Perform robustness and statistical comparisons.

## Technology Stack

- **Programming language:** R
- **Machine learning:** Gradient boosting
- **Data manipulation and analysis:** R data frames and supporting
