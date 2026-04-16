"""
Simulation of uncertainty and statistical power in small-sample evaluations.

This script reproduces the figures from:
"When ten cases are not enough: Uncertainty and statistical power in AI screening evaluations"

It simulates:
1) The uncertainty of sensitivity estimates (Clopper–Pearson CI)
2) The statistical power of McNemar's test for comparing models

Author: Gwenolé Quellec
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from scipy.stats import beta
from statsmodels.stats.contingency_tables import mcnemar
from tqdm import tqdm


# =============================================================================
# Confidence interval
# =============================================================================

def clopper_pearson_ci(k: int, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Exact Clopper–Pearson confidence interval for a binomial proportion."""
    if k == 0:
        lower = 0.0
    else:
        lower = beta.ppf(alpha / 2, k, n - k + 1)

    if k == n:
        upper = 1.0
    else:
        upper = beta.ppf(1 - alpha / 2, k + 1, n - k)

    return lower, upper


# =============================================================================
# Sensitivity instability simulation
# =============================================================================

def simulate_sensitivity_instability(
    true_sensitivity: float,
    n_values=range(1, 101),
    n_simulations: int = 10000,
    alpha: float = 0.05,
    random_state: int = 42,
) -> pd.DataFrame:
    """Simulate observed sensitivity and confidence intervals."""
    rng = np.random.default_rng(random_state)
    rows = []

    for n in tqdm(n_values, desc=f"CI simulation (p={true_sensitivity:.2f})"):

        detected = rng.binomial(n=n, p=true_sensitivity, size=n_simulations)

        for k in detected:
            ci_low, ci_high = clopper_pearson_ci(k, n, alpha)
            rows.append({
                "n_positive": n,
                "ci_width": ci_high - ci_low
            })

    df = pd.DataFrame(rows)

    return (
        df.groupby("n_positive", as_index=False)["ci_width"]
        .mean()
        .rename(columns={"ci_width": "ci_width"})
    )


# =============================================================================
# McNemar power simulation
# =============================================================================

def simulate_mcnemar_power(
    sens_A: float,
    sens_B: float,
    n_values=range(5, 101),
    n_simulations: int = 10000,
    alpha: float = 0.05,
    random_state: int = 42
) -> pd.DataFrame:
    """Estimate statistical power using McNemar's test."""
    rng = np.random.default_rng(random_state)
    results = []

    for n in tqdm(n_values, desc=f"Power ({sens_A:.2f} vs {sens_B:.2f})"):

        significant = 0

        for _ in range(n_simulations):

            A = rng.binomial(1, sens_A, size=n)
            B = rng.binomial(1, sens_B, size=n)

            table = [
                [np.sum((A == 1) & (B == 1)), np.sum((A == 1) & (B == 0))],
                [np.sum((A == 0) & (B == 1)), np.sum((A == 0) & (B == 0))]
            ]

            try:
                pval = mcnemar(table, exact=True).pvalue
                if pval < alpha:
                    significant += 1
            except:
                pass

        results.append({
            "n_positive": n,
            "power": significant / n_simulations
        })

    return pd.DataFrame(results)


# =============================================================================
# Plot utilities
# =============================================================================

def smooth_series(y, window=5):
    return pd.Series(y).rolling(window=window, center=True, min_periods=1).mean()


def setup_figure(figsize=(5, 5)):
    """Configure matplotlib figure style."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 9,
    })

    plt.figure(figsize=figsize)

    ax = plt.gca()
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1)

    ax.axvline(10, linestyle="--", color="black", alpha=0.6)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))

    ax.grid(which='both', linestyle='--', alpha=0.3)


def plot_curves(curves: dict, ylabel: str, filename: str, smooth=False, window=5, hline=None):
    """Generic plotting function."""
    setup_figure()

    for label, (x, y) in curves.items():
        if smooth:
            y = smooth_series(y, window)
        plt.plot(x, y, label=label)

    if hline is not None:
        plt.axhline(hline, linestyle="--", color="black", alpha=0.6)

    plt.xlabel("Number of positive cases (n)")
    plt.ylabel(ylabel)
    plt.legend()

    plt.tight_layout()
    plt.savefig(filename)
    plt.show()


# =============================================================================
# Main
# =============================================================================

def main():
    # --- CI simulations ---
    sensitivities = [0.95, 0.90, 0.85, 0.80]
    ci_results = {}

    for s in sensitivities:
        df = simulate_sensitivity_instability(s)
        ci_results[f"True sensitivity = {int(s*100)}%"] = (
            df["n_positive"],
            df["ci_width"]
        )

    plot_curves(
        ci_results,
        ylabel="95% confidence interval width",
        filename="fig1.pdf",
        smooth=True
    )

    # --- Power simulations ---
    comparisons = [(0.9, 0.8), (0.9, 0.7), (0.8, 0.7)]
    power_results = {}

    for a, b in comparisons:
        df = simulate_mcnemar_power(a, b)
        label = f"{int(a*100)}% vs {int(b*100)}%"
        power_results[label] = (df["n_positive"], df["power"])

    plot_curves(
        power_results,
        ylabel="Statistical power",
        filename="fig2.pdf",
        smooth=True,
        hline=0.8
    )


if __name__ == "__main__":
    main()
