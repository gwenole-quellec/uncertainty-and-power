#!/usr/bin/env python3

"""
Uncertainty analysis for sensitivity and specificity.

This script reproduces Figure 1 from:
"When ten cases are not enough: Uncertainty and statistical power in AI
screening evaluations"

Sensitivity and specificity are binomial proportions whose precision depends
on the metric-specific effective sample size: the number of positive cases for
sensitivity and the number of negative cases for specificity.

For each performance level and effective sample size, repeated binomial samples
are generated and exact two-sided Clopper--Pearson confidence intervals are
computed. Figure 1 reports their mean width as a function of sample size.

Author: Gwenolé Quellec
"""

import numpy as np
import pandas as pd

from scipy.stats import beta
from tqdm import tqdm

from plotting import plot_curves


# ============================================================================
# Configuration
# ============================================================================

ALPHA = 0.05
N_SIMULATIONS = 10000
RANDOM_STATE = 42

N_MIN = 1
N_MAX = 100
N_VALUES = range(N_MIN, N_MAX + 1)

PERFORMANCE_LEVELS = (
    0.95,
    0.90,
    0.85,
    0.80,
)

FIGURE_FILENAME = "fig1.pdf"

X_LABEL = "Metric-specific effective sample size"
Y_LABEL = "Mean width of exact 95% confidence interval"

X_TICKS = (10, 5)


# ============================================================================
# Statistical functions
# ============================================================================

def clopper_pearson_ci(
    k: int,
    n: int,
    alpha: float,
) -> tuple[float, float]:
    """
    Compute the exact two-sided Clopper--Pearson confidence interval.

    Parameters
    ----------
    k : int
        Number of successes.
    n : int
        Number of Bernoulli trials.
    alpha : float
        Significance level.

    Returns
    -------
    lower, upper : tuple of float
        Lower and upper confidence limits.
    """

    if k == 0:
        lower = 0.0
    else:
        lower = beta.ppf(
            alpha / 2,
            k,
            n - k + 1,
        )

    if k == n:
        upper = 1.0
    else:
        upper = beta.ppf(
            1 - alpha / 2,
            k + 1,
            n - k,
        )

    return lower, upper


def simulate_binomial_uncertainty(
    true_performance: float,
    n_values,
    n_simulations: int,
    alpha: float,
    random_state: int,
) -> pd.DataFrame:
    """
    Estimate mean exact confidence-interval width by simulation.

    For each sample size, repeated binomial observations are generated under
    the specified true performance. Exact Clopper--Pearson intervals are then
    computed and their widths averaged.

    Parameters
    ----------
    true_performance : float
        Underlying binomial success probability.
    n_values : iterable of int
        Effective sample sizes to evaluate.
    n_simulations : int
        Number of simulated evaluations per sample size.
    alpha : float
        Significance level.
    random_state : int
        Random seed.

    Returns
    -------
    pandas.DataFrame
        Mean confidence-interval width for each effective sample size.
    """

    rng = np.random.default_rng(random_state)
    rows = []

    for n in tqdm(
        n_values,
        desc=f"CI simulation (p={true_performance:.2f})",
    ):
        observed = rng.binomial(
            n=n,
            p=true_performance,
            size=n_simulations,
        )

        for k in observed:
            ci_low, ci_high = clopper_pearson_ci(
                k,
                n,
                alpha,
            )

            rows.append({
                "n": n,
                "ci_width": ci_high - ci_low,
            })

    df = pd.DataFrame(rows)

    return (
        df.groupby(
            "n",
            as_index=False,
        )["ci_width"]
        .mean()
    )


# ============================================================================
# Main analysis
# ============================================================================

def main():

    curves = {}

    for i, performance in enumerate(PERFORMANCE_LEVELS):

        df = simulate_binomial_uncertainty(
            true_performance=performance,
            n_values=N_VALUES,
            n_simulations=N_SIMULATIONS,
            alpha=ALPHA,
            random_state=RANDOM_STATE + i,
        )

        label = f"True performance = {int(performance * 100)}%"

        curves[label] = (
            df["n"],
            df["ci_width"],
        )

    plot_curves(
        curves,
        xlabel=X_LABEL,
        ylabel=Y_LABEL,
        filename=FIGURE_FILENAME,
        xmax=N_MAX,
        xticks=X_TICKS,
    )


if __name__ == "__main__":
    main()
