#!/usr/bin/env python3

"""
Power analysis for sensitivity and specificity comparisons.

This script reproduces Figure 3 from:
"When ten cases are not enough: Uncertainty and statistical power in AI
screening evaluations"

Sensitivity and specificity comparisons between paired systems depend not only
on their marginal performance but also on their pairwise disagreement
structure. For each system pair, the probabilities of the four paired outcomes
are estimated empirically from Messidor-2, separately among rDR-positive and
rDR-negative examinations.

For each class-specific sample size, repeated paired evaluations are generated
from the corresponding multinomial distribution. Differences between systems
are tested using the exact two-sided McNemar test, and statistical power is
estimated as the proportion of simulations rejecting the null hypothesis of
equal performance.

Author: Gwenolé Quellec
"""

import numpy as np
import pandas as pd

from scipy.stats import binom
from tqdm import tqdm

from plotting import plot_curves


# ============================================================================
# Configuration
# ============================================================================

ALPHA = 0.05
N_SIMULATIONS = 10000
RANDOM_STATE = 42

N_MIN = 1
N_MAX = 300
N_VALUES = range(N_MIN, N_MAX + 1)

SENSITIVITY_SCENARIOS = {
    "OphtAI vs RETFound": (
        140,  # both correct
        41,   # OphtAI correct, RETFound incorrect
        4,    # OphtAI incorrect, RETFound correct
        5,    # both incorrect
    ),
    "OphtAI vs ORDER-DR": (
        116,
        65,
        0,
        9,
    ),
    "RETFound vs ORDER-DR": (
        100,
        44,
        16,
        30,
    ),
}

SPECIFICITY_SCENARIOS = {
    "OphtAI vs RETFound": (
        514,  # both correct
        116,  # OphtAI correct, RETFound incorrect
        33,   # OphtAI incorrect, RETFound correct
        21,   # both incorrect
    ),
    "OphtAI vs ORDER-DR": (
        620,
        10,
        47,
        7,
    ),
    "RETFound vs ORDER-DR": (
        542,
        5,
        125,
        12,
    ),
}

FIGURE_FILENAME_SENSITIVITY = "fig3a.pdf"
FIGURE_FILENAME_SPECIFICITY = "fig3b.pdf"

X_LABEL_SENSITIVITY = "Number of rDR-positive cases"
X_LABEL_SPECIFICITY = "Number of rDR-negative cases"
Y_LABEL = "Statistical power"

X_TICKS = (100, 20)

POWER_REFERENCE = 0.80


# ============================================================================
# Statistical functions
# ============================================================================

def simulate_mcnemar_power(
    counts,
    n_values,
    n_simulations: int,
    alpha: float,
    random_state: int,
    desc: str,
) -> pd.DataFrame:
    """
    Estimate power of the exact paired McNemar test by simulation.

    Parameters
    ----------
    counts : array-like of int
        Empirical counts for the four paired outcomes, ordered as both correct,
        only system A correct, only system B correct, and both incorrect.
    n_values : iterable of int
        Class-specific sample sizes to evaluate.
    n_simulations : int
        Number of simulated paired evaluations per sample size.
    alpha : float
        Significance level.
    random_state : int
        Random seed.
    desc : str
        Description displayed by the progress bar.

    Returns
    -------
    pandas.DataFrame
        Estimated statistical power for each sample size.
    """

    counts = np.asarray(counts, dtype=float)
    probabilities = counts / counts.sum()

    rng = np.random.default_rng(random_state)
    results = []

    for n in tqdm(n_values, desc=desc):

        draws = rng.multinomial(
            n,
            probabilities,
            size=n_simulations,
        )

        n10 = draws[:, 1]
        n01 = draws[:, 2]

        discordant = n10 + n01
        minority = np.minimum(n10, n01)

        pvalues = np.ones(n_simulations)

        mask = discordant > 0

        pvalues[mask] = np.minimum(
            1.0,
            2.0 * binom.cdf(
                minority[mask],
                discordant[mask],
                0.5,
            ),
        )

        results.append({
            "n": n,
            "power": np.mean(pvalues < alpha),
        })

    return pd.DataFrame(results)


def required_n(
    df: pd.DataFrame,
    target_power: float,
):
    """
    Return the first simulated sample size reaching a target power.

    Parameters
    ----------
    df : pandas.DataFrame
        Simulated power curve.
    target_power : float
        Target statistical power.

    Returns
    -------
    int or None
        Smallest evaluated sample size reaching the target, or None if the
        target is not reached.
    """

    matches = df.loc[df["power"] >= target_power]

    if len(matches) == 0:
        return None

    return int(matches.iloc[0]["n"])


# ============================================================================
# Analysis functions
# ============================================================================

def run_analysis(
    scenarios,
    xlabel: str,
    filename: str,
):
    """
    Run and plot a paired power analysis for one performance metric.

    Parameters
    ----------
    scenarios : dict
        Empirical four-cell paired outcome tables for each system pair.
    xlabel : str
        Label of the horizontal axis.
    filename : str
        Output figure filename.
    """

    curves = {}

    print()
    print("=" * 72)
    print(filename)
    print("=" * 72)

    for i, (label, counts) in enumerate(scenarios.items()):

        counts = np.asarray(counts)
        probabilities = counts / counts.sum()

        n11, n10, n01, n00 = counts
        p11, p10, p01, p00 = probabilities

        performance_a = p11 + p10
        performance_b = p11 + p01
        discordance = p10 + p01

        print()
        print(label)
        print("-" * len(label))

        print(
            f"Empirical counts: "
            f"n11={n11}, n10={n10}, "
            f"n01={n01}, n00={n00}"
        )

        print(f"Performance A: {performance_a:.4f}")
        print(f"Performance B: {performance_b:.4f}")

        print(
            f"Difference A-B: "
            f"{performance_a - performance_b:+.4f}"
        )

        print(f"Discordance: {discordance:.4f}")

        df = simulate_mcnemar_power(
            counts=counts,
            n_values=N_VALUES,
            n_simulations=N_SIMULATIONS,
            alpha=ALPHA,
            random_state=RANDOM_STATE + i,
            desc=label,
        )

        n_required = required_n(
            df,
            POWER_REFERENCE,
        )

        print(
            f"N for {POWER_REFERENCE:.0%} power: "
            f"{n_required}"
        )

        curves[label] = (
            df["n"],
            df["power"],
        )

    plot_curves(
        curves,
        xlabel=xlabel,
        ylabel=Y_LABEL,
        filename=filename,
        xmax=N_MAX,
        xticks=X_TICKS,
        hline=POWER_REFERENCE,
    )


# ============================================================================
# Main analysis
# ============================================================================

def main():

    run_analysis(
        scenarios=SENSITIVITY_SCENARIOS,
        xlabel=X_LABEL_SENSITIVITY,
        filename=FIGURE_FILENAME_SENSITIVITY,
    )

    run_analysis(
        scenarios=SPECIFICITY_SCENARIOS,
        xlabel=X_LABEL_SPECIFICITY,
        filename=FIGURE_FILENAME_SPECIFICITY,
    )


if __name__ == "__main__":
    main()
