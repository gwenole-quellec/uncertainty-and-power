#!/usr/bin/env python3

"""
Power analysis for positive and negative predictive value comparisons.

This script reproduces Figure 4 from:
"When ten cases are not enough: Uncertainty and statistical power in AI
screening evaluations"

PPV and NPV comparisons between paired systems depend on the full joint
distribution of disease status and paired binary predictions. For each system
pair, the probabilities of the eight possible joint outcomes are estimated
empirically from Messidor-2.

For each total validation sample size, repeated paired datasets are generated
from the corresponding multinomial distribution. Differences in PPV and NPV
are tested using the generalized score test for paired predictive values
proposed by Leisenring et al. (2000), and statistical power is estimated as
the proportion of simulations rejecting the null hypothesis of equal
predictive values.

Author: Gwenolé Quellec
"""

import numpy as np
import pandas as pd

from scipy.stats import chi2
from tqdm import tqdm

from plotting import plot_curves


# ============================================================================
# Configuration
# ============================================================================

ALPHA = 0.05
N_SIMULATIONS = 10000
RANDOM_STATE = 42

N_MIN = 1
N_MAX = 2000
N_VALUES = range(N_MIN, N_MAX + 1)

SCENARIOS = {
    "OphtAI vs RETFound": (
        140, 41, 4, 5,
        21, 33, 116, 514,
    ),
    "OphtAI vs ORDER-DR": (
        116, 65, 0, 9,
        7, 47, 10, 620,
    ),
    "RETFound vs ORDER-DR": (
        100, 44, 16, 30,
        12, 125, 5, 542,
    ),
}

FIGURE_FILENAME_PPV = "fig4a.pdf"
FIGURE_FILENAME_NPV = "fig4b.pdf"

X_LABEL = "Total sample size (N)"
Y_LABEL = "Statistical power"

X_MAX_PPV = 800
X_MAX_NPV = 1600

X_TICKS_PPV = (200, 50)
X_TICKS_NPV = (200, 50)

POWER_REFERENCE = 0.80


# ============================================================================
# Statistical functions
# ============================================================================

def predictive_values(
    counts,
) -> tuple[float, float, float, float]:
    """
    Compute PPV and NPV for both systems from an eight-cell paired table.

    The cell order is:
        d11, d10, d01, d00,
        n11, n10, n01, n00

    Here, d and n denote rDR-positive and rDR-negative examinations,
    respectively. The first binary index is the prediction of system A and
    the second binary index is the prediction of system B.

    Parameters
    ----------
    counts : array-like of int
        Empirical counts for the eight joint outcomes.

    Returns
    -------
    ppv_a, ppv_b, npv_a, npv_b : tuple of float
        Positive and negative predictive values of both systems.
    """

    (
        d11, d10, d01, d00,
        n11, n10, n01, n00,
    ) = np.asarray(counts, dtype=float)

    pos_a = d11 + d10 + n11 + n10
    pos_b = d11 + d01 + n11 + n01

    ppv_a = (d11 + d10) / pos_a
    ppv_b = (d11 + d01) / pos_b

    neg_a = d01 + d00 + n01 + n00
    neg_b = d10 + d00 + n10 + n00

    npv_a = (n01 + n00) / neg_a
    npv_b = (n10 + n00) / neg_b

    return ppv_a, ppv_b, npv_a, npv_b


def generalized_score_pvalues(
    draws,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute generalized score-test p-values for paired PPV and NPV.

    This is a vectorized implementation of the generalized score statistic
    proposed by Leisenring, Alonzo and Pepe (2000).

    The cell order is:
        d11, d10, d01, d00,
        n11, n10, n01, n00

    Parameters
    ----------
    draws : array-like, shape (n_simulations, 8)
        Paired eight-cell tables.

    Returns
    -------
    p_ppv, p_npv : tuple of numpy.ndarray
        P-values for equality of PPV and NPV, respectively.

    Notes
    -----
    Non-estimable or degenerate cases are treated as non-rejections
    (p = 1).
    """

    x = np.asarray(draws, dtype=float)

    (
        d11, d10, d01, d00,
        n11, n10, n01, n00,
    ) = [x[:, i] for i in range(8)]

    n_simulations = x.shape[0]

    p_ppv = np.ones(n_simulations)
    p_npv = np.ones(n_simulations)

    # PPV
    pos_a = d11 + d10 + n11 + n10
    pos_b = d11 + d01 + n11 + n01

    total_pos = pos_a + pos_b

    valid = (
        (pos_a > 0)
        & (pos_b > 0)
        & (total_pos > 0)
    )

    zbar = np.zeros(n_simulations)
    dbar = np.zeros(n_simulations)

    zbar[valid] = pos_b[valid] / total_pos[valid]

    dbar[valid] = (
        2 * d11[valid]
        + d10[valid]
        + d01[valid]
    ) / total_pos[valid]

    numerator = (
        d11 * (1 - 2 * zbar)
        + d01 * (1 - zbar)
        - d10 * zbar
    ) ** 2

    denominator = (
        (1 - dbar) ** 2
        * (
            d11 * (1 - 2 * zbar) ** 2
            + d01 * (1 - zbar) ** 2
            + d10 * zbar ** 2
        )
        + dbar ** 2
        * (
            n11 * (1 - 2 * zbar) ** 2
            + n01 * (1 - zbar) ** 2
            + n10 * zbar ** 2
        )
    )

    valid &= (
        np.isfinite(denominator)
        & (denominator > 0)
    )

    statistic = np.zeros(n_simulations)
    statistic[valid] = numerator[valid] / denominator[valid]

    p_ppv[valid] = chi2.sf(
        statistic[valid],
        df=1,
    )

    # NPV
    neg_a = d01 + d00 + n01 + n00
    neg_b = d10 + d00 + n10 + n00

    total_neg = neg_a + neg_b

    valid = (
        (neg_a > 0)
        & (neg_b > 0)
        & (total_neg > 0)
    )

    zbar = np.zeros(n_simulations)
    dbar = np.zeros(n_simulations)

    zbar[valid] = neg_b[valid] / total_neg[valid]

    dbar[valid] = (
        2 * n00[valid]
        + n10[valid]
        + n01[valid]
    ) / total_neg[valid]

    numerator = (
        n00 * (1 - 2 * zbar)
        + n10 * (1 - zbar)
        - n01 * zbar
    ) ** 2

    denominator = (
        (1 - dbar) ** 2
        * (
            n00 * (1 - 2 * zbar) ** 2
            + n10 * (1 - zbar) ** 2
            + n01 * zbar ** 2
        )
        + dbar ** 2
        * (
            d00 * (1 - 2 * zbar) ** 2
            + d10 * (1 - zbar) ** 2
            + d01 * zbar ** 2
        )
    )

    valid &= (
        np.isfinite(denominator)
        & (denominator > 0)
    )

    statistic = np.zeros(n_simulations)
    statistic[valid] = numerator[valid] / denominator[valid]

    p_npv[valid] = chi2.sf(
        statistic[valid],
        df=1,
    )

    return p_ppv, p_npv


def simulate_power(
    counts,
    n_values,
    n_simulations: int,
    alpha: float,
    random_state: int,
    desc: str,
) -> pd.DataFrame:
    """
    Estimate power for paired PPV and NPV comparisons by simulation.

    Parameters
    ----------
    counts : array-like of int
        Empirical counts for the eight joint outcomes.
    n_values : iterable of int
        Total validation sample sizes to evaluate.
    n_simulations : int
        Number of simulated paired datasets per sample size.
    alpha : float
        Significance level.
    random_state : int
        Random seed.
    desc : str
        Description displayed by the progress bar.

    Returns
    -------
    pandas.DataFrame
        Estimated PPV and NPV power for each total sample size.
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

        p_ppv, p_npv = generalized_score_pvalues(draws)

        results.append({
            "n": n,
            "ppv_power": np.mean(p_ppv < alpha),
            "npv_power": np.mean(p_npv < alpha),
        })

    return pd.DataFrame(results)


def required_n(
    df: pd.DataFrame,
    column: str,
    target_power: float,
):
    """
    Return the first simulated sample size reaching a target power.

    Parameters
    ----------
    df : pandas.DataFrame
        Simulated power curves.
    column : str
        Column containing the power estimates.
    target_power : float
        Target statistical power.

    Returns
    -------
    int or None
        Smallest evaluated sample size reaching the target, or None if the
        target is not reached.
    """

    matches = df.loc[df[column] >= target_power]

    if len(matches) == 0:
        return None

    return int(matches.iloc[0]["n"])


# ============================================================================
# Analysis functions
# ============================================================================

def run_analysis():
    """
    Run and plot the paired PPV and NPV power analyses.
    """

    ppv_curves = {}
    npv_curves = {}

    print()
    print("=" * 72)
    print("PPV / NPV power analysis")
    print("=" * 72)

    for i, (label, counts) in enumerate(SCENARIOS.items()):

        counts = np.asarray(counts)

        ppv_a, ppv_b, npv_a, npv_b = predictive_values(counts)

        observed_ppv_p, observed_npv_p = generalized_score_pvalues(
            counts[np.newaxis, :]
        )

        print()
        print(label)
        print("-" * len(label))

        print(
            "Empirical counts: "
            + ", ".join(str(x) for x in counts)
        )

        print(f"Total N: {counts.sum()}")

        print(
            f"PPV A: {ppv_a:.4f}  "
            f"PPV B: {ppv_b:.4f}  "
            f"Difference A-B: {ppv_a - ppv_b:+.4f}"
        )

        print(
            f"NPV A: {npv_a:.4f}  "
            f"NPV B: {npv_b:.4f}  "
            f"Difference A-B: {npv_a - npv_b:+.4f}"
        )

        print(
            f"Observed generalized-score p-value (PPV): "
            f"{observed_ppv_p[0]:.6g}"
        )

        print(
            f"Observed generalized-score p-value (NPV): "
            f"{observed_npv_p[0]:.6g}"
        )

        df = simulate_power(
            counts=counts,
            n_values=N_VALUES,
            n_simulations=N_SIMULATIONS,
            alpha=ALPHA,
            random_state=RANDOM_STATE + i,
            desc=label,
        )

        print()

        ppv_n = required_n(
            df,
            "ppv_power",
            POWER_REFERENCE,
        )

        npv_n = required_n(
            df,
            "npv_power",
            POWER_REFERENCE,
        )

        print(
            f"PPV: N for {POWER_REFERENCE:.0%} power: "
            f"{ppv_n}"
        )

        print(
            f"NPV: N for {POWER_REFERENCE:.0%} power: "
            f"{npv_n}"
        )

        ppv_curves[label] = (
            df["n"],
            df["ppv_power"],
        )

        npv_curves[label] = (
            df["n"],
            df["npv_power"],
        )

    plot_curves(
        ppv_curves,
        xlabel=X_LABEL,
        ylabel=Y_LABEL,
        filename=FIGURE_FILENAME_PPV,
        xmax=X_MAX_PPV,
        xticks=X_TICKS_PPV,
        hline=POWER_REFERENCE,
    )

    plot_curves(
        npv_curves,
        xlabel=X_LABEL,
        ylabel=Y_LABEL,
        filename=FIGURE_FILENAME_NPV,
        xmax=X_MAX_NPV,
        xticks=X_TICKS_NPV,
        hline=POWER_REFERENCE,
    )


# ============================================================================
# Main analysis
# ============================================================================

def main():

    run_analysis()


if __name__ == "__main__":
    main()
