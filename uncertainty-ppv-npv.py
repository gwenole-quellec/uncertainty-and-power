#!/usr/bin/env python3

"""
Uncertainty analysis for positive and negative predictive values.

This script reproduces Figure 2 from:
"When ten cases are not enough: Uncertainty and statistical power in AI
screening evaluations"

PPV and NPV are binomial proportions whose effective denominators depend on
the numbers of positive and negative predictions, respectively. These
denominators are themselves random variables when the total validation sample
size is fixed.

For each system and total sample size, the expected width of the exact
two-sided Clopper--Pearson confidence interval is computed by averaging first
over all possible numbers of correct predictions conditional on the relevant
predictive-value denominator, and then over the distribution of that
denominator. The empirical confusion matrices observed on Messidor-2 are used
to estimate the predictive values and the probabilities of positive and
negative predictions.

Author: Gwenolé Quellec
"""

import numpy as np

from scipy.stats import beta, binom
from tqdm import tqdm

from plotting import plot_curves


# ============================================================================
# Configuration
# ============================================================================

ALPHA = 0.05

N_MIN = 1
N_MAX = 400
N_VALUES = range(N_MIN, N_MAX + 1)

N_MESSIDOR = 874

SYSTEMS = {
    "OphtAI": {
        "tp": 181,
        "fp": 54,
        "tn": 630,
        "fn": 9,
    },
    "RETFound": {
        "tp": 144,
        "fp": 137,
        "tn": 547,
        "fn": 46,
    },
    "ORDER-DR": {
        "tp": 116,
        "fp": 17,
        "tn": 667,
        "fn": 74,
    },
}

FIGURE_FILENAME_PPV = "fig2a.pdf"
FIGURE_FILENAME_NPV = "fig2b.pdf"

X_LABEL = "Total sample size (N)"
Y_LABEL = "Mean width of exact 95% confidence interval"

X_TICKS = (100, 20)


# ============================================================================
# Statistical functions
# ============================================================================

def predictive_values_from_counts(
    tp: int,
    fp: int,
    tn: int,
    fn: int,
) -> tuple[float, float, float, float]:
    """
    Compute PPV, NPV, and the probabilities of positive and negative predictions.

    Parameters
    ----------
    tp, fp, tn, fn : int
        Confusion-matrix counts.

    Returns
    -------
    ppv, npv, q_pos, q_neg : tuple of float
        Positive predictive value, negative predictive value, probability of
        a positive prediction, and probability of a negative prediction.
    """

    n = tp + fp + tn + fn

    ppv = tp / (tp + fp)
    npv = tn / (tn + fn)

    q_pos = (tp + fp) / n
    q_neg = (tn + fn) / n

    return ppv, npv, q_pos, q_neg


def precompute_clopper_pearson_widths(
    n_max: int,
    alpha: float,
) -> list:
    """
    Precompute exact Clopper--Pearson interval widths for all sample sizes.

    Parameters
    ----------
    n_max : int
        Maximum binomial denominator.
    alpha : float
        Significance level.

    Returns
    -------
    list
        Interval-width arrays indexed by sample size.
    """

    print("Precomputing Clopper--Pearson intervals...")

    widths = [None] * (n_max + 1)

    for n in tqdm(
        range(1, n_max + 1),
        desc="Clopper--Pearson intervals",
    ):
        k = np.arange(n + 1)

        lower = np.zeros(n + 1)
        upper = np.ones(n + 1)

        mask = k > 0
        lower[mask] = beta.ppf(
            alpha / 2,
            k[mask],
            n - k[mask] + 1,
        )

        mask = k < n
        upper[mask] = beta.ppf(
            1 - alpha / 2,
            k[mask] + 1,
            n - k[mask],
        )

        widths[n] = upper - lower

    return widths


def expected_width_by_denominator(
    performance: float,
    interval_widths: list,
) -> np.ndarray:
    """
    Compute expected exact confidence-interval width for each denominator.

    Parameters
    ----------
    performance : float
        Underlying PPV or NPV.
    interval_widths : list
        Precomputed Clopper--Pearson interval widths.

    Returns
    -------
    numpy.ndarray
        Expected interval width for each possible denominator.
    """

    n_max = len(interval_widths) - 1

    expected_widths = np.zeros(n_max + 1)
    expected_widths[0] = np.nan

    for n in range(1, n_max + 1):
        k = np.arange(n + 1)

        expected_widths[n] = np.sum(
            binom.pmf(k, n, performance)
            * interval_widths[n]
        )

    return expected_widths


def expected_width_vs_total_sample_size(
    n_values,
    prediction_probability: float,
    expected_widths: np.ndarray,
) -> np.ndarray:
    """
    Average interval width over the random predictive-value denominator.

    For a total validation sample size N, the denominator of PPV or NPV is
    modeled as M ~ Binomial(N, q), where q is the probability of a positive
    or negative prediction.

    Parameters
    ----------
    n_values : iterable of int
        Total validation sample sizes.
    prediction_probability : float
        Probability of a positive or negative prediction.
    expected_widths : numpy.ndarray
        Expected interval width conditional on each denominator.

    Returns
    -------
    numpy.ndarray
        Expected confidence-interval width for each total sample size.
    """

    results = []

    for n in n_values:
        m = np.arange(1, n + 1)

        weights = binom.pmf(
            m,
            n,
            prediction_probability,
        )

        p_nonzero = (
            1.0
            - binom.pmf(
                0,
                n,
                prediction_probability,
            )
        )

        width = np.sum(
            weights * expected_widths[m]
        ) / p_nonzero

        results.append(width)

    return np.asarray(results)


# ============================================================================
# Main analysis
# ============================================================================

def main():

    interval_widths = precompute_clopper_pearson_widths(
        n_max=N_MAX,
        alpha=ALPHA,
    )

    ppv_curves = {}
    npv_curves = {}

    for name, counts in tqdm(
        SYSTEMS.items(),
        desc="Systems",
    ):
        n = (
            counts["tp"]
            + counts["fp"]
            + counts["tn"]
            + counts["fn"]
        )

        if n != N_MESSIDOR:
            raise ValueError(
                f"{name}: confusion matrix sums to {n}, "
                f"expected {N_MESSIDOR}"
            )

        ppv, npv, q_pos, q_neg = predictive_values_from_counts(
            tp=counts["tp"],
            fp=counts["fp"],
            tn=counts["tn"],
            fn=counts["fn"],
        )

        print(
            f"{name}: "
            f"PPV={ppv:.4f}, "
            f"NPV={npv:.4f}, "
            f"P(pred+)={q_pos:.4f}, "
            f"P(pred-)={q_neg:.4f}"
        )

        ppv_expected_widths = expected_width_by_denominator(
            performance=ppv,
            interval_widths=interval_widths,
        )

        npv_expected_widths = expected_width_by_denominator(
            performance=npv,
            interval_widths=interval_widths,
        )

        ppv_curves[name] = (
            N_VALUES,
            expected_width_vs_total_sample_size(
                n_values=N_VALUES,
                prediction_probability=q_pos,
                expected_widths=ppv_expected_widths,
            ),
        )

        npv_curves[name] = (
            N_VALUES,
            expected_width_vs_total_sample_size(
                n_values=N_VALUES,
                prediction_probability=q_neg,
                expected_widths=npv_expected_widths,
            ),
        )

    plot_curves(
        ppv_curves,
        xlabel=X_LABEL,
        ylabel=Y_LABEL,
        filename=FIGURE_FILENAME_PPV,
        xmax=N_MAX,
        xticks=X_TICKS,
    )

    plot_curves(
        npv_curves,
        xlabel=X_LABEL,
        ylabel=Y_LABEL,
        filename=FIGURE_FILENAME_NPV,
        xmax=N_MAX,
        xticks=X_TICKS,
    )


if __name__ == "__main__":
    main()
