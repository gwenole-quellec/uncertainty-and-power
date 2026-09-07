# Experimental data and analysis parameters

Numerical values used to generate the figures in
"When ten cases are not enough: Uncertainty and statistical power
in AI screening evaluations".

All operating thresholds were fixed independently of Messidor-2.

---

## 1. Messidor-2

| Quantity | Value |
|---|---:|
| Total examinations | 874 |
| rDR-positive examinations | 190 |
| rDR-negative examinations | 684 |
| rDR prevalence | 0.2174 |

---

## 2. Systems

| System | Operating threshold |
|---|---:|
| OphtAI | -0.4 |
| RETFound | 0.487161 |
| ORDER-DR | 0.267 |

---

## 3. Confusion matrices

Counts at the prespecified operating thresholds.

| System | TP | FN | TN | FP |
|---|---:|---:|---:|---:|
| OphtAI | 181 | 9 | 630 | 54 |
| RETFound | 144 | 46 | 547 | 137 |
| ORDER-DR | 116 | 74 | 667 | 17 |

Derived performance:

| System | Sensitivity | Specificity | PPV | NPV |
|---|---:|---:|---:|---:|
| OphtAI | 0.9526 | 0.9211 | 0.7702 | 0.9859 |
| RETFound | 0.7579 | 0.7997 | 0.5125 | 0.9224 |
| ORDER-DR | 0.6105 | 0.9751 | 0.8722 | 0.9001 |

---

## 4. Pairwise tables — sensitivity

Restricted to the 190 rDR-positive examinations.

Cell order used in the Python code:

    (both correct,
     A correct / B incorrect,
     A incorrect / B correct,
     both incorrect)

| Comparison | Both correct | A only | B only | Both incorrect |
|---|---:|---:|---:|---:|
| OphtAI vs RETFound | 140 | 41 | 4 | 5 |
| OphtAI vs ORDER-DR | 116 | 65 | 0 | 9 |
| RETFound vs ORDER-DR | 100 | 44 | 16 | 30 |

---

## 5. Pairwise tables — specificity

Restricted to the 684 rDR-negative examinations.

Here, "correct" means predicted negative.

Cell order used in the Python code:

    (both correct,
     A correct / B incorrect,
     A incorrect / B correct,
     both incorrect)

| Comparison | Both correct | A only | B only | Both incorrect |
|---|---:|---:|---:|---:|
| OphtAI vs RETFound | 514 | 116 | 33 | 21 |
| OphtAI vs ORDER-DR | 620 | 10 | 47 | 7 |
| RETFound vs ORDER-DR | 542 | 5 | 125 | 12 |

---

## 6. Joint tables — PPV and NPV comparisons

Cell order used in the Python code:

    (d11, d10, d01, d00,
     n11, n10, n01, n00)

where:

- `d` = rDR-positive examination
- `n` = rDR-negative examination
- first binary index = prediction of system A
- second binary index = prediction of system B
- `1` = predicted positive
- `0` = predicted negative

| Comparison | d11 | d10 | d01 | d00 | n11 | n10 | n01 | n00 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| OphtAI vs RETFound | 140 | 41 | 4 | 5 | 21 | 33 | 116 | 514 |
| OphtAI vs ORDER-DR | 116 | 65 | 0 | 9 | 7 | 47 | 10 | 620 |
| RETFound vs ORDER-DR | 100 | 44 | 16 | 30 | 12 | 125 | 5 | 542 |

---

## 7. Figure 1 — sensitivity/specificity uncertainty

| Parameter | Value |
|---|---:|
| Significance level (`ALPHA`) | 0.05 |
| Simulations per sample size (`N_SIMULATIONS`) | 10,000 |
| Random seed (`RANDOM_STATE`) | 42 |
| Minimum effective sample size | 1 |
| Maximum effective sample size | 100 |
| True performance levels | 0.95, 0.90, 0.85, 0.80 |

Output:

    fig1.pdf

---

## 8. Figure 2 — PPV/NPV uncertainty

| Parameter | Value |
|---|---:|
| Significance level (`ALPHA`) | 0.05 |
| Minimum total sample size | 1 |
| Maximum total sample size | 400 |

Outputs:

    fig2a.pdf  PPV
    fig2b.pdf  NPV

---

## 9. Figure 3 — sensitivity/specificity power

| Parameter | Value |
|---|---:|
| Significance level (`ALPHA`) | 0.05 |
| Simulations per sample size (`N_SIMULATIONS`) | 10,000 |
| Random seed (`RANDOM_STATE`) | 42 |
| Minimum class-specific sample size | 1 |
| Maximum class-specific sample size | 300 |
| Reported target power | 0.80 |

Outputs:

    fig3a.pdf  Sensitivity
    fig3b.pdf  Specificity

---

## 10. Figure 4 — PPV/NPV power

| Parameter | Value |
|---|---:|
| Significance level (`ALPHA`) | 0.05 |
| Simulations per sample size (`N_SIMULATIONS`) | 10,000 |
| Random seed (`RANDOM_STATE`) | 42 |
| Minimum total sample size | 1 |
| Maximum simulated total sample size | 2,000 |
| PPV display range | 0–800 |
| NPV display range | 0–1,600 |
| Reported target power | 0.80 |

Outputs:

    fig4a.pdf  PPV
    fig4b.pdf  NPV

---

## 11. Common graphical parameters

| Parameter | Value |
|---|---:|
| Figure size | 5 × 5 in |
| Reference sample size | 10 |
| Y-axis range | 0–1 |
| Y major tick interval | 0.10 |
| Y minor tick interval | 0.05 |
| Power reference line | 0.80 |

