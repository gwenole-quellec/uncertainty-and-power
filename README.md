# 📊 uncertainty-and-power

Code used to quantify uncertainty and statistical power in evaluations of AI screening systems.

This repository reproduces the figures from the manuscript:

> **When ten cases are not enough: Uncertainty and statistical power in AI screening evaluations**

📝 **Preprint available on Research Square:**  
[Research Square preprint](https://www.researchsquare.com/article/rs-9601448)

---

## 🔍 What the code does

The repository contains four analysis scripts covering uncertainty and statistical power for sensitivity, specificity, positive predictive value (PPV), and negative predictive value (NPV).

### 1. Uncertainty in sensitivity and specificity 📉

`uncertainty-sens-spec.py`

Estimates the mean width of exact 95% Clopper–Pearson confidence intervals as a function of the metric-specific effective sample size.

### 2. Uncertainty in PPV and NPV 📉

`uncertainty-ppv-npv.py`

Computes the expected width of exact 95% Clopper–Pearson confidence intervals for PPV and NPV as a function of the total validation sample size, using the confusion matrices observed on Messidor-2.

### 3. Power for sensitivity and specificity comparisons ⚡

`power-sens-spec.py`

Estimates statistical power for paired sensitivity and specificity comparisons by simulation. Pairwise outcomes are generated from the empirical disagreement structure observed on Messidor-2 and compared using the exact two-sided McNemar test.

### 4. Power for PPV and NPV comparisons ⚡

`power-ppv-npv.py`

Estimates statistical power for paired PPV and NPV comparisons by simulation. The simulations preserve the empirical joint distribution of disease status and paired predictions observed on Messidor-2, and comparisons use the generalized score test for paired predictive values.

The shared plotting functions are defined in `plotting.py`.

Numerical values and analysis parameters used to generate the figures are documented in [`experimental-data.md`](experimental-data.md).

---

## 📦 Requirements

Install dependencies with:

```bash
pip install -r requirements.txt
```

---

## ▶️ Usage

Run the four analysis scripts with:

```bash
python uncertainty-sens-spec.py
python uncertainty-ppv-npv.py
python power-sens-spec.py
python power-ppv-npv.py
```

The generated figures are saved in the working directory:

- `fig1.pdf` — sensitivity/specificity uncertainty
- `fig2a.pdf`, `fig2b.pdf` — PPV/NPV uncertainty
- `fig3a.pdf`, `fig3b.pdf` — sensitivity/specificity power
- `fig4a.pdf`, `fig4b.pdf` — PPV/NPV power

---

## 🎯 Purpose

This repository accompanies the manuscript and aims to provide:

- reproducible analyses;
- transparent documentation of the experimental inputs and analysis parameters;
- intuitive visualizations of statistical uncertainty and power;
- practical insights for the evaluation of AI screening systems under limited sample sizes.

The code is intentionally lightweight and easy to adapt to other evaluation settings.

---

## 📚 Citation

If you use this code or build upon it, please cite the associated manuscript:

> Gwenolé Quellec. *When ten cases are not enough: Uncertainty and power in AI screening evaluations*, 07 May 2026, PREPRINT (Version 1), Research Square. https://doi.org/10.21203/rs.3.rs-9601448
