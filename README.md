# 📊 uncertainty-and-power

Code used to simulate uncertainty in sensitivity estimates and statistical power in small-sample evaluations of AI screening systems.

This repository reproduces the figures from the manuscript:

> **When ten cases are not enough: Uncertainty and statistical power in AI screening evaluations**

📝 **Preprint available on Research Square:**  
[Research Square preprint](https://www.researchsquare.com/article/rs-9601448/v1)

---

## 🔍 What the code does

The script `uncertainty_and_power.py` generates two figures:

### 1. Uncertainty on sensitivity estimates 📉

Simulation of the mean width of the exact 95% confidence interval  
(Clopper–Pearson interval) for sensitivity as a function of the number of positive cases.

This illustrates how small evaluation cohorts can lead to highly unstable performance estimates.

### 2. Statistical power to detect differences in sensitivity ⚡

Simulation-based estimation of the statistical power of McNemar’s test for paired comparisons between AI systems with different sensitivities.

This highlights the difficulty of reliably detecting performance differences in low-sample screening studies.

---

## 📦 Requirements

Install dependencies with:

```bash
pip install -r requirements.txt
```

---

## ▶️ Usage

Run the simulation script with:

```bash
python uncertainty_and_power.py
```

The generated figures will be saved in the working directory.

---

## 🎯 Purpose

This repository accompanies the manuscript and aims to provide:
- reproducible simulations;
- intuitive visualizations of statistical uncertainty;
- practical insights for the evaluation of AI screening systems under limited sample sizes.

The code is intentionally lightweight and easy to adapt to other evaluation settings.

---

## 📚 Citation

If you use this code or build upon it, please cite the associated manuscript:

> Gwenolé Quellec. When ten cases are not enough: Uncertainty and power in AI screening evaluations, 07 May 2026, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-9601448/v1].
