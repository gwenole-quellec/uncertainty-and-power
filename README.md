# uncertainty-and-power

Code used to simulate uncertainty in sensitivity estimates and statistical power in small-sample evaluations of AI screening systems.

This repository reproduces the figures from the manuscript:

**When ten cases are not enough: Uncertainty and statistical power in AI screening evaluations**

## What the code does

The script `uncertainty_and_power.py` generates two figures:

1. **Uncertainty on sensitivity estimates**  
   Mean width of the exact 95% confidence interval (Clopper–Pearson) for sensitivity as a function of the number of positive cases.

2. **Statistical power to detect differences in sensitivity**  
   Estimated power of McNemar’s test for paired comparisons between models with different sensitivities.

## Requirements

Install dependencies with:

```bash
pip install -r requirements.txt
