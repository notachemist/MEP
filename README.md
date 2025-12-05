# Modular Electrochromic Polyester (MEP) Color Targeting

Python tools for predicting and targeting colors obtained from blends of modular electrochromic polyesters (MEPs) using the CIE 1976 L*a*b* (CIELAB) color space.


Core capabilities:

- Convert experimental transmission spectra of MEP devices into L*a*b* coordinates.
- Predict colors from mixtures of two MEPs using linear interpolation of spectra in CIELAB space.
- Predict colors from mixtures of three MEPs using barycentric coordinates in the a*–b* plane.
- Compute ΔE color differences between predicted mixtures and a user-defined target color.
- Suggest “best starting ratios” of MEPs for experimental color targeting.

---

## Background

MEPs are electrochromic polyesters built from aromatic ester macromonomers (e.g. terephthalate, naphthalate, thiophenate) linked by flexible aliphatic segments (e.g. adipate, hexamethylene). They combine:

- Bright, high-contrast colors in the reduced state (PHAT = magenta, PHAN = green, PHATh = cyan, etc.).
- High neutral-state transparency.
- Good processability and chemical recyclability via ester hydrolysis.

In the reduced state, blends of MEPs behave approximately linearly in CIELAB space: the color of a blend lies along (2 MEPs) or inside (3 MEPs) the convex hull spanned by the individual MEP colors. This makes them ideal for model-based color targeting instead of brute-force experimental screening.

All color calculations are performed in CIE 1976 L*a*b*:

- L* – lightness (0 = black, 100 = white)
- a* – green (−) to red (+)
- b* – blue (−) to yellow (+)

Color differences are quantified using the CIE76 metric:

Delta E = sqrt{(L_2^* - L_1^*)^2 + (a_2^* - a_1^*)^2 + (b_2^* - b_1^*)^2}

---

## Methods Overview

### 2-MEP Targeting (Linear Interpolation)

For two MEPs with spectra A_1(λ) and A_2(λ), the code constructs a family of theoretical spectra:

A_th(λ, k_i) = α_i A1(λ) + (1-α_i) A2(λ)

with α_i in [0,1] sampled at user-defined resolution (e.g. 5% or 10% steps).

Each spectrum is:

1. Converted to L*a*b* (given an illuminant and observer).
2. Compared to the target L*a*b* color via ΔE.
3. Ranked to return the MEP mixing ratio that minimizes ΔE.


### 3-MEP Targeting (Barycentric Coordinates in a*–b*)
For three MEPs (for example PHAT, PHAN, PHATh) with reduced-state coordinates
A = (a_A*, b_A*), B = (a_B*, b_B*), C = (a_C*, b_C*)
and a target color
P = (a_P*, b_P*),
the code solves for barycentric weights lambda_A, lambda_B, lambda_C such that

P ≈ lambda_A A + lambda_B B + lambda_C C

with

lambda_A + lambda_B + lambda_C = 1 and lambda_i >= 0 for all i.

The weights directly correspond to the blend ratio of each MEP. In the current implementation, the targeting is carried out in the a*–b* plane (hue). The same framework naturally extends to four MEPs in full 3D Lab* space (tetrahedral barycentric coordinates) if one wishes to incorporate the L* axis in the targeting.

---

## Repository Structure


MEP/
│
├─ 2MEP targeting/
│   ├─ PHAN.csv
│   ├─ PHAT.csv
│   ├─ PHATh.csv
│   ├─ color_coor_ref.csv
│   ├─ photopic_ref.csv
│   └─ theo.py
│
├─ 3MEP targeting/
│   ├─ bary-targeting.py
│   ├─ color_coor_ref.csv
│   ├─ photopic_ref.csv
│   └─ real.csv
│
├─ README.md
└─ LICENSE
