# FlexibleJM
Flexible joint model for time-to-event and non-Gaussian longitudinal outcomes
This repository contains a **light modification of the R package JMBayes2**
to allow **flexible (non-linear) effects of baseline covariates in the survival submodel**
of joint models for longitudinal and time-to-event data.

The goal is to minimally extend the original JMBayes2 framework by adding
**penalized B-spline effects** for up to two continuous covariates in the hazard function,
while keeping the original model structure and C++ backend unchanged.

---

## Model extension

Starting from the standard joint model formulation in JMBayes2, the hazard function is extended as:

\[
h_i(t) = h_0(t)\,\exp\left\{
\boldsymbol{\gamma}^\top \mathbf{w}_i
+ \sum_{j=1}^{J} f_j(\nu_{ij})
+ \alpha\,\mu_i(t)
\right\},
\quad J \in \{0,1,2\}
\]

where:
- \(\nu_{ij}\) are continuous baseline covariates,
- \(f_j(\cdot)\) are smooth functions modeled using **penalized B-splines**,
- \(\mu_i(t)\) is the longitudinal predictor as in standard joint models.

Each smooth term is approximated as:
\[
f_j(\nu_{ij}) = \sum_{k=1}^{K_j} \gamma_{\nu,jk} B_{jk}(\nu_{ij})
\]

with a second-order difference penalty on the spline coefficients.

---

## Scope and limitations

- This code is **not a standalone R package**.
- It is intended as a **research prototype** extending JMBayes2.
- Only **0, 1 or 2 flexible effects** are supported.
- The implementation focuses on **clarity and minimal code changes**, not generality.

---

## Data convention for the survival model

For simplicity, the survival design matrix (`dataS_h`, `dataS_H`) is assumed to follow this structure:

| Column order | Content |
|-------------|---------|
| 1 | Subject ID |
| 2..p | Linear baseline covariates |
| last columns | Flexible baseline covariates (1 or 2) |

Flexible covariates must be placed **at the end** of the survival dataset.

Boundary knots for the spline basis are automatically set to the **observed min/max**
of each flexible covariate.

---

## Example usage

```r
library(survival)
library(nlme)

# Survival submodel
survObject <- coxph(
  Surv(eventtime, status) ~ Zf1 + Zf2 + Zc + Zh,
  data = d$Event,
  x = TRUE
)

# Longitudinal submodel
lmeObject <- lme(
  Yij_1 ~ tij + Znl + Zl,
  random = ~ tij | id,
  data = d$Long1
)

# Joint model with flexible effects
jmFit <- jm(
  Surv_object = survObject,
  Mixed_objects = lmeObject,
  time_var = "tij",
  nterm1 = 8,   # number of spline coefficients for first flexible covariate
  nterm2 = 8,   # number of spline coefficients for second flexible covariate
  short.output = TRUE,
  control = list(
    n_chains = 1,
    n_iter = 20500,
    n_burnin = 500,
    base_hazard_segments = 4,
    seed = 1
  )
)
