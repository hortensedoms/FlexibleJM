# FlexibleJM
This repository contains the code associated with the article:

> Doms H, Lambert P, Legrand C (2024).  
> *Flexible joint model for time-to-event and non-Gaussian longitudinal outcomes*.  
> Statistical Methods in Medical Research.

The FlexibleJM approach is a **light modification of the R package JMBayes2**
designed to allow **flexible (non-linear) effects of baseline covariates in the survival submodel**
of joint models for longitudinal and time-to-event data.

The archive `JMbayes2_0.4-5_Flexible.tar.gz` is provided for convenience.
It contains the research code used in the associated article.
This is not yet a formal R package.


## Data requirements for the survival model

The survival dataset (`data_event`) must contain **at least** the following variables,
provided in the **exact order**:

1. subject identifier (`id`),
2. flexible baseline covariates (if any),
3. linear baseline covariates,
4. event time,
5. event indicator.

Flexible covariates must appear before linear covariates and be placed consecutively.
This ordering is assumed internally when constructing the survival design matrices.

## Example usage

```r
library(survival)
library(nlme)

# Survival submodel
survObject <- coxph(Surv(eventtime, status) ~ Wf1 + Wf2 + W1 + W2, data = data_event,x = TRUE)

# Longitudinal submodel
lmeObject <- lme(Yij_1 ~ tij + X1 + X2, random = ~ tij | id,data = data_long)

# Joint model with flexible effects
jmFit <- jm(
  Surv_object = survObject,
  Mixed_objects = lmeObject,
  time_var = "tij",
  nterm1 = 8,   # number of spline coefficients for first flexible covariate Wf1
  nterm2 = 8,   # number of spline coefficients for second flexible covariate Wf2
  short.output = TRUE,
  control = list(
    n_chains = 1,
    n_iter = 10500,
    n_burnin = 500,
    base_hazard_segments = 4,
    seed = 1
  )
)
```
## Acknowledgements

This work is based on the R package **JMbayes2**, developed by Dimitris Rizopoulos
and collaborators:

https://github.com/drizopoulos/JMbayes2

All credit for the original implementation goes to the authors of JMBayes2.
This repository only presents an extension for research purposes.
