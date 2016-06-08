# About
**OmicsBMA** is a machinery to interrogate genomics / transcriptomics data via Bayesian Model Averaging. The [github repository](https://github.com/gaow/omics-bma) hosts its [source code](https://github.com/gaow/omics-bma/tree/master) and [applications to data analysis](https://github.com/gaow/omics-bma/tree/analysis). Users should refer to [this wiki](https://github.com/gaow/omics-bma/wiki) for software related information and the [OmicsBMA labnotes](http://bioinformatics.org/labnotes/bma) for information on statistical model / methods description and data analysis applications.

# Concepts
## Notations
| Notation   | Definition    |  Notes |
|----------|:-------------:|------:|
| $J$ | Number of response variables | |
| $N$ | Sample size | |
| $Y\mapsto{\rm R}^{J}$ | | |
|$X$| | |
| $\beta\mapsto{\rm R}^{J}$ |  | |
|$\Sigma\mapsto{\rm R}^{J \times J}$ | | |
|$V$| \hat{\Sigma} | | |
| | | |
| | | |
| | | |

## Basic association model
Consider a multivariate regression problem
\[Y \mid X, \beta, \Sigma \sim N_J(X\beta, \Sigma)\]
The goal is to make inference on $\beta$. 
A Bayesian model for $\beta$ is adopted
\[\beta \mid U \sim N_J(0, U)\]

```
A Bayesian hierachical model with a spike-slap prior on $\beta$ is adopted
\[\beta \mid \bar{\beta}, U_\phi, \pi_0 \sim \pi_0\delta_0 + (1 - \pi_0)N_J(\bar{\beta}, U_\phi)\]
\[\bar{\beta} \mid U_\omega \sim N_J(0, U_\omega)\]
Therefore, \[\beta \mid U_\phi, U_\omega, \pi_0 \sim \pi_0\delta_o + (1 - \pi_0) N_J(0, U_\phi + U_\omega)\]
Let $U = U_\phi + U_\omega$, then \[\beta \mid U, \pi_0 \sim \pi_0\delta_0 + (1 - \pi_0)N_J(0, U)\]
```