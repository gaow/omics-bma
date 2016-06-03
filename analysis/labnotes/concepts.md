# About
**OmicsBMA** is a machinery to interrogate genomics / transcriptomics data via Bayesian Model Averaging. The [github repository](https://github.com/gaow/omics-bma) hosts its [source code](https://github.com/gaow/omics-bma/tree/master) and [applications to data analysis](https://github.com/gaow/omics-bma/tree/analysis). Users should refer to [this wiki](https://github.com/gaow/omics-bma/wiki) for software related information and the [OmicsBMA labnotes](http://bioinformatics.org/labnotes/bma) for information on statistical model / methods description and data analysis applications.

# Concepts
## Notations
| Notation   | Definition    |  Notes |
|----------|:-------------:|------:|
| $K$ | Number of response variables | |
| $N$ | Sample size | |
| $\mathbf{Y}\mapsto{\rm R}^{N \times K}$ | | |
|$\mathbf{X}\mapsto{\rm R}^N$| | |
| $\beta\mapsto{\rm R}^{K}$ |  | |
| | | |
| | | |
| | | |
| | | |

## Basic association model
Consider the multivariate regression problem \[\mathbf{Y} = \mathbf{X}\beta^T + E,\quad E \sim N_K(0, \Sigma)\]
