# Comparison of HM Fitting Methods
## Problem
For a toy data of 63,831 gene-snp pairs (about 1/250 the volume of total data) and 515 mixture components (Bayes Factors, about 1/2 the total number of components we want to examine), the hierarchical model was fitted using interior point method as implemented in OmicsBMA, and compared against the EM algorithm implemented in `eqtlbma_hm`, using the same gene-snp pairs but 54 components.

## Implementation
To run the benchmark under `mix-model-benchmark`

```
snakemake all
```

which will generate the test data, run both methods and report summaries.

## Results
### Convex optimization on 515 components
The interior point method works on the toy data. Below are summaries and estimated performance on the full data set (assuming linear scale on both \(N\) and \(P\), to be conservative)

|  Method  |  CPU cores  |  RAM usage  |  Runtime  |
|:--------:|:-----------:|:-----------:|:---------:  |
|  IP on toy  |  8  |  1.6GB  |  63 sec  |
|  IP on full  |  8  |  800GB  |  4 hours  |

Note that if the input data with `float32` it should half the RAM usage. That will be 400GB RAM (at most) on all gene-snp pairs and >1,000 components. The results shows the potential to run the mixture model on large scale.

### Comparison between IP and EM
54 components are included in this comparison.

|  Method  |  CPU cores  |  RAM usage  |  Runtime  |  RAM relative  |  runtime relative  |
|:--------:|:-----------:|:-----------:|:---------:|:--------------:|:------------------:  |
|  EM  |  8  |  1.1GB  |  41 sec  |  -  |  -  |
|  IP  |  8  |  0.7GB  |  9.72 sec  |  64%  |  24%  |
