# Comparison of HM Fitting Methods
## Problem
For a toy data of 63,831 gene-snp pairs (about 1/250 the volume of total data) and 515 mixture components (Bayes Factors, about 1/2 the total number of components we want to examine), the hierarchical model was fitted using interior point method as implemented in OmicsBMA, and compared against the EM algorithm implemented in `eqtlbma_hm` and `ASHR`, using the same gene-snp pairs but 54 components (light toy).

## Implementation
To run the benchmark under `mix-model-benchmark`:

```
./performance.sos # light toy
./performance.sos --dataset data/test_association_l10abfs # toy
```

which will generate the test data, run both methods and report resource usage (via SoS, including time elapsed, CPU and RAM usage).

## Results
|  Method  |  CPU cores  |  RAM usage  |  Runtime  |
|:--------:|:-----------:|:-----------:|:---------:  |
|  EM eqtlbma on light toy  | 753.4% peak (575.1% avg)  | 1270.1 Mb peak (1189.4 Mb avg)  |  40.0 sec  |
|  EM ASHR on light toy  | | | |
|  IP on light toy | 796.4% peak (525.7% avg)  | 861.0 Mb peak (771.7 Mb avg)  | 10.2 sec  |
|  EM eqtlbma on toy  | | | |
|  EM ASHR on toy  | | | |
|  IP on toy  | | | |
|  IP on full (guestimate)  |  8  |  800GB  |  4 hours  |

where performance on the full data set is a guestimate assuming linear scale on both \(N\) and \(P\), to be conservative.

**On the light toy the relative RAM is 64% and runtime is 24%, for IP compared to EM.**

Note that if the input data with `float32` it should half the RAM usage. That will be 400GB RAM (at most) on all gene-snp pairs and >1,000 components. The results shows the potential to run the mixture model on large scale.
