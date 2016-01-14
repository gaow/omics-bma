# Test MOSEK Implementation of Convex Optimization via Interior Point Method
This unit-test simulates data and fit mixture model using both EM algorithm and a convex optimization approach implemented in MOSEK via Python interface. 5 replicates are evaluated in the test and any discrepancy more than 3 decimal places between the two methods will result in failure of the unit test.

To run the test, type:

```
snakemake all
```

Most of the time the test will pass. For cases when the test does not pass I notice that the difference is usually at 2 decimal places and the results still roughly agree.
