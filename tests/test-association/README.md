## Compare `eqtlbma_bf` and `OmicsBMA::test_association`
This test depends on the `library-interface` test. `library-interface` generated `output/test_association_*.h5` files from `OmicsBMA::test_association`. The goal of this test is to compute the same quantities using `eqtlbma_bf` and compare the results, at the precision level of 5 decimal places (default).
```
./test_association.sos
```
