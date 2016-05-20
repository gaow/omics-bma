## Interior Point vs. EM
This test compares MOSEK Implementation of Convex Optimization via Interior Point Method to EM algorithm. It simulates data and fit mixture model using both EM algorithm from `ashr` and a convex optimization approach implemented in MOSEK via Python interface. By default 10 replicates are evaluated in the test. For each test it starts with default 5 decimal places and check for agreement of up to 3 decimal places, then report the level of convergence at which the test pass. Any discrepancy more than 3 decimal places between the two methods will result in failure of the test.

To run the test, type:

```
./kwdual.sos compare
```
Most of the time the test will pass. For cases when the test does not pass I notice that the difference is usually at 2 decimal places and the results still roughly agree.

To generate one round and keep the data:
```
./kwdual.sos
```
