# Python Interface for Codes from Eqtlbma
## Swig Interface
`pyeqtlbma.i` defines the interface, for 3 main work-horse: Bayes factor calculation, mixture model fitting and calculation of posterior quantities. Interface is written in `swig`.

### Data structure
Data from python are organized in homogeneous dictionaries (dictionaries of data of the same type) which will be converted to `std::map<std::string, TYPE>` in C++. 5 types are implemented here:

*  Dictionary of 1D/2D float/int (4 types of data)
*  Dictionary of 1D string

## Bayes Factor Calculation
This is to essentially rewrite `eqtlbma_bf.cpp` to `pyeqtlbma::bf()` as implemented in `bf.h` and `bf.cpp`. The original `eqtlbma` program was changed to

*  Allow for directly working with prior matrices rather than rules to construct them.
*  Allow for taking in \(V\) matrix when input data is summary statistics and use diagonal matrix as default.
