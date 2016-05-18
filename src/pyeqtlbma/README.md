# Python Interface for Codes from Eqtlbma
## Overview
### Swig interface
`pyeqtlbma.i` configures the python-cpp binding, for 3 main work-horses: Bayes factor calculation, mixture model fitting and calculation of posterior quantities. It is written in `swig`. The interface functions are defined in `ibeqtlbma.hpp/.cpp` files.

### Data structure
Data from python are organized in homogeneous dictionaries (dictionaries of data of the same type) which will be converted to `std::map<std::string, TYPE>` in C++ where `TYPE` are `std::vector` implementation of vectors and matrices. 5 types are defined here:

*  Dictionary of 1D/2D float/int (4 types of data)
*  Dictionary of 1D string

### Bayes Factor Calculation
This is to essentially rewrite `eqtlbma_bf.cpp` to `pyeqtlbma::eqtlbma_bf()`. Some extension have been made to

*  Allow for directly working with prior matrices instead of providing rules (grids + configuration names as in the original `eqtlbma` program) to construct them.
*  Allow for taking user input \(V\) matrix when input data is summary statistics
 *  Appropriate \(V\) should be computed from expression data
 *  Diagonal matrix is used when such calculation is not possible
*  Allow for HDF5 format of input summary statistics
*  All output are in HDF5 format

## BF()
Here are some notes on `pyeqtlbma::eqtlbma_bf()`

### Input
Input parameters are loosely defined by using homogeneous dictionaries that contains information for the `run()` function in the `eqtlbma_bf.cpp` file. As a result, naming conventions of parameters follows from the input of the `run()` function.
