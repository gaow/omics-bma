# Source Codes
This folder contains source codes for the `OmicsBMA` [library](internal) as well as its executable workflows:

*  [generate-blocks](generate-blocks)
*  [make-prior](make-prior)
*  [test-for-association](test-for-association)
*  [consolidate-tests](consolidate-tests)
*  [fit-mixture-model](fit-mixture-model)
*  [bayesian-inference](bayesian-inference)
*  [visualization](visualization)

Quality control for coding is implemented as [unit-tests](unit-tests).

# TODO
*  Fix the temporary patch that links `__init__.py` to `pyeqtlbma.py` for package import.
*  Fix Issue 2 the bf bug.
*  Fix Issues 3 and 4 the exception handling.
*  Fix compression type and level for pandas interface via deepdish. Currently these options are ignored.
*  Allow append mode in deepdish.io.save().
*  Implement input parameter check with respect to the old omic-bma branch interface: when sumstats is used, error=ulvr, analys=sep and no permutation
*  Add column names to all BF() results
*  Implement customized priors for UVLR
