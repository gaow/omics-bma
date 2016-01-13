# External Libraries
Dependencies of OmicsBMA are:

*  The `omics-bma` branch of [eqtlbma](https://github.com/gaow/eqtlbma.git)
*  Version 1.16 of [GNU GSL](http://www.gnu.org/software/gsl/). GSL 2.x will not work.
*  [DeepDish](https://github.com/gaow/deepdish)
*  Version 3.54 of [snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home)

These libraries can be downloaded via executing command `snakemake all` under the `src/external` folder.

## Installation Notes
*  GSL is built into the pyeqtlbma library by static linking. The `libgslcblas` cannot be statically linked unless the build is configured explicitly by `-fPIC`, e.g., `./configure CFLAGS='-O3 -fPIC'`. This has already been taken care of by the installation script for pyeqtlbma library.
