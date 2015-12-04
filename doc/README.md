# Documentation
Executing command `snakemake all` under this directory deploys documentation for each folder in this repository (the `README.md` file), as well the entire [OmicsBMA Wiki](http://bioinformatics.org/bma). `README.md` documents thus generated are the same material as on the [user manual](http://bioinformatics.org/bma/manual/start) and [developer's info](http://bioinformatics.org/bma/development/start) Wiki pages, only that the Wiki is better organized and easier to browse.

**Files in this directory are the only source files to build the entire documentation**. All contributions to documentation should go into this directory using syntax of `*.notes` file, along with appropriate edits to the `Snakemake` file.
