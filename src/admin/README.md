# Administrative Programs
This directory contains administrative programs for the project managing tasks that do not fit in any steps of methods development or data analysis. These include sync resources, publishing data, etc. This folder is meant for use by developers only, not for users.

To upload contents in `data` and `output` to server (for archive and backup),

```
snakemake upload_assets
```

To download assets,

```
snakemake download_assets
```
