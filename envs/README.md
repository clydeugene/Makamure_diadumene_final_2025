## Environment Specification Files

The files in this folder are used by Snakemake in running the workflow. We installed these conda environments directly on our server, because our snakemake implementation was having issues with creating new conda environments. So you may need to manually create the environments before running the pipeline.

> **Note:**
> Please ensure that you launch the pipeline from within the snakemake conda environment as it has both Snakemake and Apptainer.
