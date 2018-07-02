
# TSS-seq analysis pipeline

## Requirements:

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)

## Instructions for use
1. Clone this repository.

```bash
git clone
```

2. Create and activate the (`snakemake_default`) virtual environment for the pipeline using conda. The virtual environment creation can take a while. If you've already created the `snakemake_default` environment from another one of my pipelines, this is the same environment, so you can skip creating the environment and just activate it.

```bash
# navigate into the pipeline directory
cd tss-seq

# create the snakemake_default environment
conda env create -v -f envs/default.yaml

# activate the environment
source activate snakemake_default

# to deactivate the environment
# source deactivate
```

3. Edit the configuration file template `config_template.yaml` and rename it to `config.yaml`.

```bash
vim config_template.yaml    # or use your favorite editor

# make your edits

# rename the configuration file
mv config_template.yaml config.yaml
```

4. Do a dry run to see what files will be created.

```bash
snakemake -p -R `snakemake --lc --li --lp` --rerun-incomplete --use-conda --dry-run
```

5. If running the pipeline on a local machine, you can run the pipeline using the above command, omitting the `--dry-run` flag. You can also use N cores by specifying the `--cores [N]` flag. The pipeline can also be run on an HPC cluster. On the HMS O2 cluster, which uses the SLURM job scheduler, entering `sbatch slurm.sh` will cause the pipeline to run as a single job which spawns individual subjobs as necessary. This can be adapted to other job schedulers and clusters by modifying `slurm.sh` and `cluster.yaml`, which specifies the resource requests for each type of job.

