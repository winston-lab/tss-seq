#!/bin/bash

bsub -q priority -W 24:00 -n 4 snakemake -p --cluster-config cluster.yaml --use-conda --cluster "bsub -q {cluster.queue} -n {cluster.n} -W {cluster.time} -R 'rusage[mem={cluster.mem}]' -J {cluster.name} -e {cluster.err} -o {cluster.log}" --jobs 40


