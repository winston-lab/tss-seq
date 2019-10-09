#!/bin/bash

snakemake -pr \
    -R `cat <(snakemake --lc --rerun-incomplete) \
            <(snakemake --li --rerun-incomplete) \
            <(snakemake --lp --rerun-incomplete) | sort -u` \
    --cores \
    --rerun-incomplete \
    --use-conda
