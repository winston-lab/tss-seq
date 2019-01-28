#!/bin/bash

tail -n +2 | \
    cut -f1-7,10,11,15 | \
    sort -k1,1 -k2,2n -u | \
    awk 'BEGIN{FS=OFS="\t"}{if($5=="NA"){$5=0}; if($9=="NA"){$9=0}; print $0}' | \
    sort -k9,9nr


