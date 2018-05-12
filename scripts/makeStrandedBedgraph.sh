#!/bin/bash

awk 'BEGIN{FS=OFS="\t"}{$1=$1"-plus; print $0}' $1 | cat - <(awk 'BEGIN{FS=OFS="\t"}{$1=$1"-minus"; print $0}' $2) | LC_COLLATE=C sort -k1,1 -k2,2n
