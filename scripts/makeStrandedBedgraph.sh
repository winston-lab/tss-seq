#!/bin/bash

awk 'BEGIN{FS=OFS="\t"}{print $1"-plus", $2, $3, $4}' $1 | cat - <(awk 'BEGIN{FS=OFS="\t"}{print $1"-minus", $2, $3, $4}' $2) | LC_COLLATE=C sort -k1,1 -k2,2n
