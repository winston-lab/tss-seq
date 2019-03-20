#!/bin/bash

awk 'BEGIN{FS=OFS="\t"}{start=$2+$10; print $1, start, start+1, $4, $5, $6}'
