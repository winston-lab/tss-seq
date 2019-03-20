#!/bin/bash

awk 'BEGIN{FS=OFS="\t"}{print $1, $2+$15, $2+$15+1, $4, $5, $6, $0}'
