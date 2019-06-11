#!/bin/bash

awk 'BEGIN{FS=OFS="\t"}{print $1, $2+$10, $2+$10+1, $4, $5, $6, $0}'
