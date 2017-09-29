#!/bin/bash

awk 'BEGIN{FS=OFS="\t"}$6=="+"{print $1"-plus", $2, $3, $4, $5, $6} $6=="-"{print $1"-minus", $2, $3, $4, $5, $6}' $1
