#!/bin/bash

awk 'BEGIN{FS=OFS="\t"}$6=="+"{$1=$1"-plus"; print $0} $6=="-"{$1=$1"-minus"; print $0}' $1
