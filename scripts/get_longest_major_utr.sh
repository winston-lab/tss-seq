#!/bin/bash

awk 'BEGIN{FS=OFS="\t"}
{if($3=="+" && $4<$10 && $5<$10) {
    end = $10;
    $6>0 ? start=$4 : start=$5;
    print $2, start, end, $1, 0, $3
    }
 if($3=="-" && $4>$11 && $5>$11) {
    start = $11;
    $6>0 ? end=$4+1 : end=$5+1;
    print $2, start, end, $1, 0, $3
    }
}'
