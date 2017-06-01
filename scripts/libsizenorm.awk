#!/bin/awk -f

## normalize two.bedgraph by million counts present in one.bedgraph
##
##  
## usage: libsizenorm.awk one.bedgraph two.bedgraph > two-normalized-by-one.bedgraph
##
BEGIN{
     FS=OFS="\t";
     sum=0
     }
{
    if(NR==FNR){
               sum+=$4
               }
    else{
        print $1, $2, $3, $4*1000000/sum
        }
}
