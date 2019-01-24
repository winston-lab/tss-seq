#!/bin/bash

tail -n +2 | \
    cut -f1-7,10,11,15 | \
    sort -k1,1 -k2,2n -u | \
    sort -k9,9nr


