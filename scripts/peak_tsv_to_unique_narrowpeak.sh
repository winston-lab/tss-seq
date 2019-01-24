#!/bin/bash

tail -n + 2 | \
    sort -k1,1 -k2,2n -u | \
    cut -f1-10

