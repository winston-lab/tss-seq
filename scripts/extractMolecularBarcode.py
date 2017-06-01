#!/usr/bin/env python

"""
Date : April 24th 2014

script used to extract molecular barcodes from the reads and check their raw distribution before alignment

use : python extractMolecularBarcode.py inFastq outFastq outBarcodes outLigation

"""

import sys, itertools


iFastq=open(sys.argv[1], 'r')
oFastq=open(sys.argv[2], 'w')
oBarcode=open(sys.argv[3], 'w')
oLigation=open(sys.argv[4], 'w')

dicoBarcode={}
dicoLigation={}
nct='ACTGN'
for barcode in list(itertools.product(nct, repeat=6)):
    dicoBarcode["".join(barcode)] = 0
    dicoLigation["".join(barcode)] = 0


header= iFastq.readline().rstrip()
while header != '':
    totseq= iFastq.readline()
    plus = iFastq.readline()
    qual = iFastq.readline()
    barcode = totseq[0:6]
    ligation = totseq[3:9]
    seq = totseq[6:]
    oFastq.write(header.split(" ")[0]+'_MolecularBarcode:'+barcode+' '+header.split(" ")[1]+'\n')
    oFastq.write(seq)
    oFastq.write(plus)
    oFastq.write(qual[6:])
    header= iFastq.readline().rstrip()

    dicoBarcode[barcode] += 1
    if len(seq) >= 4 :
        dicoLigation[ligation] += 1

for barcode, times in dicoBarcode.items():
    oBarcode.write("%s\t%s\n" % (barcode, str(times)))

for ligation, times in dicoLigation.items():
    oLigation.write("%s\t%s\n" % (ligation, str(times)))

iFastq.close()
oFastq.close()
oBarcode.close()
oLigation.close()

