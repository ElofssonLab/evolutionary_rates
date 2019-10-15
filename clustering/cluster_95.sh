#!/bin/bash -l

#CATH v4_2_0
#Cluster sequences on  95%
#CD-HIT version 4.8.1
cdhit/cd-hit -i cath-domain-seqs-S95.fa -o reallybelow95.fa -c 0.95
