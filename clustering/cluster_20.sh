#!/bin/bash -l

#CATH v4_2_0
#Cluster sequences on 20 %
#MMseqs2 Version: f05f8c51d6e9c7c0b15fbd533e4b678303f50b3e
MMseqs2/build/bin/mmseqs easy-cluster cath-domain-seqs-S95.fa cath-domain-seqs-0.2 tmp --min-seq-id 0.2 -c 0.8 --cov-mode 1
