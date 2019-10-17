#!/bin/bash -l


#A bash script for getting CATH files for a certain H-group.
#Sequence and structural alignments are then performed.
#Different metrics are thereafter calculated: secondary structure,
#surface accessibility, RMSD, Maximum likelihood amino acid distance,
#TMscore, LDDT-score

#Path to git directory
GITDIR=/home/p/pbryant/pfs/evolutionary_rates
#Filter the clustered sequences according to experimental method used for
#structural determination and X-ray resolution
CLUSTERED_SEQUENCES='cath-domain-seqs-0.2_rep_seq.fasta' #Clustered sequences
RESULTS_DIR=/home/p/pbryant/pfs/results/CATH/20191017 #RESULTS DIRECTORY
$GITDIR/pdb_filter.py $CLUSTERED_SEQUENCES Xray_res2.6Å.txt ./home/pbryant/data/CATH/mmseqs2/ $RESULTS_DIR/
wait
#Select H-groups with at least two sequences and group them according to their H-group
HGROUPS='H_group.tsv'
SEQUENCES='clustering/cath-domain-seqs-0.2_rep_seq.fasta'
FAILED_PDB_FILTER='failed_pdb_filter_2.6Å.txt'
FASTADIR=/home/p/pbryant/pfs/data/CATH/below95_2.6Å_above2_grouped #where to write the grouped fasta sequnces
$GITDIR/above2.py --H_groups $HGROUPS --sequences $SEQUENCES --failed_pdb_filter $FAILED_PDB_FILTER --outdir $FASTADIR/


#Get data, make alignments and calculate scores
#best done in parallel
wait
$GITDIR/get_parallel.sh

wait
#Get all results into a unified df
$GITDIR/merge_results.py
