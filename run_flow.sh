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
CLUSTERED_SEQUENCES=$GITDIR'/clustering/cath-domain-seqs-0.2_rep_seq.fasta' #Clustered sequences
RESULTS_DIR=/home/p/pbryant/pfs/results/CATH/20191017 #RESULTS DIRECTORY
$GITDIR/pdb_filter.py $CLUSTERED_SEQUENCES $GITDIR/Xray_res2.6Å.txt $RESULTS_DIR/
wait
#Select H-groups with at least two sequences and group them according to their H-group
HGROUPS=$GITDIR'/H_groups.tsv'
SEQUENCES=$GITDIR'/clustering/cath-domain-seqs-0.2_rep_seq.fasta'
FAILED_PDB_FILTER=$RESULTS_DIR'/failed_pdb_filter'
FASTADIR=/home/p/pbryant/pfs/data/CATH/complete_flow/distant #where to write the grouped fasta sequnces
singularity exec /pfs/nobackup/home/p/pbryant/singularity/bio.sif $GITDIR/above2.py --H_groups $HGROUPS --sequences $SEQUENCES --failed_pdb_filter $FAILED_PDB_FILTER --outdir $FASTADIR/


#Get data, make alignments and calculate scores
#best done in parallel
wait
#$GITDIR/get_parallel.sh

#wait
#Get all results into a unified df
#$GITDIR/merge_results.py
