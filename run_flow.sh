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

#Select H-groups with at least two sequences and group them according to their H-group
HGROUPS='H_group.tsv'
SEQUENCES='clustering/cath-domain-seqs-0.2_rep_seq.fasta'
FAILED_PDB_FILTER='failed_pdb_filter_2.6Å.txt'
FASTADIR=/home/p/pbryant/pfs/data/CATH/below95_2.6Å_above2_grouped#where to write the grouped fasta sequnces
$GITDIR/above2.py --H_groups $HGROUPS --sequences $SEQUENCES --failed_pdb_filter $FAILED_PDB_FILTER --outdir $FASTADIR/




#The next step of aligning sequences and structures as well as getting
#pdb files is to be done in parallel for each H-group

move into get_parallel.sh
#Align and get pdb files
$FILE_NAME="H-group"

#Go to results directory
TEMPDIR='/scratch/'$FILE_NAME #Temporary dir to write the output to
cd $RESULTS_DIR
#Program input paths
CATH_API=www.cathdb.info/version/v4_2_0/api/rest/id/
HHBLITS=/home/p/pbryant/pfs/hh-suite/build/bin/hhblits #Path to hhblits version 3.1.0
HHALIGN=/home/p/pbryant/pfs/hh-suite/build/bin/hhalign #Path to hhalign version 3.1.0
UNIPROT=/home/p/pbryant/pfs/uniprot20_2016_02/data/uniprot20_2016_02 #Path to uniprot. Here version 20 was used.

PUZZLE=/home/p/pbryant/pfs/tree-puzzle-5.3.rc16-linux/src/puzzle #Path to tree-puzzle
TMALIGN=/home/p/pbryant/pfs/TMalign #Path to TMalign Version 20170708
TMSCORE=/home/p/pbryant/pfs/TMscore #Path to TMscore Version 2016/03/23
LDDT_IMAGE=/home/p/pbryant/pfs/singularity/ost.img #VERSION 1.9.0

#Get ids and run hhalign
./get_data.py --input_dir FASTADIR/$FILE_NAME/ --output_dir $RESULTS_DIR/$FILE_NAME/ --hhblits $HHBLITS $HHALIGN $UNIPROT 15 $CATH_API

#Run dssp
DSSP=/home/p/pbryant/pfs/dssp
wait
mkdir $RESULTS_DIR/dssp
wait
#Run dssp for pdbs
for file in $RESULTS_DIR/*.pdb
do
if [[ $file != *aln* ]] && [[ $file != *rf_* ]];  #If no aln or rf_ in name
then
$DSSP $file > $file'.dssp'
fi
done
wait
#Move all files to dssp dir
mv $RESULTS_DIR/*.dssp $RESULTS_DIR/dssp/

wait
#Make structure alignment dir
mkdir $RESULTS_DIR/structure/
wait
#Run TMalign and tree-puzzle
./run_tmalign_treepuzzle_ind.py $RESULTS_DIR/ $RESULTS_DIR/structure/ $RESULTS_DIR/$FILE_NAME/ $FILE_NAME $PUZZLE $TMALIGN

#Go to where the files are
cd $RESULTS_DIR/structure/
wait
#Run lddt for pdbs
./run_lddt.py $RESULTS_DIR/structure/ $RESULTS_DIR/structure/ guide $LDDT_IMAGE

wait
#Make sequence alignment dir
mkdir $RESULTS_DIR/sequence/
#Move all phylip files from sequence alignment to sequnce directory
mv $RESULTS_DIR/*.phy $RESULTS_DIR/sequence
cd $RESULTS_DIR/sequence
wait
#Run TMscore and tree-puzzle
./run_tmscore_treepuzzle.py $RESULTS_DIR/ $RESULTS_DIR/sequence/ $RESULTS_DIR/$FILE_NAME/ $FILE_NAME $PUZZLE $TMSCORE

wait
#Run lddt for aligned pdbs
./run_lddt.py $RESULTS_DIR/ $RESULTS_DIR/sequence/ guide $LDDT_IMAGE


#Get all results into a unified df
./merge_results.py
