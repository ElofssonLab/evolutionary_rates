#!/bin/bash -l


#A bash script for getting CATH files for a certain H-group.
#Sequence and structural alignments are then performed.
#Different metrics are thereafter calculated: secondary structure,
#surface accessibility, RMSD, Maximum likelihood amino acid distance,
#TMscore, LDDT-score

#Filter the clustered sequences according to experimental method used for
#structural determination and X-ray resolution
CLUSTERED_SEQUENCES='cath-domain-seqs-0.2_rep_seq.fasta' #Clustered sequences
RESULTS_DIR='.' #RESULTS DIRECTORY
./pdb_filter.py $CLUSTERED_SEQUENCES Xray_res2.6Å.txt ./home/pbryant/data/CATH/mmseqs2/ $RESULTS_DIR/

#Select H-groups with at least two sequences and group them according to their H-group
HGROUPS='H_group.tsv'
SEQUENCES='clustering/cath-domain-seqs-0.2_rep_seq.fasta'
FAILED_PDB_FILTER='failed_pdb_filter_2.6Å.txt'
./above2.py $HGROUPS $SEQUENCES $FAILED_PDB_FILTER $RESULTS_DIR/




#The next step of aligning sequences and structures as well as getting
#pdb files is to be done in parallel for each H-group

#Align and get pdb files
$FILE_NAME="H-group"

#Go to results directory
cd $RESULTS_DIR
#Program input paths
CATH_API=www.cathdb.info/version/v4_2_0/api/rest/id/
HHBLITS=/hh-suite/build/bin/hhblits #Path to hhblits version 3.1.0
HHALIGN=/hh-suite/build/bin/hhalign #Path to hhalign version 3.1.0
UNIPROT=/uniprot20_2016_02/data/uniprot20_2016_02 #Path to uniprot. Here version 20 was used.

PUZZLE=/tree-puzzle-5.3.rc16-linux/src/puzzle #Path to tree-puzzle
TMALIGN=/TMalign #Path to TMalign Version 20170708
TMSCORE=/TMscore #Path to TMscore Version 2016/03/23
LDDT_IMAGE=/home/singularity/ost.img #VERSION 1.9.0
#Get ids and run hhalign
./get_data.py $RESULTS_DIR/$FILE_NAME/ $RESULTS_DIR/ $HHBLITS $HHALIGN $UNIPROT 15 $CATH_API

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
#Make TMalign dir
mkdir $RESULTS_DIR/TMalign/
wait
#Run TMalign and tree-puzzle
./run_tmalign_treepuzzle_ind.py $RESULTS_DIR/ $RESULTS_DIR/TMalign/ $RESULTS_DIR/$FILE_NAME/ $FILE_NAME $PUZZLE $TMALIGN

#Go to where the files are
cd $RESULTS_DIR/TMalign/
wait
#Run lddt for pdbs
./run_lddt.py $RESULTS_DIR/TMalign/ $RESULTS_DIR/TMalign/ guide $LDDT_IMAGE

wait
#Make TMScore dir
mkdir $RESULTS_DIR/TMscore/
cd $RESULTS_DIR
wait
#Run TMscore and tree-puzzle
./run_tmscore_treepuzzle.py $RESULTS_DIR/ $RESULTS_DIR/TMscore/ $RESULTS_DIR/$FILE_NAME/ $FILE_NAME $PUZZLE $TMSCORE

wait
#Run lddt for aligned pdbs
./run_lddt.py $RESULTS_DIR/ $RESULTS_DIR/TMscore/ guide $LDDT_IMAGE


#Get all results into a unified df
