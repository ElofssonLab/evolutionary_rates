#!/bin/bash -l
#SBATCH -A SNIC2019-3-319
#SBATCH -c 1
#SBATCH -t 15:00:00
#SBATCH --array=1-1000
#SBATCH --error=/home/p/pbryant/pfs/results/CATH/20190904/error/%A_%a.error
#SBATCH --output=/home/p/pbryant/pfs/results/CATH/20190904/out/%A_%a.out

#The next step of aligning sequences and structures as well as getting
#pdb files is to be done in parallel for each H-group
#The array on the cluster used has a maximum of 1000 entries per job submission

#Path to git directory
GITDIR=/home/p/pbryant/pfs/evolutionary_rates
#Path to fasta files
FASTADIR=/home/p/pbryant/pfs/data/CATH/below95_2.6Å_above2_grouped
#Get filename from line index
HGROUPS=/home/p/pbryant/pfs/runs/CATH/20190904/over2_1.txt #This file has to contain the selected newline separated H-groups
FILE_NAME=$(sed -n $SLURM_ARRAY_TASK_ID'p' $HGROUPS)
echo $FILE_NAME

$FILE_NAME="H-group"

#Go to results directory
TEMPDIR='/scratch/'$FILE_NAME #Temporary dir to write the output to
mkdir $TEMPDIR
cd $TEMPDIR
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
$GITDIR/get_data.py --input_dir FASTADIR/$FILE_NAME/ --output_dir $TEMPDIR/ --hhblits $HHBLITS --hhalign $HHALIGN --uniprot $UNIPROT --get_max 15 --address $CATH_API

#Run dssp
DSSP=/home/p/pbryant/pfs/dssp
wait
mkdir $TEMPDIR/dssp
wait
#Run dssp for pdbs
for file in $TEMPDIR/*.pdb
do
if [[ $file != *aln* ]] && [[ $file != *rf_* ]];  #If no aln or rf_ in name
then
$DSSP $file > $file'.dssp'
fi
done
wait
#Move all files to dssp dir
mv $TEMPDIR/*.dssp $TEMPDIR/dssp/

wait
#Make structure alignment dir
mkdir $TEMPDIR/structure/
wait
#Run TMalign and tree-puzzle
$GITDIR/run_tmalign_treepuzzle_ind.py --indir $TEMPDIR/ --outdir $TEMPDIR/structure/ --fastadir $FASTADIR/$FILE_NAME/ --hgroup $FILE_NAME --puzzle $PUZZLE --TMalign $TMALIGN

#Go to where the files are
cd $TEMPDIR/structure/
wait
#Run lddt for pdbs
$GITDIR/run_lddt.py --indir $TEMPDIR/structure/ --outdir $TEMPDIR/structure/ --mode guide --lddt_image $LDDT_IMAGE

wait
#Make sequence alignment dir
mkdir $TEMPDIR/sequence/
#Move all phylip files from sequence alignment to sequnce directory
mv $TEMPDIR/*.phy $RESULTS_DIR/sequence
cd $RESULTS_DIR/sequence
wait
#Run TMscore and tree-puzzle
./run_tmscore_treepuzzle.py $RESULTS_DIR/ $RESULTS_DIR/sequence/ $RESULTS_DIR/$FILE_NAME/ $FILE_NAME $PUZZLE $TMSCORE

wait
#Run lddt for aligned pdbs
./run_lddt.py $RESULTS_DIR/ $RESULTS_DIR/sequence/ guide $LDDT_IMAGE
