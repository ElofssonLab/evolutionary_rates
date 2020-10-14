# evolutionary_rates
Workflow for comparing evolutionary rates of homologous protein domain pairs from CATH in terms of sequence and structure.

run_flow.sh
1. Cluster the data on 95 % sequence identity using CD-HIT. clustering/cluster_95.sh
2. Download all pdbids with X-ray crystallography and a resolution of less than 2.6 Å (Xray_res2.6Å.txt) (from pdb website)
3. Filter out the fasta sequences from the 95 % clustering based on the downloaded pdbids: 
./pdbfilter reallybelow95.fa Xray_res2.6Å.txt outputdir
Writes the fasta sequences that failed the filter into a file.
4. Select all H-groups with above two entries remaining after the filter and write to individual folders. 
./above2.py --H_groups H_groups.tsv --sequences clustering/reallybelow95.fa --failed_pdb_filter (out from previous) --outdir
5. get_parallel.sh
-Convert fasta files to HMMs using hhblits. All pairs of x randomly selected entries are then aligned with hhalign (/get_data.py)
-Download the selected pdb files. (download_pdb.sh)
-Run DSSP
-Run TMalign and tree-puzzle on pdb files (run_tmalign_treepuzzle_ind.py)
-Run lDDT on the TMalign alignmens (run_lddt.py)
-Run TMscore and tree-puzzle on structural information extracted according to the sequence alignments.
-Run lDDT on the structural information extracted according to the sequence alignments.
6. Merge all results
./merge_results.py 
  --indir INDIR         path to input directory.
  --threshold THRESHOLD Threshold for percent aligned of the shortest sequence in each pair.
  --outdir OUTDIR       path to output directory.
  --dssp_path DSSP_PATH path to dssp.
  --fastadir FASTADIR   path to fastadir.
  --gitdir GITDIR       path to gitdir.
  --logfile LOGFILE     logfile for stats.


