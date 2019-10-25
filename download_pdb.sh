#!/bin/bash -l

#A script for downloading .pdb files from cath using wget

#File with uids to download
SELECTED_UIDS=selected_uids.txt
#CATH api
CATH_API=www.cathdb.info/version/v4_2_0/api/rest/id/
while read p
do
  wget $CATH_API$p'.pdb'
done <$SELECTED_UIDS
