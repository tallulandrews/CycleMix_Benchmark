#!/bin/bash
# Run this file in bash with this command:  ./filename
HOST=ftp.ebi.ac.uk
USER=anonymous
ftp -pinv $HOST <<EOF
user $USER
cd biostudies/fire/E-MTAB-/805/E-MTAB-2805/Files
binary
mget "G2M_singlecells_counts.txt"
mget "G1_singlecells_counts.txt"
mget "S_singlecells_counts.txt"
disconnect
bye
EOF