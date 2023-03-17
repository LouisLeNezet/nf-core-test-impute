#!/bin/sh
REGION=38
MAP_SNP=./data/illumina/chr38.snp.map

# Select only chr 21 in the region of interest, filter out only SNP and concatenate it as chrX:X
cat /groups/dog/data/canFam4/liftover/CanFam3.1_GDF_1.0.sites.csv | \
    awk -F',' '$2 == "\"chr38\"" { print $2":"$3}' | \
    sed 's/\"//g' | \
    sed 's/chr//g' \
    > $MAP_SNP