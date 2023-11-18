#!/bin/bash

SAMPLE=$1
DATA=$2
DTYPE=$3
GBUILD=$4
NorT=$5
script_path=$6

cmd="python $script_path/load_high20.py $SAMPLE $DATA $DTYPE $NorT"
echo -e "******************** load high_20 ********************\n$cmd"
$cmd

cmd="bedtools intersect -a ${SAMPLE}.${DTYPE}.snpRAW.txt -b $script_path/../ENCODE_blacklist/${GBUILD}-blacklist.v2.bed.nochr -v | cut -f 1,3- > ${SAMPLE}.${DTYPE}.snp.txt"
echo -e "******************** clean high_20 ********************\n$cmd"

#$cmd did not work
bedtools intersect -a ${SAMPLE}.${DTYPE}.snpRAW.txt -b $script_path/../ENCODE_blacklist/${GBUILD}-blacklist.v2.bed.nochr -v | cut -f 1,3- > ${SAMPLE}.${DTYPE}.snp.txt

rm ${SAMPLE}.${DTYPE}.snpRAW.txt


