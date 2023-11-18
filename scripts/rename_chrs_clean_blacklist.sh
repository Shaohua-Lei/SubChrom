#!/bin/bash

SAMPLE=$1
VCF=$2
GBUILD=$3
DTYPE=$4
script_path=$5

# rename
cmd="bcftools annotate --rename-chrs $script_path/../data/chr_name_conv.txt $VCF -Oz -o ${SAMPLE}.${DTYPE}.gatkHC_nochr.vcf.gz"
echo $cmd
$cmd

# index
bcftools index ${SAMPLE}.${DTYPE}.gatkHC_nochr.vcf.gz

# clean with ENCODE blacklist
cmd="bcftools view -v snps -T ^$script_path/../data/${GBUILD}-blacklist.v2.bed.nochr ${SAMPLE}.${DTYPE}.gatkHC_nochr.vcf.gz | bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AD{0}\t%AD{1}]\n' -o ${SAMPLE}.${DTYPE}.snp.txt"
echo $cmd

# $cmd did not work
bcftools view -v snps -T ^$script_path/../data/${GBUILD}-blacklist.v2.bed.nochr ${SAMPLE}.${DTYPE}.gatkHC_nochr.vcf.gz | bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AD{0}\t%AD{1}]\n' -o ${SAMPLE}.${DTYPE}.snp.txt
