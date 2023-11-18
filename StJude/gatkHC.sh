#!/bin/bash

SAMPLE=$1
DATA=$2
DTYPE=$3
REF=$4
WCHR=$5

################ Split gatkHC runs ################
for chrN in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}; do
    if [[ $WCHR == 'N' ]]; then
        CHR=${chrN}
    elif [[ $WCHR == 'Y' ]]; then
        CHR=chr${chrN}
    fi

    # each chr
    cmd="gatk --java-options '-Xmx4G' HaplotypeCaller -R $REF -I $DATA -O chr${chrN}.vcf.gz -L $CHR"
    echo -e "******************** run gatkHC for chr${chrN} ********************\n$cmd\n"
    bsub -P SubChrom -J ${SAMPLE}.${chrN} -q standard -M 4000 -n 2 -R "span[hosts=1]" -eo chr${chrN}.err -oo chr${chrN}.out gatk --java-options "-Xmx4G" HaplotypeCaller -R $REF -I $DATA -O chr${chrN}.vcf.gz -L $CHR
done

# concat all chr to one file
cmd="bcftools concat -o ${SAMPLE}.${DTYPE}.gatkHC.vcf.gz chr1.vcf.gz chr2.vcf.gz chr3.vcf.gz chr4.vcf.gz chr5.vcf.gz chr6.vcf.gz chr7.vcf.gz chr8.vcf.gz chr9.vcf.gz chr10.vcf.gz chr11.vcf.gz chr12.vcf.gz chr13.vcf.gz chr14.vcf.gz chr15.vcf.gz chr16.vcf.gz chr17.vcf.gz chr18.vcf.gz chr19.vcf.gz chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz chrX.vcf.gz chrY.vcf.gz"
echo -e "******************** concat gatkHC calls ********************\n$cmd\n"
bsub -P SubChrom -J ${SAMPLE}.HCx -w "done(${SAMPLE}.1)&&done(${SAMPLE}.2)&&done(${SAMPLE}.3)&&done(${SAMPLE}.4)&&done(${SAMPLE}.5)&&done(${SAMPLE}.6)&&done(${SAMPLE}.7)&&done(${SAMPLE}.8)&&done(${SAMPLE}.9)&&done(${SAMPLE}.10)&&done(${SAMPLE}.11)&&done(${SAMPLE}.12)&&done(${SAMPLE}.13)&&done(${SAMPLE}.14)&&done(${SAMPLE}.15)&&done(${SAMPLE}.16)&&done(${SAMPLE}.17)&&done(${SAMPLE}.18)&&done(${SAMPLE}.19)&&done(${SAMPLE}.20)&&done(${SAMPLE}.21)&&done(${SAMPLE}.22)&&done(${SAMPLE}.X)&&done(${SAMPLE}.Y)" -q standard -M 3000 -eo HCx.err -oo HCx.out $cmd

# remove intermediate files
bsub -P SubChrom -J ${SAMPLE}.HCy -w "done(${SAMPLE}.HCx)" -q standard -M 3000 -eo HCy.err -oo HCy.out "rm chr*"
