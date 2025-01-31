#!/bin/bash

#################################################################################
# Default setting
#################################################################################

SAMPLE='NA'
BAM='NA'
DTYPE='NA'
PANELBIN='NA'
OUTPUT=$(pwd)

#################################################################################
# Display usage
#################################################################################

usage()
{
    echo ""
    echo "****************************************************************************************************"
    echo "Program: generate the coverage bedGraph file for SubChrom"
    echo "****************************************************************************************************"
    echo ""
    echo "Usage: $0 -s sample_ID -i sample.bam -d WES -p panel_bin.bed"
    echo ""
    echo -e "   -h  --help                Show this help message and exit"
    echo ""
    echo -e "   -s  --sample             Unique sample ID"
    echo -e "   -i  --input              /path/to/sample.bam"
    echo -e "   -d  --data_type          Sequencing data type. Options: WGS, WES, cfDNA, panel, custom, etc."
    echo -e "   -p  --panel_bin          Bed file of your panel bins. /path/to/panel_bin.bed"
    echo -e "   -o  --output             Ouput directory. Default: . (current directory)"
    echo ""
}

#################################################################################
# Read inputs
#################################################################################

while [ $# -gt 0 ]; do
    case "$1" in
        -h | --help)
            usage
            exit 0
            ;;
        -s | --sample)
            SAMPLE=$2
            shift
            ;;
        -i | --input)
            BAM=$(realpath $2)
            shift
            ;;
        -d | --data_type)
            DTYPE=$2
            shift
            ;;
        -p | --panel_bin)
            PANELBIN=$2
            shift
            ;;
        -o | --output)
            OUTPUT=$(realpath $2)
            shift
            ;;
        *)
            echo "ERROR: unknown option \"$1\""
            usage
            exit 1
            ;;
    esac
    shift
done

#################################################################################
# Setup
#################################################################################

if [[ $SAMPLE == 'NA' ]]; then
    echo "Need to specify: -s --sample"
    usage
    exit 1
elif [[ $BAM == 'NA' ]]; then
    echo "Need to specify: -i --input"
    usage
    exit 1    
elif [[ ! -f $BAM ]]; then
    echo "ERROR: Input bam file ($BAM) not found!"
    usage
    exit 1
elif [[ $PANELBIN != 'WGS' ]] && [[ ! -f $PANELBIN ]]; then
    echo "ERROR: panel bin bed file ($PANELBIN) not found!"
    usage
    exit 1
fi

#################################################################################
# Analysis
#################################################################################

if [[ ! -f ${SAMPLE}.${DTYPE}.bw ]]; then
    bamCoverage -b $BAM -o ${SAMPLE}.${DTYPE}.bw
fi

if [[ $DTYPE == 'WGS' ]]; then
    echo 0 > ${SAMPLE}.${DTYPE}.bedGraph
else
    multiBigwigSummary BED-file -b ${SAMPLE}.${DTYPE}.bw --BED $PANELBIN --outRawCounts ${SAMPLE}.${DTYPE}.bedGraph -o deeptools.npz
    rm deeptools.npz
fi
