#!/bin/bash

#################################################################################
# Default setting
#################################################################################

SAMPLE='NA'
DATA='NA'
DTYPE='NA'
REF='NA'
PANELBIN='NA'
SNPmarker='NA'

OUTPUT=$(pwd)
GBUILD='hg38'
PON='none'
COVseg='False'
ROHseg='False'
minTF=0.1
minSize=1000000
minBins=100
GENDER='auto'
dipDep='auto'
covWindow='auto'
plotTF='True'
preFile='True'
GENES='cfDNA'

VERSION='0.1.0 (October 17, 2024)'

#################################################################################
# Display usage
#################################################################################

usage()
{
    echo ""
    echo "****************************************************************************************************"
    echo "Program: SubChrom CNV/cnLOH calling"
    echo -e "Version: $VERSION"
    echo "****************************************************************************************************"
    echo ""
    echo "Usage: $0 -s filename -i input.data -d cfDNA -r reference.fa -p panel_bin.bed -md hg38"
    echo ""
    echo -e "   -h  --help                Show this help message and exit"
    echo ""
    echo -e "Required arguments:"
    echo -e "   -s  --sample             Unique sample name for saving files"
    echo -e "   -i  --input              /path/to/input.data"
    echo -e "                               Format: BAM from WES and panel sequencing, BAM/vcf/high20 from WGS"
    echo -e "   -d  --data_type          Sequencing data type. Options: WGS, WES, cfDNA, panel, custom, etc."
    echo -e "   -r  --reference          Genome reference (same for bam files). /path/to/reference.fa"
    echo -e "   -p  --panel_bin          Bed file of your panel bins. /path/to/panel_bin.bed"
    echo -e "                               Options: 'WES' to use files from /SubChrom/data/, 'WGS' to skip"
    echo -e "   -md --marker_dir         Directory of SNP marker database from SubChrom. /path/to/SNPmarker"
    echo -e "                               Options: 'hg19' or 'hg38' for /SubChrom/data/SNPmarker_*/"
    echo ""
    echo -e "Optional arguments:"
    echo -e "   -o  --output             Ouput directory. Default: ./ (current directory)"
    echo -e "   -g  --genome_build       Options: hg38 (default), hg19"
    echo -e "   -n  --normal             Panel of Normals, or bedGraph file of a normal sample. /path/to/PoN.txt"
    echo -e "                               Default: none"
    echo -e "   -cs --coverage_seg       Perform coverage segmentation. Options: True, False (default)"
    echo -e "                               If --normal is none above for WES/panel data, this is False by default"
    echo -e "   -rs --ROH_seg            Perform ROH (runs of homozygosity) segmentation. Options: True, False (default)"
    echo -e "                               Required for high purity samples, in which no heterozygotes in cnLOH and loss"
    echo -e "   -mf --minTF              Minimal tumor fraction to report a CNV event. Default: 0.1. Minimum: 0.01"
    echo -e "   -ms --minSize            Minimal size (bp) to report a CNV event. Default: 1000000. Minimum: 10000"
    echo -e "   -mb --minBins            Minimal number of bins to report a CNV event. Default: 100. Minimum: 10"
    echo -e "   -sg --sample_gender      Options: Male/M, Female/F, auto. Default: auto for automatic detection"
    echo -e "   -dd --diploid_depth      How to compute the diploid depth."
    echo -e "                               Options: auto, chr1, ..., chr22, chrX, a specific value such as 500"
    echo -e "                               Default: auto for automatic optimization"
    echo -e "   -cw --covWinSize         Coverage window size (bp) for visualization"
    echo -e "                               Default: auto (WGS, 1000000; WES/cfDNA, 2000000. Minimum: 500000"
    echo -e "   -pf --plotTF             Plot tumor fraction or not. Options: True (default), False"
    echo -e "   -if --intermediate_file  Use intermediate files from the previous run, such as vcf file and segements"
    echo -e "                               Options: True (default), False (remove previous files and generate new ones)"
    echo -e "   -gl --gene_list          Gene list of interest for visualization. /path/to/geneList.bed"
    echo -e "                               Default: /SubChrom/data/geneList.bed"
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
            DATA=$(realpath $2)
            shift
            ;;
        -d | --data_type)
            DTYPE=$2
            shift
            ;;
        -r | --reference)
            REF=$2
            shift
            ;;
        -p | --panel_bin)
            PANELBIN=$2
            shift
            ;;
        -md | --marker_dir)
            SNPmarker=$2
            shift
            ;;
        -o | --output)
            OUTPUT=$(realpath $2)
            shift
            ;;
        -g | --genome_build)
            GBUILD=$2
            shift
            ;;  
        -n | --pon)
            PON=$2
            shift
            ;;
        -cs | --coverage_seg)
            COVseg=$2
            shift
            ;;
        -rs | --ROH_seg)
            ROHseg=$2
            shift
            ;;
        -mf | --minTF)
            minTF=$2
            shift
            ;;
        -ms | --minSize)
            minSize=$2
            shift
            ;;
        -mb | --minBins)
            minBins=$2
            shift
            ;;
        -sg | --sample_gender)
            GNEDER=$2
            shift
            ;;
        -dd | --diploid_depth)
            dipDep=$2
            shift
            ;;
        -cw | --covWinSize)
            covWindow=$2
            shift
            ;;
        -pf | --plotTF)
            plotTF=$2
            shift
            ;;
        -if | --intermediate_file)
            preFile=$2
            shift
            ;;
        -gl | --gene_list)
            GENES=$2
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

echo ""
# code & directory
script_path="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

########## check required arguments ##########
if [[ $SAMPLE == 'NA' ]]; then
    echo "Need to specify: -s --sample"
    usage
    exit 1
elif [[ $DATA == 'NA' ]]; then
    echo "Need to specify: -i --input"
    usage
    exit 1    
elif [[ ! -f $DATA ]]; then
    echo "ERROR: Input file ($DATA) not found!"
    usage
    exit 1
elif [[ $DTYPE == 'NA' ]]; then
    echo "Need to specify: -d --data_type"
    usage
    exit 1
elif [[ $REF == 'NA' ]]; then
    echo "Need to specify: -r --reference"
    usage
    exit 1
elif [[ $PANELBIN == 'NA' ]]; then
    echo "Need to specify: -p --panel_bin"
    usage
    exit 1
elif [[ $SNPmarker == 'NA' ]]; then
    echo "Need to specify: -md --marker_dir"
    usage
    exit 1
fi

# check reference & chr
if [ ! -f $REF ]; then
    echo "ERROR: genome reference ($REF) not found!"
    usage
    exit 1
else
    LINE=$(head -n1 $REF)
    if [[ $LINE == '>chr'* ]]; then
        G_CHR='Y'
    else
        G_CHR='N'
    fi
fi

# default panel bin
if [[ $PANELBIN == 'WES' ]] || [[ $PANELBIN == 'cfDNA' ]]; then 
    if [[ $G_CHR == 'Y' ]]; then
        PANELBIN=$script_path/../data/${GBUILD}.${PANELBIN}.bed.chr
    elif [[ $G_CHR == 'N' ]]; then
        PANELBIN=$script_path/../data/${GBUILD}.${PANELBIN}.bed.nochr
    fi
fi

# check panel bin
if [[ $PANELBIN == 'WGS' ]]; then
    PANELBIN='WGS'
elif [[ ! -f $PANELBIN ]]; then
    echo "ERROR: panel bin bed file ($PANELBIN) not found!"
    usage
    exit 1
else
    LINE=$(head -n1 $PANELBIN)
    if [[ $LINE == 'chr'* ]]; then
        P_CHR='Y'
    else
        P_CHR='N'
    fi
    
    # check if chr is in reference AND panel
    if [[ $G_CHR != $P_CHR ]]; then
        echo "ERROR: the presences of 'chr' in reference and panel bed are inconsistent"
        usage
        exit 1
    fi
fi

# check SNP marker database
if [[ $SNPmarker == 'hg19' ]]; then
    SNPmarker="$script_path/../data/SNPmarker_hg19"
elif [[ $SNPmarker == 'hg38' ]]; then
    SNPmarker="$script_path/../data/SNPmarker_hg38"
fi

if [[ ! -f $SNPmarker/chr1 ]]; then
    echo "ERROR: SNPmarker/chr1 ($SNPmarker) not found!"
    usage
    exit 1
fi

########## check optional arguments ##########
# check output directory
if [[ $OUTPUT != '.' ]] && [[ ! -d $OUTPUT ]]; then
    echo "ERROR: output directory ($OUTPUT) not found!"
    usage
    exit 1
fi

# check genome build
if [[ $GBUILD != 'hg38' ]] && [[ $GBUILD != 'hg19' ]]; then
    echo "ERROR: Genome build ($GBUILD) not supported!"
    usage
    exit 1
fi

# default PoN
if [[ $PON == 'WES' ]] || [[ $PON == 'cfDNA' ]]; then    
    PON=$script_path/../data/${PON}_PoN.txt
fi

# check PoN
if [[ $PON != 'none' ]] && [[ ! -f $PON ]]; then
    echo "ERROR: Panel of Normals ($PON) not found!"
    usage
    exit 1
fi

# covWindow
if [[ $covWindow == 'auto' ]]; then
    if [[ $DTYPE == 'WGS' ]]; then
        covWindow=1000000
    else
        covWindow=2000000
    fi
fi

# check gene list
if [[ $GENES == 'cfDNA' ]]; then    
    GENES=$script_path/../data/geneList.bed
fi

if [[ $GENES != 'none' ]] && [[ ! -f $GENES ]]; then
    echo "ERROR: Gene list ($GENES) not found!"
    usage
    exit 1
fi

#################################################################################
# Analysis
#################################################################################

cd $OUTPUT

WORK_DIR=${SAMPLE}.${DTYPE}.SubChrom        
# rename existing folder or mkdir a new one
if [[ -d ${WORK_DIR}.unpaired ]]; then
    mv ${WORK_DIR}.unpaired ${WORK_DIR}
elif [[ $preFile != 'True' ]]; then
    rm -rf $WORK_DIR
else
    mkdir -p $WORK_DIR
fi

cd $WORK_DIR


if [[ $DATA == *'bam' ]]; then
    # bedGraph
    echo -e "******************** run bedGraph ********************"
    if [[ ! -f ${SAMPLE}.${DTYPE}.bedGraph ]] ; then
        cmd="$script_path/bedGraph.sh -s $SAMPLE -i $DATA -d $DTYPE -p $PANELBIN"
        echo -e "$cmd\n"
        $cmd
    else
        echo -e "Found bedGraph file ...\n"
    fi
    
    # gatk_HaplotypeCaller
    echo -e "******************** run gatkHC ********************"
    if [[ ! -f ${SAMPLE}.${DTYPE}.gatkHC.vcf.gz ]]; then
        echo -e "gatk --java-options '-Xmx10G' HaplotypeCaller -R $REF -I $DATA -O ${SAMPLE}.${DTYPE}.gatkHC.vcf.gz"
        echo ""
        gatk --java-options "-Xmx10G" HaplotypeCaller -R $REF -I $DATA -O ${SAMPLE}.${DTYPE}.gatkHC.vcf.gz
    else
        echo -e "Found gatkHC file ...\n"
    fi

    # rename chrs and clean ENCODE_blacklist
    echo -e "******************** run SNP cleaning ********************"
    if [[ ! -f ${SAMPLE}.${DTYPE}.snp.txt ]]; then
        cmd="$script_path/rename_chrs_clean_blacklist.sh $SAMPLE ${SAMPLE}.${DTYPE}.gatkHC.vcf.gz $GBUILD $DTYPE $script_path"
        echo -e "$cmd\n"
        $cmd
    else
        echo -e "Found cleaned variant file ...\n"
    fi

elif [[ $DATA == *'vcf'* ]]; then
    echo -e "******************** run SNP cleaning ********************"    
    if [[ ! -f ${SAMPLE}.${DTYPE}.snp.txt ]]; then
        cmd="$script_path/rename_chrs_clean_blacklist.sh $SAMPLE $DATA $GBUILD $DTYPE $script_path"
        echo -e "$cmd\n"
        $cmd
    else
        echo -e "Found cleaned variant file ...\n"
    fi

elif [[ $DATA == *'.out' ]] && [[ ! -f ${SAMPLE}.${DTYPE}.snp.txt ]]; then
    if [[ $SAMPLE == *'_G'* ]]; then
        NorT="normal"
    else
        NorT="tumor"
    fi
    cmd="$script_path/load_clean_high20.sh $SAMPLE $DATA $DTYPE $GBUILD $NorT $script_path"
    echo -e "******************** load & clean high_20 ********************\n$cmd"
    $cmd
fi

# SubChrom
cmd="python $script_path/SubChrom_segmentation.py -s $SAMPLE -d $DTYPE -g $GBUILD -n $PON -cs $COVseg -rs $ROHseg -md $SNPmarker -minTF $minTF -minSize $minSize -minBins $minBins -gender $GENDER -dd $dipDep -covWin $covWindow -plotTF $plotTF -genes $GENES"
echo -e "******************** run SubChrom ********************\n$cmd"
$cmd
