#!/bin/bash

#################################################################################
# Default setting
#################################################################################

SAMPLE_LIST='NA'
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

QUEUE='standard'
LOG='SubChrom.run.log'
TASK_S=1
TASK_E=999999

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
    echo "Usage: $0 -s sample.list.txt -d cfDNA -r hg38 -p cfDNA -md hg38"
    echo ""
    echo -e "   -h  --help                Show this help message and exit"
    echo ""
    echo -e "Required arguments:"
    echo -e "   -s  --sample_list        Sample list with names and data paths. /path/to/sample.list.txt"
    echo -e "                               Format: <folder/patient_name><tab><sample_name><tab><data_path>"
    echo -e "   -d  --data_type          Sequencing data type. Options: WGS, WES, cfDNA, panel, custom, etc."
    echo -e "   -r  --reference          Genome reference (same for bam files). /path/to/reference.fa"
    echo -e "                               Options: 'hg38' for the CAB hg38, 'hg19' for the CompBio hg19"
    echo -e "   -p  --panel_bin          Bed file of your panel bins. /path/to/panel_bin.bed"
    echo -e "                               Options: 'WES' or 'cfDNA' to use pre-computed files, 'WGS' to skip"
    echo -e "   -md --marker_dir         Directory of SNP marker database from SubChrom. /path/to/SNPmarker"
    echo -e "                               Options: 'hg19' or 'hg38' for /SubChrom/data/SNPmarker_*/"
    echo ""
    echo -e "Optional arguments:"
    echo -e "   -o  --output             Ouput directory. Default: ./ (current directory)"
    echo -e "   -g  --genome_build       Options: hg38 (default), hg19"
    echo -e "   -n  --normal             Panel of Normals, or bedGraph file of a normal sample. /path/to/PoN.txt"
    echo -e "                               Options: 'cfDNA' to use pre-computed PoN. Default: none"
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
    echo -e "   -q  --queue               St. Jude HPC queue (Default: rhel8_standard)"
    echo -e "   -l  --log                 Analysis log file name (Default: $LOG)"
    echo -e "   -S  --task_start_line     Task N start to run (Default $TASK_S)"
    echo -e "   -E  --task_end_line       Task N stop to run  (Default $TASK_E)"
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
        -s | --sample_list)
            SAMPLE_LIST=$(realpath $2)
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
        -q | --queue)
            QUEUE=$2
            shift
            ;;
        -l | --log)
            LOG=$2
            shift
            ;;
        -S | --task_start_line)
            TASK_S=$2
            shift
            ;;
        -E | --task_end_line)
            TASK_E=$2
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
if [[ $SAMPLE_LIST == 'NA' ]]; then
    echo "Need to specify: -s --sample_list"
    usage
    exit 1
elif [[ ! -f $SAMPLE_LIST ]]; then
    echo "ERROR: sample_list ($SAMPLE_LIST) not found!"
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

# check the sample list
COLUMN=$(head -n1 $SAMPLE_LIST | awk '{print NF}')
if [[ $COLUMN != 3 ]]; then
    echo "ERROR: sample_list ($SAMPLE_LIST) need 3 columns!"
    usage
    exit 1
fi

ROW1=$(wc -l < $SAMPLE_LIST)
ROW2=$(sort -u $SAMPLE_LIST | wc -l)
if [[ $ROW1 != $ROW2 ]]; then
    echo "ERROR: Duplicate rows in sample_list ($SAMPLE_LIST)!"
    usage
    exit 1
fi

ROW1=$(cut -f1-2 $SAMPLE_LIST | wc -l)
ROW2=$(cut -f1-2 $SAMPLE_LIST | sort -u | wc -l)
if [[ $ROW1 != $ROW2 ]]; then
    echo "ERROR: Duplicate folder & sample names in sample_list ($SAMPLE_LIST)!"
    usage
    exit 1
fi

# default reference
if [[ $REF == 'hg19' ]]; then
    REF="/research/rgs01/applications/hpcf/authorized_apps/cab/Automation/REF/Homo_sapiens/NCBI/GRCh37-lite/bwa-index/GRCh37-lite_wochr/GRCh37-lite.fa"
elif [[ $REF == 'hg19MP' ]]; then
    REF="/home/slei/Group_dir/MolPath/hg19.fa"
elif [[ $REF == 'hg38' ]]; then
    REF="/research/rgs01/applications/hpcf/authorized_apps/cab/Automation/REF/Homo_sapiens/NCBI/GRCh38_no_alt/bwa-index/0.7.17-r1188/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
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
        PANELBIN=$script_path/../data/${PANELBIN}.${GBUILD}.chr.bed
    elif [[ $G_CHR == 'N' ]]; then
        PANELBIN=$script_path/../data/${PANELBIN}.${GBUILD}.nochr.bed
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
# loading modules
#################################################################################

module purge
source /research_jude/rgs01_jude/groups/cab/projects/automapper/common/slei/python_env/bin/activate
module load gatk/4.1.8.0
module load bcftools/1.14
module load deeptools/3.5.0
module load bedtools/2.30.0

#################################################################################
# Records
#################################################################################

echo ""
echo -e "[VERSION] SubChrom $VERSION"
echo -e "[DATE]    `date`"
echo -e "[SAMPLE]  $SAMPLE_LIST"
echo -e "[DTYPE]   $DTYPE"
echo -e "[REF]     $REF"
echo -e "[PANEL]   $PANELBIN"
echo -e "[SNP]     $SNPmarker"
echo ""
echo -e "[OUTPUT]  $OUTPUT"
echo -e "[GBUILD]  $GBUILD"
echo -e "[PON]     $PON"
echo -e "[COVseg]  $COVseg"
echo -e "[ROHseg]  $ROHseg"
echo -e "[minTF]   $minTF"
echo -e "[minSize] $minSize"
echo -e "[minBins] $minBins"
echo -e "[GNEDER]  $GENDER"
echo -e "[dipDep]  $dipDep"
echo -e "[covWin]  $covWindow"
echo -e "[plotTF]  $plotTF"
echo -e "[preFile] $preFile"
echo -e "[GENES]   $GENES"
echo ""
echo -e "[QUEUE]   $QUEUE"
echo -e "[LOGfile] $LOG"
echo -e "[TaskS]   $TASK_S"
echo -e "[TaskE]   $TASK_E"
echo -e "[G_CHR]   $G_CHR"
echo -e "[P_CHR]   $P_CHR"

echo -e "[Version] $VERSION" >>$LOG
echo -e "[DATE]    `date`" >> $LOG
echo -e "[CMD]     $0 -s $SAMPLE_LIST -d $DTYPE -r $REF -p $PANELBIN -md $SNPmarker -o $OUTPUT -g $GBUILD -n $PON -cs $COVseg -rs $ROHseg -mf $minTF -ms $minSize -mb $minBins -sg $GENDER -dd $dipDep -cw $covWindow -pf $plotTF -if $preFile -gl $GENES -q $QUEUE -l $LOG -S $TASK_S -E $TASK_E\n" >> $LOG

#################################################################################
# Analysis
#################################################################################

cd $OUTPUT

line_num=0
while IFS=$'\t' read FOLDER SAMPLE DATA; do
    
    line_num=$(( $line_num + 1 ))
    
    if [[ $line_num -ge $TASK_S ]] && [[ $line_num -le $TASK_E ]]; then
        echo ""
        echo -e "TASK#:   $line_num"
        echo -e "SAMPLE:  $SAMPLE"
        echo -e "DATA:    $DATA"

        mkdir -p $FOLDER
        cd $FOLDER
        
        WORK_DIR=${SAMPLE}.${DTYPE}.SubChrom
        # rename existing folder or mkdir a new one
        if [[ -d ${WORK_DIR}.unpaired ]]; then
            mv ${WORK_DIR}.unpaired ${WORK_DIR}
        else
            mkdir -p $WORK_DIR
        fi

        ################ CNV calling ################
        cmd="$script_path/../scripts/SubChrom.sh -s $SAMPLE -i $DATA -d $DTYPE -r $REF -p $PANELBIN -md $SNPmarker -g $GBUILD -n $PON -cs $COVseg -rs $ROHseg -mf $minTF -ms $minSize -mb $minBins -sg $GENDER -dd $dipDep -cw $covWindow -pf $plotTF -if $preFile -gl $GENES"
        echo -e "******************** run SubChrom ********************\n$cmd"
        bsub -P SubChrom -J ${SAMPLE}.SC -q $QUEUE -M 10000 -n 2 -R "span[hosts=1]" -eo ${WORK_DIR}/SC.err -oo ${WORK_DIR}/SC.out $cmd

        cd $OUTPUT
    fi
done < $SAMPLE_LIST
