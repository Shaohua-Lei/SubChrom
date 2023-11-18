#!/bin/bash

#################################################################################
# Default setting
#################################################################################

SAMPLE_LIST='NA'
DTYPE='NA'
REF='NA'
PANELBIN='NA'

GBUILD='hg38'
QUEUE='standard'
LOG='SubChrom.run.log'
TASK_S=1
TASK_E=999999

VERSION='0.1.0 (September 29, 2023)'

#################################################################################
# Display usage
#################################################################################

usage()
{
    echo ""
    echo "Program: SubChrom variant calling"
    echo -e "Version: $VERSION"
    echo ""
    echo "Usage: $0 -s sample.list.txt -d cfDNA -r hg38 -p cfDNA"
    echo ""
    echo -e "\t-h --help"
    echo -e "\t-s --sample_list         [Required] Sample list with sample names and data paths."
    echo -e "\t                            Formats: <Folder_name><tab><Sample_name><tab><Sample_path>"
    echo -e "\t-d --data_type           [Required] Options: WGS, WES, cfDNA, panel, etc"
    echo -e "\t-r --reference           [Required] /path/to/reference.fa"
    echo -e "\t                            Provide genome reference used for the original bam file"
    echo -e "\t                            Options: 'hg38' to use the CAB version, 'hg19' to use the CompBio version"
    echo -e "\t-p --panel_bin           [Required] /path/to/PanelBin.bed"
    echo -e "\t                            Options: 'WES' or 'cfDNA' to use pre-computed files, 'WGS' to skip"
    echo ""
    echo -e "\t-g --genome_build        Options: hg38 (default), hg19"
    echo -e "\t-q --queue               St. Jude HPC queue (default: standard)"
    echo -e "\t-l --log                 Analysis log file name (default: $LOG)"
    echo -e "\t-S --task_start_line     Task N start to run (default $TASK_S)"
    echo -e "\t-E --task_end_line       Task N stop to run  (default $TASK_E)"
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
            SAMPLE_LIST=$2
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
        -g | --genome_build)
            GBUILD=$2
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

# code & work directory
script_path="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
main_dir=$(pwd)

# check required inputs
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
elif [[ $GBUILD != 'hg38' ]] && [[ $GBUILD != 'hg19' ]]; then
    echo "ERROR: Genome build ($GBUILD) not supported!"
    usage
    exit 1
fi

# default reference
if [[ $REF == 'hg19' ]]; then
    REF="/research/rgs01/applications/hpcf/authorized_apps/cab/Automation/REF/Homo_sapiens/NCBI/GRCh37-lite/bwa-index/GRCh37-lite_wochr/GRCh37-lite.fa"
elif [[ $REF == 'hg19MP' ]]; then
    REF="/research/rgs01/applications/hpcf/authorized_apps/cab/Automation/REF/Homo_sapiens/NCBI/GRCh37-lite/bwa-index/0.7.17-r1188/GRCh37-lite_wchr.fa"
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

#################################################################################
# loading modules
#################################################################################

module load gatk/4.1.8.0
module load deeptools/3.5.0 # bw and bG files
module load bcftools/1.14   # concat vcf files

#################################################################################
# Records
#################################################################################

echo -e "[VERSION] $VERSION"
echo -e "[DATE]    `date`"
echo -e "[SAMPLE]  $SAMPLE_LIST"
echo -e "[DTYPE]   $DTYPE"
echo -e "[REF]     $REF"
echo -e "[PANEL]   $PANELBIN"
echo -e "[GBUILD]  $GBUILD"
echo -e "[LOGfile] $LOG"
echo -e "[TaskS]   $TASK_S"
echo -e "[TaskE]   $TASK_E"
echo -e "[G_CHR]   $G_CHR"
echo -e "[P_CHR]   $P_CHR"

echo -e "[Version] $VERSION" >>$LOG
echo -e "[DATE]    `date`" >> $LOG
echo -e "[CMD]     $0 -s $SAMPLE_LIST -d $DTYPE -r $REF -p $PANELBIN -g $GBUILD -q $QUEUE -l $LOG -taskS $TASK_S -taskE $TASK_E\n" >> $LOG

#################################################################################
# Analysis
#################################################################################

# number of columns in the sample list
COLUMN=$(head -n1 $SAMPLE_LIST | awk '{print NF}')

line_num=0
while IFS=$'\t' read FOLDER SAMPLE DATA DIP_CHR; do      
    # 3 or 4 columns in the sample list
    if [[ $COLUMN == 3 ]]; then
        DIP_CHR='auto'
    fi

    line_num=$(( $line_num + 1 ))

    if [[ $line_num -ge $TASK_S ]] && [[ $line_num -le $TASK_E ]]; then
        echo ""
        echo -e "TASK:    $line_num"
        echo -e "SAMPLE:  $SAMPLE"
        echo -e "DATA:    $DATA"

        if [[ ! -f $DATA ]] ; then
            echo "Input file $DATA [BAM] not found!"
            exit 1
        fi

        mkdir -p $FOLDER
        cd $FOLDER

        WORK_DIR=${SAMPLE}.${DTYPE}.SubChrom
        
        # rename existing folder or mkdir a new one
        if [[ -d ${WORK_DIR}.unpaired ]]; then
            mv ${WORK_DIR}.unpaired ${WORK_DIR}
        else
            mkdir -p $WORK_DIR
        fi
        
        cd $WORK_DIR

        ################ BAM to bigwig and bedGraph ################
        if [[ $DATA == *'bam' ]] && [[ ! -f ${SAMPLE}.${DTYPE}.bedGraph ]]; then
            cmd="$script_path/../scripts/bedGraph.sh -s $SAMPLE -i $DATA -d $DTYPE -p $PANELBIN"
            echo -e "******************** run bedGraph ********************\n$cmd"
            bsub -P SubChrom -J ${SAMPLE}.bG -q standard -M 3000 -eo bG.err -oo bG.out $cmd
        fi

        ################ Variant calling ################
        if [[ $DATA == *'bam' ]] && [[ ! -f ${SAMPLE}.${DTYPE}.gatkHC.vcf.gz ]]; then
            cmd="$script_path/gatkHC.sh $SAMPLE $DATA $DTYPE $REF $G_CHR"
            echo -e "******************** run gatkHC ********************\n$cmd"
            bsub -P SubChrom -J ${SAMPLE}.HC -q $QUEUE -M 4000 -n 2 -R "span[hosts=1]" -eo HC.err -oo HC.out $cmd
        fi
        cd $main_dir
    fi
done < $SAMPLE_LIST
