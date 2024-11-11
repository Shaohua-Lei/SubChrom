# *SubChrom*
A tool for dectecting ***Sub***clonal ***Chrom***osomal aberrations and estimating tumor fraction in cell-free DNA from next generation sequencing data.

## Introduction
SubChrom performs segmenation on coverage and variant allele frequency (VAF) for the detection of copy number variations (CNV), such as copy gain, copy loss, and copy neutral loss of heterozygosity (cnLOH). SubChrom is optimized to work on different types of data, including whole genome sequencing (WGS >15X), whole exome sequencing (WES), and especially custom panel sequencing. The estimate of tumor fraction is computed from the coverage and VAF changes of CNV events.
![SubChrom_fig](https://github.com/MullighanLab/SubChrom/assets/146115901/332c18c9-19b9-4297-be65-7a62af6356b5)



## Installation
1. Checkout the latest release of SubChrom from GitHub
    ```
    git clone git@github.com:MullighanLab/SubChrom.git
    ```
2. Download the SubChrom [database](https://zenodo.org/records/10155688) to the `/path/to/SubChrom/data` directory and unzip it.
3. Install dependencies
   - [gatk (>=4.1.8.0)](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)
   - [bcftools (>=1.14)](https://samtools.github.io/bcftools/bcftools.html)
   - [deeptools (>=3.5.0)](https://deeptools.readthedocs.io/en/develop/content/installation.html)
   - [bedtools (>=2.30.0)](https://bedtools.readthedocs.io/en/latest/content/installation.html)
   - Python3 with packages: [pybedtools](https://daler.github.io/pybedtools/main.html), pandas, numpy, matplotlib, scipy, copy, argparse, os, shutil

## Panel Bins
For custom panel sequencing data, the first step is to generate a set of target bins of your panel. The default size of target bins is 500bp, and it can be tuned using the option `-s/--size`. A smaller bin size will yield more target bins, presumably resulting in higher resolution for coverage segmentation with noisier bins.
```
python /path/to/SubChrom/scripts/PanelBin.py -i input.bed -o output.bed -g hg38 -s 500
```
For WES data, pre-computed files of exon bins are available in the database directory as `/data/*WES.bed*` if downloaded above.

For WGS data, panel bins are not required.

## Running SubChrom
Simply use the `SubChrom.sh` script provided in the `SubChrom/scripts/` directory. Here is a short example of how to launch the script from the command line:
```
/path/to/SubChrom/scripts/SubChrom.sh -s filename -i sample.bam -d cfDNA -r hg38.fa -p panel_bin.bed -md hg38
```

Invoking the `-h/--help` flag will show the list of options.
```
Usage: ./SubChrom.sh -s filename -i sample.bam -d cfDNA -r hg38.fa -p panel_bin.bed -md hg38

   -h  --help                Show this help message and exit

Required arguments:
   -s  --sample             Unique sample name
   -i  --input              /path/to/sample.data
                               Format: BAM from WES and panel sequencing, BAM/vcf/high20 from WGS
   -d  --data_type          Sequencing data type. Options: WGS, WES, panel, custom, etc.
   -r  --reference          Genome reference (same for bam files). /path/to/reference.fa
   -p  --panel_bin          Bed file of your panel bins. /path/to/panel_bin.bed
                               Options: 'WES' to use files from /SubChrom/data/, 'WGS' to skip
   -md --marker_dir         Directory of SNP marker database from SubChrom. /path/to/SNPmarker
                               Options: 'hg19' or 'hg38' for /SubChrom/data/SNPmarker_*/

Optional arguments:
   -o  --output             Ouput directory. Default: ./ (current directory)
   -g  --genome_build       Options: hg38 (default), hg19
   -n  --normal             Panel of Normals, or bedGraph file of a normal sample. /path/to/PoN.txt
                               Default: none
   -cs --coverage_seg       Perform coverage segmentation. Options: True, False (default)
                               If --normal is none above for WGS/panel data, this is False by default
   -mf --minTF              Minimal tumor fraction to report a CNV event. Default: 0.1. Minimum: 0.01
   -ms --minSize            Minimal size (bp) to report a CNV event. Default: 1000000. Minimum: 10000
   -mb --minBins            Minimal number of bins to report a CNV event. Default: 100. Minimum: 10
   -sg --sample_gender      Options: Male/M, Female/F, auto. Default: auto for automatic detection
   -dd --diploid_depth      How to compute the diploid depth.
                               Options: auto, chr1...chr22, chrX, a specific value such as 500
                               Default: auto for automatic optimization
   -cw --covWinSize         Coverage window size (bp) for visualization. Default: 2000000. Minimum: 500000
   -pf --plotTF             Plot tumor fraction or not. Options: True (default), False
   -if --intermediate_file  Use intermediate files from the previous run, such as vcf file and segements
                               Options: True (default), False (remove previous files and generate new ones)
   -gl --gene_list          Gene list of interest for visualization. /path/to/geneList.bed
                               Default: /SubChrom/data/geneList.bed
```

For local users at St. Jude, SubChrom [wrapper](./StJude/README.md) scripts provided in the `StJude` directory are available for batch job submissions on the HPC cluster.

## Normal or Panel of Normals
Normal or panel of normals (PoN) is not needed for WGS, while it is not required but highly recommended for WES and custom panel sequencing. To enable coverage segmentation and improve estimation of tumor fraction, normal/PoN is used to normalize coverage data from WES or custom panel.

To use a normal for your target/tumor sample, the coverage bedGraph file of your normal sample can be added by option `-n/--normal` when running SubChrom above. The bedGraph file of any sample can be found after running SubChrom. Alternatively, the `bedGraph.sh` script can generate one (see usage by `-h/--help`). Please note that your target and normal samples do not have to be an exact pair, they can come from different individuals. In addition, your normal can be a tumor sample, if it has no copy loss/gain events.
```
/path/to/SubChrom/scripts/bedGraph.sh -s sample_ID -i sample.bam -d WES -p panel_bin.bed -o .
```

To create a PoN file, a bedGraph file for each sample in your PoN should be generated as above. Then, use the `PanelOfNormals.py` script to take your list of bedGraph files (see format in usage by `-h/--help`).
```
python /path/to/SubChrom/scripts/PanelOfNormals.py -s samples.txt -g hg38 -o PoN.txt
```

## Output
### 1. Fields of *SubChrom.txt*
| Field    | Description     |
|----------|-------------|
| CHR   | Chromosome   |
| START   | Start position   |
| END   | End position   |
| SIZE   | Event size   |
| covMean   | Mean coverage   |
| covLog2Ratio   | Log2 ratio of coverage   |
| EVENT   | CNV event   |
| Allelic_Imbalance   | Presence of allelic imbalance   |
| chrA_freq   | Event frequency of chrA's CNV   |
| chrA_CNV   | Integer value of chrA's CNV   |
| chrA_TF   | Estimated tumor fraction based on chrA's CNV   |
| chrB_freq   | Event frequency of chrB's CNV   |
| chrB_CNV   | Integer value of chrB's CNV   |
| chrB_TF   | Estimated tumor fraction based on chrB's CNV   |

Note: To differentiate the potential CNV event(s) of each chromosome from the diploid pair, chrA & chrB are annotated separately. For example, the CNV events of chrA/chrB in a cnLOH event are -1/+1, and a two-copy gain event can be +1/+1 or NA/+2.

### 2. Visualization in *CNV.png*
![SJBALL032020_F1 cfDNA CNV](https://github.com/Shaohua-Lei/SubChrom/assets/146115901/bfc1e4f5-3b4d-4aa4-9f0b-1b7ee02da066)
  - First panle: Mean coverage with a window size of 2M bp. Dotted line indicates the sequencing depth of diploid, which can be tuned using the option `-dd/--diploid_depth`.
  - Second panel: VAF distribution across the genome.
  - Third panel: Copy ratio of SNP markers (WGS) or panel bins (WES/custom panel) in grey, and CNV events in different colors.
  - Fourth panel: Estimated tumor fraction based on CNV events.

### 3. Visualization in *focal.png*
An enlarged view of the second & third panels of *CNV.png* to visualize focal events, and genes of interest can be annotated in this figure using the option `-gl/--gene_list`.

## Citation
Lei, S., Jia, S., Takalkar, S. et al. Genomic profiling of circulating tumor DNA for childhood cancers. Leukemia (2024). https://doi.org/10.1038/s41375-024-02461-x

## Contacts
For questions and feedback, please go to the *Discussions* or *Issues* section, or email <shaohua.lei@stjude.org>

## Acknowledgements
This tool was developed in collaboration with:
 - Center of Excellence for Leukemia Studies, St. Jude Children's Research Hospital
 - Clinical Biomarkers Laboratory, St. Jude Children's Research Hospital
 - Center for Applied Bioinformatics, St. Jude Children's Research Hospital

## License
SubChrom Copyright (C) 2023 St. Jude Children's Research Hospital

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/. For additional questions about the License, please contact the St. Jude Office of Technology Licensing at scott.elmer@stjude.org.
