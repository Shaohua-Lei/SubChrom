# Batch submission
To submit SubChrom jobs to HPC for a batch of samples, a sample list with names and data paths is required as input. Invoking the `-h/--help` flag will show options and formats.
```
Usage: ./StJude/SubChromSJ.sh -s sample.list.txt -d cfDNA -r hg38 -p cfDNA -n PoN.txt

   -h  --help                Show this help message and exit

Required arguments:
   -s  --sample_list        Sample list with names and data paths. /path/to/sample.txt
                               Format: <folder/patient_name><tab><sample_name><tab><data_path>
   -d  --data_type          Sequencing data type. Options: WGS, WES, panel, custom, etc.
   -r  --reference          Genome reference (same for bam files). /path/to/reference.fa
                               Options: 'hg38' for the CAB hg38, 'hg19' for the CompBio hg19
   -p  --panel_bin          Bed file of your panel bins. /path/to/panel_bin.bed
                               Options: 'WES' or 'cfDNA' to use pre-computed files, 'WGS' to skip
   -md --marker_dir         Directory of SNP marker database from SubChrom. /path/to/SNPmarker
                               Options: 'hg19' or 'hg38' for /SubChrom/data/SNPmarker_*/

Optional arguments:
   -o  --output             Ouput directory. Default: current directory
   -g  --genome_build       Options: hg38 (default), hg19
   -n  --normal             Panel of Normals, or bedGraph file of a normal sample. /path/to/PoN.txt
                               Options: 'cfDNA' to use pre-computed PoN. Default: none
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

   -q  --queue               St. Jude HPC queue (Default: standard)
   -l  --log                 Analysis log file name (Default: SubChrom.run.log)
   -S  --task_start_line     Task N start to run (Default 1)
   -E  --task_end_line       Task N stop to run  (Default 999999)
```

# Variant Calling
The most time comsuming step of SubChrom is the germline variant calling (gatk HaplotypeCaller). The default setting above takes 24-48h for WGS and 3-6h for WES/custom panel, and it varies by sequencing depth and panel size. To expedite this process, the `VariantCalling.sh` script splits the genome, submits variant calling jobs for individual chromomose, and then merges them into one file for SubChrom. This stategy reduces the time of variant calling by ~90%, and the bedGraph file of your samples will be generated with separate jobs as well.

  - Run `VariantCalling.sh` with your sample list. 
  - Wait and complete those jobs. 
  - Run `SubChromSJ.sh` with the same sample list, parameters, and output directory (if not current).

```
Usage: ./StJude/VariantCalling.sh -s sample.list.txt -d cfDNA -r hg38 -p cfDNA

	-h --help
	-s --sample_list         [Required] Sample list with sample names and data paths.
	                            Formats: <Folder_name><tab><Sample_name><tab><Sample_path>
	-d --data_type           [Required] Options: WGS, WES, cfDNA, panel, etc
	-r --reference           [Required] /path/to/reference.fa
	                            Provide genome reference used for the original bam file
	                            Options: 'hg38' to use the CAB version, 'hg19' to use the CompBio version
	-p --panel_bin           [Required] /path/to/PanelBin.bed
	                            Options: 'WES' or 'cfDNA' to use pre-computed files, 'WGS' to skip

	-g --genome_build        Options: hg38 (default), hg19
	-q --queue               St. Jude HPC queue (default: standard)
	-l --log                 Analysis log file name (default: SubChrom.run.log)
	-S --task_start_line     Task N start to run (default 1)
	-E --task_end_line       Task N stop to run  (default 999999)
```
