#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import pybedtools
import sys

def main():
    parser = argparse.ArgumentParser(description="SubChrom panel binning")
    parser.add_argument("-i", "--input", dest='input', type=str, required=True,
                       help="/path/to/Panel.bed. BED format, no header")
    parser.add_argument("-o", "--output", dest='output', type=str, required=True,
                       help='/path/to/PanelBin.bed')
    parser.add_argument("-g", "--genome", dest='genome', type=str, default='hg38',
                       help="Genome build. Options: hg38 (default), hg19")
    parser.add_argument("-s", "--size", dest='binsize', type=int, default=500,
                       help='Size of bin (bp) for the input panel. Default: 500')
    
    args = parser.parse_args()
    Panel = args.input
    PanelBin  = args.output
    binSize = args.binsize
    build = args.genome
    
    print('***** SubChrom panel bining *****')
    print('Input file:', Panel)
    print('Bin size (bp):', binSize)
    print('Genome build:', build)

    df = pd.read_csv(Panel, sep='\t', header=None, dtype={0: str})

    # with or without CHR
    if df.iloc[3,0][:3] == 'chr':
        ENCODE = '../data/' + build + '-blacklist.v2.bed.gz'
        print("'chr' detected in the input file")
    else:
        ENCODE = '../data/' + build + '-blacklist.v2.bed.nochr'
        print("'chr' not detected in the input file")
    print('Blacklist:', ENCODE)

    ################ subtract ENCODE blacklist ################
    PANEL = pybedtools.BedTool(Panel)
    BLACK = pybedtools.BedTool(ENCODE)
    clean = PANEL.subtract(BLACK, A=False)
    df = clean.to_dataframe()
    # remove extra columns
    df = df[['chrom', 'start', 'end']]
    df.columns = ['CHR','START','END']

    ################ bining ################
    dfBED = pd.DataFrame()
    for i in df.index:
        CHR, START, END = df.loc[i, 'CHR'], df.loc[i, 'START'], df.loc[i, 'END']
        # small region, append the row
        if END - START <= binSize*1.2:
            dfROW = df.loc[[i]].copy()     
        else:
            bins = np.arange(START, END+1, step=binSize)
            # the last bin
            if END - bins[-1] >= binSize*0.2:
                bins = np.append(bins, END+1)
            else:
                bins = np.append(bins[:-1], END+1)
            # convert into df
            dfROW = pd.DataFrame([bins[:-1], bins[1:]-1], index=['START','END']).T
            dfROW['CHR'] = CHR
            dfROW = dfROW[['CHR','START','END']]  
        dfBED = pd.concat([dfBED, dfROW], axis=0)

    # save
    dfBED.to_csv(PanelBin, index=False, header=False, sep='\t')
    
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nScript terminated by user.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
