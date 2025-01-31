#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import sys
import os

def bedGraph(path):
    df = pd.read_csv(path, sep='\t')
    df.columns = ["CHR", "POS", "END", "COV"]
    # zero low coverage regions
    lowCOV = df[df.COV <= 3].index
    df.loc[lowCOV, "COV"] = 0
    return df

def computeCOV(df1):
    df = df1[(df1.CHR != 'chrX') & (df1.CHR != 'chrY')].copy()
    covCount = len(df)
    data = df.COV.values
    data = np.sort(data)
    quantile = int(covCount/4)
    mm50 = (data[quantile:-quantile]).mean()
    return round(mm50, 5)

def normalization(sample, gender, path, build):
    df = bedGraph(path)
    mm50 = computeCOV(df)
    
    df[sample] = df.COV/mm50
    if gender in ['Male', 'M']:
        if build == 'hg38':
            chrX = df[(df.CHR == 'chrX') & (df.POS > 2781479) & (df.POS < 155701383)].index
        elif build == 'hg19':
            chrX = df[(df.CHR == 'chrX') & (df.POS > 2699520) & (df.POS < 154931044)].index
        else:
            print("Genome build", build, 'not supported !!!')
            return
        df.loc[chrX, sample] = df.loc[chrX, 'COV']*2 /mm50
        chrY = df[(df.CHR == 'chrY')].index
        df.loc[chrY, sample] = df.loc[chrY, 'COV']*2 /mm50
    else:
        chrY = df[(df.CHR == 'chrY')].index
        df.loc[chrY, sample] = np.nan
        
    return mm50, df

class NewlineHelpFormatter(argparse.RawTextHelpFormatter):
    def _split_lines(self, text, width):
        # Preserve newlines in the help message
        if text.startswith("R|"):
            return text[2:].splitlines()
        return argparse.RawTextHelpFormatter._split_lines(self, text, width)

def main():
    # Create an ArgumentParser with the custom formatter
    parser = argparse.ArgumentParser(
        description="Create a Panel of Normals for coverage normalization in SubChrom.",
        formatter_class=NewlineHelpFormatter
    )

    # Define arguments and options
    parser.add_argument("-s", "--samples", type=str, required=True, 
        help="""R|A list of samples with sample name, gender, path to bedGraph file. No header.
    Format: <Sample_name><tab><Gender><tab></path/to/bedGraph>
    Gender options: Male/M, Female/F"""
    )
    parser.add_argument("-g", "--genome", type=str, required=True, 
                        help="Genome build. Options: hg38, hg19")
    parser.add_argument("-o", "--output", type=str, required=True, 
                        help="Output path. /path/to/PoN.txt")

    # Parse command-line arguments
    args = parser.parse_args()
    dfSamples = pd.read_csv(args.samples, sep='\t', names=['Sample', 'Gender', 'bedGraph'])
    build = args.genome
    
    if build not in ['hg19', 'hg38']:
        print("Stopping the code execution, genome build not supported:", build)
        sys.exit()

    # start analysis
    depth = []
    df = pd.DataFrame()

    for i in dfSamples.index:
        sample, gender, path = dfSamples.loc[i, 'Sample'], dfSamples.loc[i, 'Gender'], dfSamples.loc[i, 'bedGraph']
        mm50, dfCOV = normalization(sample, gender, path, build)
        
        if gender not in ['Male', 'M', 'Female', 'F']:
            print("Stopping the code execution, gender not supported:", gender)
            sys.exit(0)
        elif not os.path.exists(path):
            print("Stopping the code execution, bedGraph not found:", path)
            sys.exit(0)
        
        # keep bin coordinates
        if i == 0:
            dfCOV0 = dfCOV[["CHR", "POS", "END"]].copy()

        # exclude samples of low depth
        if mm50 > 30: 
            print('Computing:', sample, gender, '| Sequencing depth:', mm50)   
            depth.append(mm50)
            df = pd.concat([df, dfCOV[[sample]]], axis=1)
     
    df['MEDIAN'] = round(df.median(axis=1), 5)   
    df['MEAN'] = df.mean(axis=1)
    df['SEM'] = df.std(axis=1)/df.MEAN

    # combine df
    df = pd.concat([dfCOV0, df], axis=1)
   
    # remove low coverage regions
    lowCOV = df[df.SEM > 3].index
    df.loc[lowCOV, "MEDIAN"] = np.nan

    # save index file
    df[["CHR", "POS", "END", 'MEDIAN']].to_csv(args.output,sep='\t', index=False, header=True)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nScript terminated by user.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
