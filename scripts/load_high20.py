#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys

# load and pre-process high20 file
# sample is normal or tumor
def load_high20(path, sample='tumor'):
    high20 = pd.read_csv(path,sep = '\t',
                           usecols = ["Chr", 'Pos', "Type",'Chr_Allele', 'Alternative_Allele', 
                                      'reference_normal_count', 'reference_tumor_count',
                                      'alternative_normal_count', 'alternative_tumor_count'])
    high20.columns = ["CHR", 'POS', "Type",'REF', 'ALT', 'REF_G', 'REF_D', 'ALT_G', 'ALT_D']
    df = high20[high20.Type == 'SNP'].copy()
    # remove 'chr' from chromosome names
    if df.iloc[0,0][:3] == 'chr': 
        df.CHR = pd.DataFrame(df.CHR.str[3:])

    df['POS_0'] = df['POS'] - 1
    if sample=='normal':
        df['COV'] = df.REF_G + df.ALT_G
        df = df[df.COV >= 8]
        df = df[["CHR", 'POS_0', 'POS', 'REF', 'ALT', "COV", 'REF_G', 'ALT_G']]
    elif sample =='tumor':
        df['COV'] = df.REF_D + df.ALT_D
        df = df[df.COV >= 8]
        df = df[["CHR", 'POS_0', 'POS', 'REF', 'ALT', "COV", 'REF_D', 'ALT_D']]

    return df

if __name__ == "__main__":
    SAMPLE = sys.argv[1]
    DATA  = sys.argv[2]
    DTYPE = sys.argv[3]
    NorT  = sys.argv[4]
    
    df = load_high20(DATA, NorT)
    df.to_csv(SAMPLE + '.' + DTYPE + '.snpRAW.txt', sep='\t', index=False, header=False)
