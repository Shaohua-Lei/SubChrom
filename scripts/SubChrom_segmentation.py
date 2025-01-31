#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.text import Text
import matplotlib.patches as patches
from matplotlib.ticker import FuncFormatter
import scipy.stats as sts
import scipy.optimize as opt
import copy
import argparse
import os
import shutil
import pybedtools
import sys

################################################################
# load variant file and keep good SNP markers
################################################################
def genomeDatabase(build='hg38'):
    # build: only hg19 or hg38
    # chromosome length of hg19, hg38, and color
    genomeSize={
        '1':  [249250621, 248956422, '#3c5fb0'],
        '2':  [243199373, 242193529, '#ff823e'],
        '3':  [198022430, 198295559, '#368019'],
        '4':  [191154276, 190214555, '#cc4f43'],
        '5':  [180915260, 181538259, '#874aa8'],
        '6':  [171115067, 170805979, '#70521b'],
        '7':  [159138663, 159345973, '#cc5bc8'],
        '8':  [146364022, 145138636, '#757575'],
        '9':  [141213431, 138394717, '#b0b558'],
        '10': [135534747, 133797422, '#33b8b4'],
        '11': [135006516, 135086622, '#3c5fb0'],
        '12': [133851895, 133275309, '#ff823e'],
        '13': [115169878, 114364328, '#368019'],
        '14': [107349540, 107043718, '#cc4f43'],
        '15': [102531392, 101991189, '#874aa8'],
        '16': [90354753,  90338345,  '#70521b'],
        '17': [81195210,  83257441,  '#cc5bc8'],
        '18': [78077248,  80373285,  '#757575'],
        '19': [59128983,  58617616,  '#b0b558'],
        '20': [63025520,  64444167,  '#33b8b4'],
        '21': [48129895,  46709983,  '#3c5fb0'],
        '22': [51304566,  50818468,  '#ff823e'],
        'X':  [155270560, 156040895, '#368019'],
        'Y':  [59373566,  57227415,  '#cc4f43']}
    # convert into a dataframe
    genoInfo = pd.DataFrame(genomeSize, index=['hg19', 'hg38', 'Color']).T
    # compute accumulative length
    accumLen = 0
    for CHR in genoInfo.index:
        genoInfo.loc[CHR, 'prevLen'] = accumLen
        # 6e6 to separate two chr, 6e6 may be the best (not too small not too large)
        accumLen = accumLen + genoInfo.loc[CHR, build] + 6e6
    # clean up and standardize df
    if build == 'hg19':
        genoInfo.drop('hg38', axis=1, inplace=True)
        genoInfo.rename({'hg19':'chrSize'}, axis=1, inplace=True)
    elif build == 'hg38':
        genoInfo.drop('hg19', axis=1, inplace=True)
        genoInfo.rename({'hg38':'chrSize'}, axis=1, inplace=True)
    # plot positions of chr labels on x-axis
    genoInfo['plotPos'] = genoInfo['prevLen'] + genoInfo.chrSize/2
    return genoInfo

def load_variants(path):
    print("********** Loading variants raw data ... **********")
    df = pd.read_csv(path, sep='\t', names=['CHR','POS','REF','ALT','COV','COV_REF','COV_ALT'],
                     dtype={'CHR':str,'POS':int,'REF':str,'ALT':str,'COV':int,'COV_ALT':int})
    df = df[(df.COV>=8) & (df.COV_ALT>0)]
    # remove 'chr' from chromosome names
    if df.iloc[0,0][:3] == 'chr': 
        df.CHR = pd.DataFrame(df.CHR.str[3:])
    # variant full name
    df['variant'] = df.CHR + '.' + df.POS.astype(str) + '.' + df.REF + '.' + df.ALT
    return df

def filter_variants(variants, build, SNPmarker):
    print("********** Filtering variants ... **********")
    dictSNVall = {}
    genoInfo = genomeDatabase(build)
    for CHR in genoInfo.index:
        #print("Preparing chr" + CHR + " ...")
        marker = pd.read_csv(SNPmarker + '/chr' + CHR, names=['variant'])
        dictSNVall[CHR] = pd.merge(variants[variants.CHR == CHR], marker)
    return dictSNVall

def SNPmarkers(path, build, workDir, SNPmarker):
    SNP_markers_chr = workDir + '/SNP_markers_chr'
    # load complete markers
    if os.path.exists(SNP_markers_chr + '/complete.txt'):
        print("********** Loading markers ... **********")
        dictSNVall = {}
        genoInfo = genomeDatabase(build)
        for CHR in genoInfo.index:
            dictSNVall[CHR] = pd.read_csv(SNP_markers_chr + '/chr' + CHR, sep='\t', dtype={'CHR':str})
        return dictSNVall
    
    # remove incomplete markers
    if os.path.exists(SNP_markers_chr):
        shutil.rmtree(SNP_markers_chr)
    
    # generate markers
    os.mkdir(SNP_markers_chr)
    variants = load_variants(path)
    dictSNVall = filter_variants(variants, build, SNPmarker)
    
    # save markers
    print("********** Saving markers ... **********")
    for CHR in dictSNVall:
        dictSNVall[CHR].to_csv(SNP_markers_chr + '/chr' + CHR, sep='\t', index=False)
    with open(SNP_markers_chr + '/complete.txt', 'w') as complete:
        complete.write("0\n")
    return dictSNVall
        
def filter_markers(dictSNVall, build, dataType, bedGraph):
    print("********** Filtering markers ... **********")
    # minCOV: minimal coverage to keep a marker.
    # minMAC: minimal minor allele count to keep a heterozygote.
    if dataType == 'WGS': # WGS variants has been filtered before
        dictRawDepth = []
        for CHR in range(1,23):
            dictRawDepth.append(dictSNVall[str(CHR)].COV.median())
        rawDepth = np.median(dictRawDepth)
        minCOV,minMAC = 8,3
    else:
        # compute rawDepth
        dfBG = pd.read_csv(bedGraph, sep='\t',header=0, names=["CHR", "POS", "END", "COV"], dtype={'CHR':str})
        if dfBG.iloc[0,0][:3] == 'chr': 
            dfBG.CHR = dfBG.CHR.str[3:]
        dfBG = dfBG[(dfBG.CHR != 'X') & (dfBG.CHR != 'Y')]
        rawDepth = dfBG.COV.median()
        minCOV = np.max([15, round(rawDepth/30)]) # >=15 for WES and cfDNA data
        minMAC = round(minCOV/3)
        
    if dataType == 'WGS' and rawDepth < 15:
        print("Low coverage data (<15x) failed SubChrom:", rawDepth)
        sys.exit(1)
            
    print('dataType:', dataType)
    print('rawDepth:', rawDepth)
    print('minCOV to keep a SNP marker:', minCOV)
    print('minMAC to keep a heterozygotes:', minMAC)
     
    # filter markers
    dictSNV = {}
    genoInfo = genomeDatabase(build)
    for CHR in dictSNVall:
        # filter COV
        df = dictSNVall[CHR].copy()
        df['VAF'] = np.nan # define VAF for empty df, such as chrY
        df = df[df.COV >= minCOV]
        # compute VAF
        lowVAF = df[df.COV_ALT < minMAC].index
        df.loc[lowVAF, 'VAF'] = 0
        highVAF = df[df.COV_REF < minMAC].index
        df.loc[highVAF, 'VAF'] = 1
        normVAF = df[(df.COV_ALT >= minMAC) & (df.COV_REF >= minMAC)].index
        df.loc[normVAF, 'VAF'] = round(df.loc[normVAF, 'COV_ALT'] / df.loc[normVAF, 'COV'], 5)
        # PlotPos
        df = pd.merge(df, genoInfo[['prevLen']], left_on='CHR', right_index=True, how='left')
        df['PlotPos'] = df.POS + df.prevLen
        df.reset_index(inplace=True)
        df.drop(['prevLen', 'index'],axis=1,inplace=True)
        dictSNV[CHR] = df
        
        # genome info
        # WGS: median of snpCount_1Mb = 1136, hetCount_1Mb = 683, hetRatio = 0.6
        if len(df)>0:
            chrSize = genoInfo.loc[CHR, 'chrSize']
            genoInfo.loc[CHR, 'snpCount'] = len(df)
            genoInfo.loc[CHR, 'snpCount_1Mb'] = round(len(df)*1e6/chrSize,3)
            genoInfo.loc[CHR, 'hetCount'] = ((df.VAF>0) & (df.VAF<1)).sum()
            genoInfo.loc[CHR, 'hetCount_1Mb'] = round(genoInfo.loc[CHR, 'hetCount']*1e6/chrSize, 3)
            vafMeanCHR = df[(df.VAF>0) & (df.VAF<1)].VAF.mean()
            # vafMeanCHR can be >>0.5, such as in Loss or cn_LOH events
            genoInfo.loc[CHR, 'vafMeanCHR'] = np.min([vafMeanCHR, 0.5]) 
            genoInfo.loc[CHR, 'hetRatio'] = genoInfo.loc[CHR, 'hetCount'] / genoInfo.loc[CHR, 'snpCount']
        # chrY can be empty
        else:
            genoInfo.loc[CHR, 'snpCount'] = 0
            genoInfo.loc[CHR, 'snpCount_1Mb'] = 0
            genoInfo.loc[CHR, 'hetCount'] = 0
            genoInfo.loc[CHR, 'hetCount_1Mb'] = 0
            genoInfo.loc[CHR, 'hetRatio'] = 0
    return dictSNV, genoInfo, rawDepth

def get_gender(GENDER, genoInfo, dataType):
    print("********** Estimating gender and chrY status ... **********")
    if GENDER in ['Male', 'M']:
        gender, chrY = 'M', 'Unknown'
    elif GENDER in ['Female', 'F']:
        gender, chrY = 'F','no_Y'
    
    # Estimates based on data
    elif dataType == 'WGS':
        if genoInfo.loc['Y','snpCount'] > 30:
            gender, chrY = 'M', 'Y'
        # most likely a male, but also possibly a female with X_loss
        elif genoInfo.loc['X','hetRatio'] < 0.2: 
            gender, chrY = 'M', 'Y_loss'
        else:
            gender, chrY = 'F','no_Y'
    else:
        if genoInfo.loc['Y','snpCount'] > 30:
            gender, chrY = 'M', 'Y'
        # most likely a male, but also possibly a female with X_loss
        elif genoInfo.loc['X','hetRatio'] < 0.2:
            gender, chrY = 'M','Unknown'
        else:
            gender, chrY = 'F','no_Y'
    print('Gender:', gender)
    print('chrY status:', chrY)
    return gender, chrY

################################################################
# compute coverage
################################################################
def load_coverage(dictSNV, bedGraph, genoInfo, dataType, PON):
    print("********** Loading coverage raw data ... **********")
    # no need to filter high coverage regions, which may be local amplications such as myc/mycn
    if dataType == 'WGS':
        print('WGS: dictCOV = dictSNV')
        return dictSNV
        
    if PON == 'none':
        print('Returning raw coverage ...')
        df = pd.read_csv(bedGraph, sep='\t',header=0, names=["CHR", "POS", "END", "COV"], dtype={'CHR':str})
        df = df[df.COV >= 8]
        if df.iloc[0,0][:3] == 'chr': 
            df.CHR = df.CHR.str[3:]
    else:
        print('Normalizing coverage ...')
        df = pd.read_csv(bedGraph, sep='\t',header=0, names=["CHR", "POS", "END", "COV_raw"], dtype={'CHR':str})
        if df.iloc[0,0][:3] == 'chr': 
            df.CHR = df.CHR.str[3:]
        
        dfIndex = pd.read_csv(PON, sep='\t',header=0, names=["CHR", "POS", "END", "MEDIAN"], dtype={'CHR':str})
        print('Median coverage of PON:', round(dfIndex.MEDIAN.median()))
        dfIndex.MEDIAN = dfIndex.MEDIAN/dfIndex.MEDIAN.median()
        # remove low coverage regions
        lowCOV = dfIndex[(dfIndex.MEDIAN < 0.001)].index
        dfIndex.loc[lowCOV, "MEDIAN"] = np.nan
        
        df = pd.concat([df, dfIndex[["MEDIAN"]]], axis=1)
        df['COV'] = df.COV_raw/df.MEDIAN
        df = df[df.COV_raw >= 8]
    
    dictCOV = {}
    for i in genoInfo.index:
        dictCOV[i] = df[(df.CHR == i)].copy()
        dictCOV[i]['PlotPos'] = dictCOV[i].POS + genoInfo.loc[i, 'prevLen']
    return dictCOV


def computeCOV(dfSNV, START, END, dataType):
    # compute mean of middle 50, while WGS median is an integer
    if 'POS' in dfSNV.columns: # dfSNV
        df = dfSNV[(dfSNV.POS>=START) & (dfSNV.POS<=END)]
        data = df.COV.values
    else: # dfCOV
        df = dfSNV[(dfSNV.START>=START) & (dfSNV.START<=END)]
        data = df.covMean.values
        
    data = data[data != 0] # remove na   
    covCount = len(data)
    # covered or not, such as gaps. chrY has very limited SNVs. Keep 0 for plotting coverage
    if dataType == 'WGS' and covCount <= 10:
        mm50, SEM_mm50 = 0, 0
    elif covCount <= 1:
        mm50, SEM_mm50 = 0, 0
    elif covCount <= 10: # WES/cfDNA may have limited data in a bin
        mm50 = np.mean(data)
        if mm50 == 0:
            SEM_mm50 = 0
        else:
            SEM_mm50 = np.std(data)/mm50
    else:
        data = np.sort(data)
        quantile = int(covCount/4)
        mm50 = (data[quantile:-quantile]).mean()
        SEM_mm50 = np.std(data[quantile:-quantile])/mm50
    return covCount, round(mm50, 5), round(SEM_mm50, 5)

def genoInfo2build(genoInfo):
    if genoInfo.loc['1', 'chrSize'] == 249250621:
        return 'hg19'
    elif genoInfo.loc['1', 'chrSize'] == 248956422:
        return 'hg38'

def coverage_chr(dictCOV, dictSNV, CHR, genoInfo, gender, covWindow, dataType):
    build = genoInfo2build(genoInfo)
    dfCOV = dictCOV[CHR]
    dfSNV = dictSNV[CHR]
    # sliding window
    chrSize = genoInfo.loc[CHR,'chrSize']
    # no data, female chrY, male with chrY_loss
    if len(dfCOV) <= 10:
        return pd.DataFrame()
    # seperate chrX in male
    elif gender == 'M' and CHR == 'X':
        [PAR1_S,PAR1_E] = PAR_Xpos(build)[0]
        [PAR2_S,PAR2_E] = PAR_Xpos(build)[1]
        windowPAR1 = np.arange(0, PAR1_E, step=covWindow)
        windowX    = np.arange(PAR1_E, PAR2_S, step=covWindow)
        windowPAR2 = np.arange(PAR2_S, PAR2_E, step=covWindow)
        window = np.concatenate((windowPAR1, windowX, windowPAR2, chrSize), axis=None)
    # one bin for chrY
    elif gender == 'M' and CHR == 'Y':
        window = np.array([dfCOV.iloc[0,1],dfCOV.iloc[-1,1]])
    else:
        window = np.arange(0, chrSize, step=covWindow)
        window = np.concatenate((window, chrSize), axis=None)

    #compute coverage
    df = pd.DataFrame([window[:-1]+1, window[1:]], index=['START','END']).T
    df['CHR'] = CHR
    for i in df.index:
        START, END = df.loc[i,'START'], df.loc[i,'END']
        #print(START, END)
        covCount, covMean, covSEM = computeCOV(dfCOV, START, END, dataType)
        snpCount = ((dfSNV.POS>=START) & (dfSNV.POS<=END)).sum()
        
        df.loc[i,'covCount'], df.loc[i,'covMean'], df.loc[i,'covSEM'] = covCount, covMean, covSEM
        df.loc[i,'snpCount'] = snpCount
    
    # covSmooth
    if CHR == 'Y': # only 1 bin
        df['covSmooth'] = df['covMean']
    else:
        df['covSmooth'] = df['covMean'].rolling(3,center=True).median()       
    # PAR2 should not be smoothed
    if gender == 'M' and CHR == 'X':
        df.loc[df.index[-1], 'covSmooth'] = df.loc[df.index[-1], 'covMean']
        
    df['plotSTART'] = df.START + genoInfo.loc[CHR, 'prevLen']
    df['plotEND'] = df.END + genoInfo.loc[CHR, 'prevLen']
    return df

def coverage_genome(dictCOV, dictSNV, genoInfo, gender, dataType, chrY_status, covWindow=2e6):
    print("********** Working on coverage bins ... **********")
    dfCOV = pd.DataFrame()
    for CHR in dictCOV:
        if CHR == 'Y' and chrY_status in ['no_Y']:
            print("Skipping chrY, because of", chrY_status)
        else:
            #print("Preparing chr" + CHR + " ...")
            dfCOV = dfCOV.append(coverage_chr(dictCOV, dictSNV, CHR, genoInfo, gender, covWindow, dataType))
    dfCOV.reset_index(inplace=True)
    dfCOV.drop('index',axis=1,inplace=True)
    return dfCOV

def PAR_Xpos(build):
    # start and end position of PAR1 and PAR2
    if build=='hg19':
        return [60001, 2699520], [154931044, 155260560]
    elif build=='hg38':
        return [10001, 2781479], [155701383, 156030895]
    
################################################################
# Circular Binary Segmentation (CBS)
################################################################
def cbs_stat(x):
    '''Given x, Compute the subinterval x[i0:i1] with the maximal segmentation statistic t. 
    Returns t, i0, i1'''
    
    x0 = x - np.mean(x)
    n = len(x0)
    y = np.cumsum(x0)
    e0, e1 = np.argmin(y), np.argmax(y)
    i0, i1 = min(e0, e1), max(e0, e1)
    s0, s1 = y[i0], y[i1]
    return (s1-s0)**2*n/(i1-i0+1)/(n+1-i1+i0), i0, i1+1

def tstat(x, i):
    '''Return the segmentation statistic t testing if i is a (one-sided)  breakpoint in x'''
    n = len(x)
    s0 = np.mean(x[:i])
    s1 = np.mean(x[i:])
    return (n-i)*i/n*(s0-s1)**2

def cbs(x, shuffles=1000, p=.05):
    '''Given x, find the interval x[i0:i1] with maximal segmentation statistic t. Test that statistic against
    given (shuffles) number of random permutations with significance p.  Return True/False, t, i0, i1; True if
    interval is significant, false otherwise.'''

    max_t, max_start, max_end = cbs_stat(x)
    if max_end-max_start == len(x):
        return False, max_t, max_start, max_end
    if max_start < 5:
        max_start = 0
    if len(x)-max_end < 5:
        max_end = len(x)
    thresh_count = 0
    alpha = shuffles*p
    xt = x.copy()
    for i in range(shuffles):
        np.random.shuffle(xt)
        threshold, s0, e0 = cbs_stat(xt)
        if threshold >= max_t:
            thresh_count += 1
        if thresh_count > alpha:
            return False, max_t, max_start, max_end
    return True, max_t, max_start, max_end

def rsegment(x, start, end, L=[], shuffles=1000, p=.05):
    '''Recursively segment the interval x[start:end] returning a list L of pairs (i,j) where each (i,j) is a significant segment.
    '''
    threshold, t, s, e = cbs(x[start:end], shuffles=shuffles, p=p)
    #log.info('Proposed partition of {} to {} from {} to {} with t value {} is {}'.format(start, end, start+s, start+e, t, threshold))
    if (not threshold) | (e-s < 5) | (e-s == end-start):
        L.append((start, end))
    else:
        if s > 0:
            rsegment(x, start, start+s, L)
        if e-s > 0:
            rsegment(x, start+s, start+e, L)
        if start+e < end:
            rsegment(x, start+e, end, L)
    return L

def segment(x, shuffles=1000, p=.05):
    '''Segment the array x, using significance test based on shuffles rearrangements and significance level p
    '''
    start = 0
    end = len(x)
    L = []
    rsegment(x, start, end, L, shuffles=shuffles, p=p)
    return L

def validate(x, L, shuffles=1000, p=.01):
    S = [x[0] for x in L]+[len(x)]
    SV = [0]
    left = 0
    for test, s in enumerate(S[1:-1]):
        t = tstat(x[S[left]:S[test+2]], S[test+1]-S[left])
        #log.info('Testing validity of {} in interval from {} to {} yields statistic {}'.format(S[test+1], S[left], S[test+2], t))
        threshold = 0
        thresh_count = 0
        site = S[test+1]-S[left]
        xt = x[S[left]:S[test+2]].copy()
        flag = True
        for k in range(shuffles):
            np.random.shuffle(xt)
            threshold = tstat(xt, site)
            if threshold > t:
                thresh_count += 1
            if thresh_count >= p*shuffles:
                flag = False
                #log.info('Breakpoint {} rejected'.format(S[test+1]))
                break
        if flag:
            #log.info('Breakpoint {} accepted'.format(S[test+1]))
            SV.append(S[test+1])
            left += 1
    SV.append(S[-1])
    return SV

################################################################
# segmentation function
################################################################
def segmentation(dictSNV, CHR, genoInfo, workDir, START=0, END=1e9, dataCol='VAF', dataType='WGS'):
    #print('***** Working on '  + dataCol + ' segmentation for chr' + CHR + ' *****')
    dfData = dictSNV[str(CHR)]
    dfData = dfData[(dfData.POS>=START) & (dfData.POS<=END)].copy()
    if len(dfData) <= 10:
        print('Too few data points (<=10) for chr' + CHR + ' ...')
        return pd.DataFrame(), [], pd.DataFrame()
    
    if dataCol == 'VAF':
        vafMeanCHR = genoInfo.loc[CHR, 'vafMeanCHR']
        # use both upper and lower half to generate accurate breakpoints
        upper = dfData[dfData.VAF>vafMeanCHR].index
        lower = dfData[dfData.VAF<=vafMeanCHR].index
        dfData.loc[upper, 'MAF'] = 2*vafMeanCHR - dfData.loc[upper, 'VAF']
        dfData.loc[lower, 'MAF'] = dfData.loc[lower, 'VAF']
        # remove homozygotes
        dfData = dfData[dfData.MAF > 0].copy()
        dfData.reset_index(inplace=True)
        dfData.drop(['index'],axis=1,inplace=True)
        # no heterozygotes
        if len(dfData) <= 10:
            print('Too few heterozygotes (<=10) for chr' + CHR + ' ...')
            return pd.DataFrame(), [], pd.DataFrame()
        data = dfData.MAF.values
        Segments_folder = workDir + '/Segments_VAF'
        
    elif dataCol == 'ROH':
        # %5 of tolerance
        if dataType == 'WGS':
            dfData['MAF'] = dfData['VAF'].rolling(100,center=True).apply(\
                lambda x: (0.5 if ((x!=1) & (x!=0)).sum() >= 5 else 1), raw=True)
        else:
            dfData['MAF'] = dfData['VAF'].rolling(20,center=True).apply(\
                lambda x: (0.5 if ((x!=1) & (x!=0)).sum() >= 2 else 1), raw=True)
        dfData = dfData[dfData.MAF > 0].copy()
        dfData.reset_index(inplace=True)
        dfData.drop(['index'],axis=1,inplace=True)
        data = dfData.MAF.values
        Segments_folder = workDir + '/Segments_ROH'
        
    elif dataCol == 'COV':
        dfData = dfData[dfData.COV > 0].copy()
        dfData.COV = np.log2(dfData.COV)
        dfData.reset_index(inplace=True)
        dfData.drop(['index'],axis=1,inplace=True)
        data = dfData.COV.values 
        Segments_folder = workDir + '/Segments_COV'
    
    # load pre-generated segments
    if os.path.exists(Segments_folder + '/chr' + CHR):
        if CHR == '1':
            print('Loading segments ...')
        SV = pd.read_csv(Segments_folder + '/chr' + CHR, sep='\t', dtype={'BreakPoint':int})
        SV = SV.BreakPoint.values
    else:
        if os.path.exists(Segments_folder) == False:
            os.mkdir(Segments_folder)
        print('Generating segments for chr' + CHR + ' ...')
        L = segment(data)
        print('Validating segments for chr' + CHR + ' ...')
        # runs may return slightly different segments
        SV = validate(data, L)
        
        # save segments
        with open(Segments_folder + '/chr' + CHR, 'w') as out:
            out.write('BreakPoint\n')
            for i in SV:
                out.write(str(i) + '\n')

    # retrieve segment positions for visualization
    dfSV = pd.DataFrame([SV[:-1], SV[1:]], index=['BreakPt1','BreakPt2']).T       
    dfPOS = dfData[['POS']]
    # discard break points, because it is ambiguous which event it belongs to
    dfSV['BP1new'] = dfSV.BreakPt1 + 1
    dfSV['BP2new'] = dfSV.BreakPt2 - 1
    dfSV = pd.merge(dfSV, dfPOS, left_on='BP1new', right_index=True, how='left')
    dfSV = pd.merge(dfSV, dfPOS, left_on='BP2new', right_index=True, how='left')
    dfSV['SIZE'] = dfSV.POS_y - dfSV.POS_x
    # compute metadata, such as mafMean
    for i in dfSV.index:
        POS1 = dfSV.loc[i,'POS_x']
        POS2 = dfSV.loc[i,'POS_y']
        # median is more robust mean regarding outliers
        dfSV.loc[i,'MEAN'] = dfData[(dfData.POS>=POS1)&(dfData.POS<=POS2)].COV.mean()
        dfSV.loc[i,'SEM'] = dfData[(dfData.POS>=POS1)&(dfData.POS<=POS2)].COV.std()/dfSV.loc[i,'MEAN']
        dfSV.loc[i,'MEDIAN'] = dfData[(dfData.POS>=POS1)&(dfData.POS<=POS2)].COV.median()
    
    # keep the entire chrX, instead of PAR1 for male
    if dataCol == 'VAF' and CHR == 'X':
        dfData = dictSNV[str(CHR)].copy()
        vafMeanCHR = genoInfo.loc[CHR, 'vafMeanCHR']
        # use both upper and lower half to generate accurate breakpoints
        upper = dfData[dfData.VAF>vafMeanCHR].index
        lower = dfData[dfData.VAF<=vafMeanCHR].index
        dfData.loc[upper, 'MAF'] = 2*vafMeanCHR - dfData.loc[upper, 'VAF']
        dfData.loc[lower, 'MAF'] = dfData.loc[lower, 'VAF']
        # remove homozygotes
        dfData = dfData[dfData.MAF > 0].copy()
        dfData.reset_index(inplace=True)
        dfData.drop(['index'],axis=1,inplace=True)
    
    return dfData, SV, dfSV

################################################################
# VAF segmentation for heterozygotes
################################################################
def segment_genome(dictSNV, gender, chrY_status, genoInfo, workDir, dataCol='VAF', dataType='WGS'):
    print('*********', dataCol, 'segmentation ... *********')
    dictMAF, dictSV, dictDfSV  = {},{},{}
    build = genoInfo2build(genoInfo)
    
    for CHR in dictSNV:
        #print('Working on segments for chr' + CHR + ' ...')
        # no need to segment vaf for chrY
        if CHR == 'Y':
            #print('***** Skipping chrY *****')
            dfMAF, SV, dfSV = pd.DataFrame(), np.nan, pd.DataFrame()
        elif CHR == 'X' and gender == 'M':
            # keep PAR1, while PAR2 is too small
            [PAR1_S,PAR1_E] = PAR_Xpos(build)[0]
            dfMAF, SV, dfSV = segmentation(dictSNV, CHR, genoInfo, workDir, PAR1_S, PAR1_E, dataCol, dataType)  
        else:
            dfMAF, SV, dfSV = segmentation(dictSNV, CHR, genoInfo, workDir, dataCol=dataCol, dataType=dataType)
        dictMAF[CHR], dictSV[CHR], dictDfSV[CHR] = dfMAF, SV, dfSV
    return dictMAF, dictSV, dictDfSV

def segment_smooth(dfSNV, dfMAF, SV, dataType, minSize=1e5, minHet=30, dataCol='VAF'):
    # WGS: median of snpCount_1Mb = 1136, hetCount_1Mb = 683, hetRatio = 0.6   
    # To have enough heterozygotes for model fitting, min(minSize)=1e4, min(minHET)=10
    # Due to off-target SNPs, minHET = minBins should be OK
    minSize = np.max([1e4, minSize])
    minHet  = np.max([10, minHet])
    # VAFcutoff = 0.01 indicates 0.02 of cnLOH, 0.028 of Loss, or 0.062 of Gain
    # VAFcutoff = 0.01 will generate a lot more sub segments
    vafCutoff = 0.01
    
    # no segments due to no heterozygous SNP markers
    if len(SV) == 0:
        return pd.DataFrame(), [], pd.DataFrame()
    
    smoothed = False
    round_num = 0
    dfSmooth = pd.DataFrame()
    while smoothed == False:
        smoothed = True
        round_num += 1
        
        # retrieve segment positions
        df = pd.DataFrame([SV[:-1], SV[1:]], index=['BreakPt1','BreakPt2']).T       
        dfPOS = dfMAF[['POS']]
        # discard break points, because it is ambiguous which event it belongs to
        df['BP1new'] = df.BreakPt1 + 1
        df['BP2new'] = df.BreakPt2 - 1
        df = pd.merge(df, dfPOS, left_on='BP1new', right_index=True, how='left')
        df = pd.merge(df, dfPOS, left_on='BP2new', right_index=True, how='left')
        
        df['Round'] = round_num
        df['mafCount'] = df.BreakPt2 - df.BreakPt1 - 1 # exlude both ends
        df['minMCount'] = df.mafCount.rolling(window=2).min() # min maf/het count
        df['SIZE'] = df.POS_y - df.POS_x # event size
        df['minSIZE'] = df.SIZE.rolling(window=2).min() # min event size
        
        # compute metadata, such as mafMean
        for i in df.index:
            POS1 = df.loc[i,'POS_x']
            POS2 = df.loc[i,'POS_y']
            # median is more robust mean regarding outliers
            df.loc[i,'mafMean'] = dfMAF[(dfMAF.POS>=POS1)&(dfMAF.POS<=POS2)].MAF.median()
            df.loc[i,'mafSTD'] = dfMAF[(dfMAF.POS>=POS1)&(dfMAF.POS<=POS2)].MAF.std()
            df.loc[i,'snpCount'] = ((dfSNV.POS >= POS1) & (dfSNV.POS <= POS2)).sum()
            df.loc[i,'hetCount'] = ((dfMAF.POS >= POS1) & (dfMAF.POS <= POS2) & (dfMAF.MAF>0)).sum()
            
        df['meanDiff'], df['DROP'] = np.abs(df.mafMean.diff()), np.nan
      
        if round_num == 1:
            maxSTD = df.mafSTD.mean() + 2*df.mafSTD.std()
        
        # one segment left or to begin with, such as germline samples
        if len(df)==1:
            df = df[['BreakPt1','BreakPt2','POS_x', 'POS_y','SIZE','snpCount', \
                'mafCount', 'hetCount', 'mafMean','mafSTD','meanDiff']]
            return dfSmooth, SV, df
        
        for i in df.index:
            # vafCutoff for each row
            vafCutoffRow = size2cutoff(vafCutoff, dataType, df.loc[i,'minSIZE'], df.loc[i,'minMCount'])
            
            # if dropped a previous row, need to skip this row
            if i != 0 and df.loc[i-1,'DROP'] == df.loc[i-1,'DROP']:
                continue 
            
            # drop small segments FIRST
            if df.loc[i,'mafCount']<minHet or df.loc[i,'SIZE']<minSize:
                smoothed = False
                df.loc[i,'DROP'] = 'Small'
                df.loc[i,'dropBP'] = drop_BP(df, i)
            # similar mafMean with the previous segment
            elif df.loc[i,'meanDiff'] < vafCutoffRow:
                smoothed = False
                df.loc[i,'DROP'] = 'mafMean'
                df.loc[i,'dropBP'] = df.loc[i,'BreakPt1']        
            # segments with excess of heterozygosity
            elif df.loc[i,'mafCount'] > df.loc[i,'snpCount']*0.9 and dataCol=='VAF':
                smoothed = False
                df.loc[i,'DROP'] = 'EOH'
                df.loc[i,'dropBP'] = drop_BP(df, i)
            # large mafSTD indicates poor regions, but mafSTD of WES/cfDNA can be very big
            elif dataType == 'WGS' and df.loc[i,'mafSTD'] > maxSTD and dataCol=='VAF' and df.loc[i,'SIZE']<minSize*10:
                smoothed = False
                df.loc[i,'DROP'] = 'mafSTD'
                df.loc[i,'dropBP'] = drop_BP(df, i)
            # too few SNP in WGS, probably a poor sequencing gap such as centromere
            elif dataType=='WGS' and df.loc[i,'snpCount']/df.loc[i,'SIZE'] < 0.0002:
                smoothed = False
                df.loc[i,'DROP'] = 'Gap'
                df.loc[i,'dropBP'] = drop_BP(df, i)
                
        # drop breakpoints from SV after going over df
        if smoothed == False:
            SV = [x for x in SV if x not in df.dropBP.values]
        dfSmooth = pd.concat([dfSmooth, df], sort=True, ignore_index=True)

    # to avoid return SV below
    SVs = SV
    dfSVs = df[['BreakPt1','BreakPt2','POS_x', 'POS_y','SIZE','snpCount', \
                'mafCount', 'hetCount', 'mafMean','mafSTD','meanDiff']]

    return dfSmooth, SVs, dfSVs

def size2cutoff(vafCutoff, dataType, eventSize, mafCount):
    # small segments are more likely to be merged
    if dataType == 'WGS':
        if eventSize > 1e7:     return vafCutoff
        elif eventSize > 1e6:   return vafCutoff + 0.02
        elif eventSize > 1e5:   return vafCutoff + 0.03
        else: return vafCutoff + 0.04
    else:
        if mafCount > 400: return vafCutoff
        if mafCount > 200: return vafCutoff + 0.01
        elif mafCount > 100: return vafCutoff + 0.02
        elif mafCount > 50: return vafCutoff + 0.03
        else: return vafCutoff + 0.04

def drop_BP(df, i, dataCol = 'mafMean'):
    BreakPt1 = df.loc[i,'BreakPt1']
    BreakPt2 = df.loc[i,'BreakPt2']
    # first segment
    if i == 0:
        return BreakPt2
    # last segment
    elif i == df.index[-1]:
        return BreakPt1
    # combine it into the larger neighbor to minimize its side effects
    else:
        value1 = df.loc[i-1,dataCol]
        value  = df.loc[i,dataCol]
        value2 = df.loc[i+1,dataCol]
        if value <= value1 and value <= value2:
            if value1<value2:
                return BreakPt1
            else:
                return BreakPt2
        elif value >= value1 and value >= value2:
            if value1<value2:
                return BreakPt2
            else:
                return BreakPt1
        else:
            if abs(value - value1) < abs(value - value2):
                return BreakPt1
            else: 
                return BreakPt2

def smooth_genome(dictSNV, dictMAF, dictSV, gender, dataType, minSize, minHet, build, dataCol='VAF', PRINT=False):
    print('********* ' + dataCol + ' segments smoothing ... *********')
    dictSmooth, dictSVs, dictDfSVs = {},{},{}
    for CHR in dictSNV:
        if PRINT == True:
            print('***** VAF segment smoothing for chr' + CHR + ' *****')
        dfSNV, dfMAF, SV = dictSNV[CHR], dictMAF[CHR], dictSV[CHR]
        
        if CHR != 'Y': # no Y segments
            dfSmooth, SVs, dfSVs = segment_smooth(dfSNV, dfMAF, SV, dataType, minSize, minHet, dataCol=dataCol)
        
        # handle X and Y separately
        # no need to smooth Y in female, or no SNP markers
        if (gender == 'F' and CHR == 'Y') or len(dfMAF) == 0:
            dfSmooth, SVs, dfSVs = pd.DataFrame(), [], pd.DataFrame()
        elif gender == 'M' and CHR == 'Y':
            dfSmooth, SVs = np.nan, np.nan
            dfSVs = pd.DataFrame([[0, 3e7]], columns=['POS_x','POS_y'])
            dfSVs['snpCount'] = len(dfSNV)
            dfSVs['hetCount'] = np.nan
            dfSVs['mafMean'] = 0
            dfSVs['mafSTD'] = np.nan
        elif gender == 'M' and CHR == 'X':  
            POS_xy = []
            [PAR1_S,PAR1_E] = PAR_Xpos(build)[0]
            [PAR2_S,PAR2_E] = PAR_Xpos(build)[1]
            for i in dfSVs.index:
                if dfSVs.loc[i,'POS_y'] < PAR1_E:
                    POS_xy.append([dfSVs.loc[i,'POS_x'], dfSVs.loc[i,'POS_y']])
                elif dfSVs.loc[i,'POS_x'] < PAR1_E:
                    POS_xy.append([dfSVs.loc[i,'POS_x'], PAR1_E])
            POS_xy.append([PAR1_E, PAR2_S])
            POS_xy.append([PAR2_S, PAR2_E])
                
            # replace dfSVs
            dfSVs = pd.DataFrame(POS_xy, columns=['POS_x','POS_y'])
            dfSVs['BreakPt1'] = 'np.nan'
            dfSVs['BreakPt2'] = 'np.nan'

            for i in dfSVs.index:
                POS1 = dfSVs.loc[i,'POS_x']
                POS2 = dfSVs.loc[i,'POS_y']
                dfSVs.loc[i,'snpCount'] = ((dfSNV.POS >= POS1) & (dfSNV.POS <= POS2)).sum()
                dfSVs.loc[i,'hetCount'] = ((dfMAF.POS >= POS1) & (dfMAF.POS <= POS2)).sum()
                MAFs = dfMAF[(dfMAF.POS >= POS1) & (dfMAF.POS <= POS2)].MAF.values
                if len(MAFs) >= 3:
                    dfSVs.loc[i,'mafMean'] = np.mean(MAFs)
                    dfSVs.loc[i,'mafSTD'] = np.std(MAFs)
                else: # For example, WES has no heterozygotes
                    dfSVs.loc[i,'mafMean'], dfSVs.loc[i,'mafSTD'] = np.nan, np.nan
            # samll sub-segments
            #dfSVs = dfSVs[dfSVs.covCount > 10]
            
        dictSmooth[CHR] = dfSmooth
        dictSVs[CHR] = SVs
        dictDfSVs[CHR] = dfSVs
        
    return dictSmooth, dictSVs, dictDfSVs

################################################################
# Merge VAF and ROH segments
################################################################
def merge_VAF_ROH(dictSNV, dictMAF, dictDfSVs, dictDfSVs_ROH, CHR, minSize, minBins):
    dfSNV = dictSNV[CHR]
    dfMAF = dictMAF[CHR]
    dfSVs = dictDfSVs[CHR]
    dfSVs_ROH = dictDfSVs_ROH[CHR]
    
    # no VAF segments, such as chrY
    if dfSVs.empty:
        return dfSVs

    # dfROH
    dfROH = dfSVs_ROH[dfSVs_ROH.mafMean == 1].copy() 
    
    # no ROH events
    if dfROH.empty:
        return dfSVs
    
    print('Found ROH event(s) on chr' + CHR + '...')
    
    BP = list(dfSVs.POS_x.values)+list(dfSVs.POS_y.values) + list(dfROH.POS_x.values) + list(dfROH.POS_y.values)
    for i in dfSVs.index:
        START1, END1 = dfSVs.loc[i, 'POS_x'], dfSVs.loc[i, 'POS_y']
        for j in dfROH.index:
            START2, END2 = dfROH.loc[j, 'POS_x'], dfROH.loc[j, 'POS_y']
            # With intersection, remove VAF breakpoint(s)
            if START2 <= START1 <= END2:
                BP.remove(START1)
            if START2 <= END1 <= END2:
                BP.remove(END1)
    SVs = np.sort(BP)
    df = pd.DataFrame([SVs[:-1], SVs[1:]], index=['POS_x','POS_y']).T
    df['SIZE'] = df.POS_y - df.POS_x
    
    for i in df.index:
        POS1, POS2 = df.loc[i, 'POS_x'], df.loc[i, 'POS_y']
        df.loc[i,'snpCount'] = ((dfSNV.POS >= POS1) & (dfSNV.POS <= POS2)).sum()  
        df.loc[i,'mafCount'] = ((dfMAF.POS >= POS1) & (dfMAF.POS <= POS2)).sum() 
        df.loc[i,'mafMean'] = dfMAF[(dfMAF.POS>=POS1)&(dfMAF.POS<=POS2)].MAF.mean()
        df.loc[i,'mafSTD'] = dfMAF[(dfMAF.POS>=POS1)&(dfMAF.POS<=POS2)].MAF.std()
    # drop breakpoints and gaps
    df = df[(df.SIZE > minSize) & (df.snpCount > minBins)]
    
    return df

def merge_genome(dictSNV, dictMAF, dictDfSVs, dictDfSVs_ROH, minSize, minBins):
    print('********* Merging VAF and ROH segments ... *********')
    dictDfSVs_Merged = {}
    for CHR in dictDfSVs:
        #print('Working on chr' + CHR + ' ...')
        dictDfSVs_Merged[CHR] = merge_VAF_ROH(dictSNV, dictMAF, dictDfSVs, dictDfSVs_ROH, CHR, minSize, minBins)
    return dictDfSVs_Merged

################################################################
# SubChrom algorithm
################################################################

# model for one normal distribution
def OneNorm(x,mean,sd):
    w = sts.norm.cdf(x,mean,sd)
    return w

# model for two normal distributions
# 'mean' is the mean value of the distribution
# 'frac' is the fraction of 1 distribution, as the two separate distributions are not necessarty the same
# 'vafMean' is the mean value of all vaf data
def TwoNorm(x,mean,sd,frac,vafMean):
    w = frac*sts.norm.cdf(x,mean,sd) + (1-frac)*sts.norm.cdf(x,vafMean*2-mean,sd)
    return w

# fit data into two models
def fitModel(dfSNV, START, END, mafMean, mafSTD, vafMeanGenome):
    # vafCount, vafMean, oneMean, oneSTD, oneProb, twoMeanUp, twoMeanDown, twoSTD, twoProb, vafEvent
    df = dfSNV[(dfSNV.POS>=START) & (dfSNV.POS<=END)].copy()
    snpCount = len(df)
    df = df[(df.VAF > 0) & (df.VAF < 1)] # remove homozygotes
    vafCount = len(df)
    vafMean = df.VAF.mean()
    
    # too few heterozygotes
    if vafCount <= snpCount/25:
        return [vafCount, 0.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 'ROH']
    
    # clearly TwoNorm, fit into OneNorm model, which returns accurate vafMean
    if mafMean < 0.2:
        if (df.VAF > 0.5).sum() >= (df.VAF < 0.5).sum():
            df = df[df.VAF > 0.5]
        else:
            df = df[df.VAF < 0.5]
        try:    
            popt1, pcov1 = opt.curve_fit(OneNorm,np.sort(df.VAF.values),
                               np.linspace(0,1,len(df)),p0=[0.9, 0.1])
            twoMean = round(popt1[0],6)
            twoMeanUp = np.max([twoMean, 1-twoMean])
            twoMeanDown = np.min([twoMean, 1-twoMean])
            twoSTD = round(popt1[1],6)
            return [vafCount, vafMean, np.nan, np.nan, 0, \
                    twoMeanUp, twoMeanDown, twoSTD, 1, 'TwoNorm']
        except:
            return [vafCount, 0.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 'ROH']
    
    # detecting subclonal events. NOT working well for cfDNA WGS, because VAF is dispearsed
    #if mafMean > 0.43:
        # remove lowMAF SNV to improve model fitting
    #    df = df[(df.VAF >= 0.2) & (df.VAF <= 0.8)]

    # OneNorm
    try:
        popt1, pcov1 = opt.curve_fit(OneNorm, np.sort(df.VAF.values), 
                                     np.linspace(0,1,len(df)),p0=[0.5, 0.1])
        oneMean = round(popt1[0],6)
        oneSTD = round(popt1[1],6)
        oneProb = sts.kstest(df.VAF.values, OneNorm, args=popt1)[1]
    except:
        oneMean,oneSTD,oneProb = vafMeanGenome, np.nan, -1
    
    # unbalanced VAF distribution, oneMean can't be > 0.5
    if oneMean > 0.501 or oneMean < vafMeanGenome/1.05: # add some noise
        print('Unbalanced VAF oneMean on chr'+dfSNV.iloc[0,0]+':', oneMean, "vafMeanGenome:", vafMeanGenome)
        try:
            popt1, pcov1 = opt.curve_fit(lambda x, sd: OneNorm(x, vafMeanGenome, sd), np.sort(df.VAF.values), 
                                             np.linspace(0,1,len(df)),p0=[0.1])
            oneMean = vafMeanGenome
            oneSTD = round(popt1[0],6)
            oneProb = sts.kstest(df.VAF.values, OneNorm, args=[vafMeanGenome, popt1])[1]
        except:
            oneMean,oneSTD,oneProb = vafMeanGenome, np.nan, -1
    
    # TwoNorm
    try:
        popt2, pcov2 = opt.curve_fit(lambda x,mean,sd,frac: TwoNorm(x,mean,sd,frac,oneMean),
                                    np.sort(df.VAF.values),np.linspace(0,1,len(df)),p0=[0.3,0.1,0.5])
        twoMean = round(popt2[0],6)
        if twoMean > oneMean:
            twoMeanUp = twoMean
            twoMeanDown = oneMean * 2 - twoMeanUp
        else:
            twoMeanDown = twoMean
            twoMeanUp = oneMean * 2 - twoMeanDown
        twoSTD = round(popt2[1],6)
        twoProb = sts.kstest(df.VAF.values, TwoNorm, args=np.append(popt2, oneMean))[1]
    except:
        # small failed events: very likely not normal distribution at all, and thus poor regions
        # large failed events: very likely fit into one normal distribution only
        return [vafCount, vafMean, oneMean, oneSTD, oneProb, \
                np.nan, np.nan, np.nan, -1, 'OneNorm']
  
    # compare two models
    # TwoNorm fits better, adding some noise
    if oneProb <= twoProb/1.05 and oneSTD > twoSTD:
        vafEvent = 'TwoNorm'
    # both fit similar & TwoNorm STD is smaller
    elif oneProb/1.1 <= twoProb and oneSTD/1.1 > twoSTD:
        vafEvent = 'TwoNorm'

    else:
        vafEvent = 'OneNorm'
    return [vafCount, vafMean, oneMean, oneSTD, oneProb, twoMeanUp, twoMeanDown, twoSTD, twoProb, vafEvent]  
    
def subchrom_chr(dictSNV, dictCOV, dfCOV, dictDfSVs, CHR, vafMeanGenome, dataType):
    # working dataframe of chr
    CHR = str(CHR)
    dfCHR = dictSNV[CHR]
    dfSVs = dictDfSVs[CHR]
    dfCOVCHR = dictCOV[CHR]
    dfCOVbin = dfCOV[dfCOV.CHR == CHR]
    EVENTS = []
    for i in dfSVs.index:
        START, END = dfSVs.loc[i, 'POS_x'], dfSVs.loc[i, 'POS_y']
        mafMean, mafSTD = dfSVs.loc[i, 'mafMean'], dfSVs.loc[i, 'mafSTD']
        # compute coverage
        snpCount, covMean, covSEM = computeCOV(dfCHR, START, END, dataType)
        # to avoid local sequencing, such as chr14q
        if END-START > 2e7:
            covCountCOV1, covMean1, covSEM1 = computeCOV(dfCOVbin, START, END, dataType)
            covCountCOV2, covMean2, covSEM2 = computeCOV(dfCOVCHR, START, END, dataType)
            if covSEM1 > covSEM2:
                covCountCOV, covMean, covSEM = covCountCOV2, covMean2, covSEM2
            else:
                covCountCOV, covMean, covSEM = covCountCOV1, covMean1, covSEM1
        else: 
            covCountCOV, covMean, covSEM = computeCOV(dfCOVCHR, START, END, dataType)
        # compute vaf
        if CHR != 'Y':
            vafEvent = fitModel(dfCHR, START, END, mafMean, mafSTD, vafMeanGenome)
        else:
            vafEvent = [np.nan, 0.5, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 'Y']
        EVENTS.append([CHR, START, END, END-START, snpCount, covMean, covSEM, mafMean, mafSTD] + vafEvent) 
    df = pd.DataFrame(EVENTS, 
            columns=['CHR', 'START', 'END', 'SIZE', 'snpCount', 'covMean', 'covSEM',
                     'mafMean', 'mafSTD',
                     'vafCount', 'vafMean', 'oneMean', 'oneSTD', 'oneProb', 
                     'twoMeanUp', 'twoMeanDown', 'twoSTD', 'twoProb', 'vafEvent'])
    return df

def subchrom_genome(dictSNV, dictCOV, dfCOV, dictDfSVs, genoInfo, chrY_status, dataType):
    print('********** Model fitting for VAF segments **********')
    final = pd.DataFrame()
    for CHR in dictSNV:
        #print('***** fitModel for chr' + CHR + ' *****')
        # no need to fit chrY
        if CHR=='Y' and chrY_status in ['no_Y']:
            continue
        vafMeanGenome = genoInfo.vafMeanCHR.median()
        df = subchrom_chr(dictSNV, dictCOV, dfCOV, dictDfSVs, CHR, vafMeanGenome, dataType)
        final = final.append(df)
    final.reset_index(inplace=True)
    final.drop('index',axis=1,inplace=True)
    return final

################################################################
# Diploid depth
################################################################
def depth(dictSNV, dataType, depthMode):
    if depthMode == 'auto': # consider all autosomes
        meanCovChr = []
        for CHR in range(1,23):
            dfSNV = dictSNV[str(CHR)]
            covCount, covMean, covSEM = computeCOV(dfSNV, 1, dfSNV.POS.max(), dataType)
            meanCovChr.append(covMean)
        return np.median(meanCovChr)
    else:
        CHR = str(depthMode) # consider a specific chromosome, such as 7
        dfSNV = dictSNV[str(CHR)]
        covCount, covMean, covSEM = computeCOV(dfSNV, 1, dfSNV.POS.max(), dataType)
        return covMean
    
def diploidDepth(dictCOV, vafSeg, dataType, dipDep):
    print("********** Optimizing diploid depth ... **********")
    print("Input diploid mode:", dipDep)
    
    # handle both str and int
    try:
        return int(dipDep)  # Try to convert to int
    except ValueError:
        if dipDep in ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",\
                      "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                      "chr21","chr22","chrX"]:
            dipDep = dipDep[3:]
        elif 'chr' in str(dipDep):
            print('Not supported chromosome:', dipDep)
            dipDep = 'auto'
        elif dipDep == 'auto':
            dipDep = 'auto'
        else:
            print('Not supported chromosome:', dipDep)
            dipDep = 'auto'

        if dipDep == 'auto':
            OneNormSeg = vafSeg[vafSeg.vafEvent == 'OneNorm']
            if len(OneNormSeg) >= 3:
                seqDepth = OneNormSeg.covMean.median()
            else:
                seqDepth = depth(dictCOV, dataType, 'auto')
        else:
            seqDepth = depth(dictCOV, dataType, dipDep)
        print('Final diploid mode:', dipDep)
        print('SeqDepth:', seqDepth)
        return seqDepth
    
################################################################
# Post process VAF events
################################################################
def coverageEvent(seqDepth, coverage, gender='F', CHR='1', minTF=0.05):
    if gender == 'M' and CHR in ['X', 'Y']:
        covFrac = round((coverage - 0.5 * seqDepth) / (0.5 * seqDepth),5)
        log2ratio = round(np.log2(coverage / (0.5 * seqDepth)), 5)
    else:
        covFrac = round((coverage - seqDepth) / (0.5 * seqDepth),5)
        log2ratio = round(np.log2(coverage / seqDepth), 5)
        
    if coverage == 0: # uncovered regions
        covEvent = 'NoCoverage'
        covFrac = 0
    elif covFrac <= -minTF:
        covEvent = 'Loss'
    elif covFrac >=  minTF:
        covEvent = 'Gain'
    else:
        covEvent = 'Diploid'        
    return covEvent, covFrac, log2ratio

# VAF to cell fraction conversion function        
def VAFtoFrac(vafMean, mode='Diploid'):
    lowerMean = np.min([vafMean, 1-vafMean])
    if mode=='Diploid':
        return 1 - lowerMean*2
    elif mode == 'Loss':
        return (lowerMean - 0.5) / (0.5*lowerMean - 0.5) # OneLoss
    elif mode == 'Gain':
        return 1/lowerMean -2 # OneGain

def VAFtoAB(seqDepth, covMean, vafMean):
    # x gain of chrA , y gain of chrB
    # seqDepth * 'x' + seqDepth * 'y' = 2*(covMean-seqDepth)
    # (1-vafMean) * 'x' - vafMean * 'y' = 2*vafMean - 1
    left  = [[seqDepth, seqDepth],[1-vafMean, -vafMean]]
    right = [2*(covMean-seqDepth), 2*vafMean-1]
    results = np.linalg.solve(left, right)
    chrA_freq = round(results[0], 3)
    chrB_freq = round(results[1], 3)
    return (chrA_freq, chrB_freq) if chrA_freq<chrB_freq else (chrB_freq, chrA_freq)

def vafEVENTS(vafSeg, dictSNV, genoInfo, seqDepth, gender, dataType, minTF):
    print('********** Post processing VAF events ... **********')
    # exclude no coverage regions, probably due to off-target SNPs 
    #df = vafSeg[((vafSeg.vafEvent == 'TwoNorm') | (vafSeg.vafEvent == 'ROH')) & (vafSeg.covMean != 0)].copy()
    df = vafSeg[(vafSeg.covMean != 0)].copy()
    
    build = genoInfo2build(genoInfo)
    [PAR1_S,PAR1_E], [PAR2_S,PAR2_E] = PAR_Xpos(build)

    if len(df) == 0:
        return pd.DataFrame([], columns=['CHR', 'START', 'END', 'SIZE', 'covMean',
       'covEvent', 'log2ratio', 'cnvEvent', 'vafEvent', 'vafFrac',
       'vaf_cov_Ratio', 'chrA_num', 'chrA_freq', 'chrB_num', 'chrB_freq'])
    
    df['chrA_num'], df['chrA_freq'], df['chrB_num'], df['chrB_freq'] = np.nan, np.nan, np.nan, np.nan
    
    for i in df.index:
        CHR = df.loc[i,'CHR']
        # segments are autosomes, PAR, XY
        if gender == 'M' and CHR=='X' and (df.loc[i, 'END']<=PAR1_E or df.loc[i, 'START']>=PAR2_S):
            covEvent, covFrac, log2ratio = coverageEvent(seqDepth, df.loc[i,'covMean'], gender, 'PAR', minTF)     
        else:
            covEvent, covFrac, log2ratio = coverageEvent(seqDepth, df.loc[i,'covMean'], gender, CHR, minTF)
        df.loc[i,'covEvent'], df.loc[i,'covFrac'], df.loc[i,'log2ratio'] = covEvent, covFrac, log2ratio
        
        # vaf evetns
        vafShift = 0.5 - (df.loc[i,'twoMeanUp'] + df.loc[i,'twoMeanDown'])/2
        vafMeanNEW = df.loc[i,'twoMeanDown'] + vafShift
        vafFrac = VAFtoFrac(vafMeanNEW, covEvent)
        df.loc[i,'cnvEvent'] = df.loc[i,'covEvent']
        df.loc[i,'vafFrac'] = vafFrac
        if covFrac != covFrac or covFrac == 0 or vafFrac != vafFrac: # not covered
            vaf_cov_Ratio = 1
        else:
            vaf_cov_Ratio = abs(vafFrac/covFrac)
        df.loc[i,'vaf_cov_Ratio'] = vaf_cov_Ratio
        
        # chrA, chrB, cnvEvent
        if df.loc[i, 'vafEvent'] == 'OneNorm':
            if covEvent == 'Diploid':
                continue
            elif covEvent == 'Gain':
                df.loc[i,'chrA_num'], df.loc[i,'chrA_freq'] = 1, np.nan
                df.loc[i,'chrB_num'], df.loc[i,'chrB_freq'] = 1 + covFrac, covFrac
            else: # Loss
                df.loc[i,'chrA_num'], df.loc[i,'chrA_freq'] = 1, np.nan 
                df.loc[i,'chrB_num'], df.loc[i,'chrB_freq'] = 1 + covFrac, covFrac
                
        elif df.loc[i, 'vafEvent'] == 'ROH':
            if gender == 'M' and df.loc[i,'CHR'] in ['X','Y']:
                df.loc[i, 'vafEvent'] = 'maleXY'
                if covEvent == 'Diploid':
                    df.loc[i,'cnvEvent'] = 'cnLOH'
                    df.loc[i,'chrA_num'], df.loc[i,'chrA_freq'] = 0, np.nan
                    df.loc[i,'chrB_num'], df.loc[i,'chrB_freq'] = 1, np.nan
                else:
                    df.loc[i,'cnvEvent'] = covEvent
                    df.loc[i,'chrA_num'], df.loc[i,'chrA_freq'] = 0, np.nan
                    df.loc[i,'chrB_num'], df.loc[i,'chrB_freq'] = 1+covFrac, abs(covFrac)
            elif df.loc[i, 'covMean'] < seqDepth*0.7:
                df.loc[i,'cnvEvent'] = 'Loss'
                df.loc[i,'vafFrac'] = 1
                df.loc[i,'chrA_num'], df.loc[i,'chrA_freq'] = 1, np.nan
                df.loc[i,'chrB_num'], df.loc[i,'chrB_freq'] = 0, 1
            else:
                df.loc[i,'cnvEvent'] = 'cnLOH'
                df.loc[i,'vafFrac'] = 1
                df.loc[i,'chrA_num'], df.loc[i,'chrA_freq'] = 0, np.nan
                df.loc[i,'chrB_num'], df.loc[i,'chrB_freq'] = 2, 1
                
        else: # TwoNorm
            if covEvent == 'Diploid': # cn_LOH
                df.loc[i,'cnvEvent'] = 'cnLOH'
                df.loc[i,'chrA_num'], df.loc[i,'chrA_freq'] = 1-vafFrac, -vafFrac
                df.loc[i,'chrB_num'], df.loc[i,'chrB_freq'] = 1+vafFrac, vafFrac
            elif covEvent == 'Gain':
                if vaf_cov_Ratio < 0.7 or vaf_cov_Ratio > 1.43:
                    chrA_gain, chrB_gain = VAFtoAB(seqDepth, df.loc[i,'covMean'], vafMeanNEW)
                    df.loc[i,'cnvEvent'] = 'complexGain'
                    df.loc[i,'chrA_num'], df.loc[i,'chrA_freq'] = 1+chrA_gain, chrA_gain
                    df.loc[i,'chrB_num'], df.loc[i,'chrB_freq'] = 1+chrB_gain, chrB_gain
                else:
                    df.loc[i,'cnvEvent'] = covEvent
                    df.loc[i,'chrA_num'], df.loc[i,'chrA_freq'] = 1, np.nan
                    df.loc[i,'chrB_num'], df.loc[i,'chrB_freq'] = 1 + vafFrac, vafFrac
            else: # Loss
                if vaf_cov_Ratio < 0.7 or vaf_cov_Ratio > 1.43 :
                    chrA_gain, chrB_gain = VAFtoAB(seqDepth, df.loc[i,'covMean'], vafMeanNEW)
                    df.loc[i,'cnvEvent'] = 'complexLoss'
                    df.loc[i,'chrA_num'], df.loc[i,'chrA_freq'] = 1+chrA_gain, chrA_gain
                    df.loc[i,'chrB_num'], df.loc[i,'chrB_freq'] = 1+chrB_gain, chrB_gain
                else:
                    df.loc[i,'cnvEvent'] = covEvent
                    df.loc[i,'chrA_num'], df.loc[i,'chrA_freq'] = 1, np.nan 
                    df.loc[i,'chrB_num'], df.loc[i,'chrB_freq'] = 1 - vafFrac, -vafFrac
    
    df = df[['CHR', 'START', 'END', 'SIZE', 'covMean',
       'covEvent', 'covFrac', 'log2ratio', 'cnvEvent', 'vafEvent', 'vafFrac',
       'vaf_cov_Ratio', 'chrA_num', 'chrA_freq', 'chrB_num', 'chrB_freq']]
    
    # remove lowFrac cnLOH
    LOWcnLOH = ((df.cnvEvent == 'cnLOH') & (df.chrB_freq < 0.08))
    df = df[~LOWcnLOH]
    
    # remove germline cnLOH
    cnLOH = df[(df.cnvEvent == 'cnLOH') & (df.vafEvent == 'ROH')]
    for i in cnLOH.index:
        CHR = cnLOH.loc[i, 'CHR']
        CHRevents = df[df.CHR == CHR]
        # cnLOH is not the first and the last event, most likely germline events
        if i != CHRevents.index[0] and i != CHRevents.index[-1]:
            df.drop(i, axis=0, inplace=True)

    TwoNorms = (df.vafEvent == 'TwoNorm')|(df.vafEvent == 'ROH')
    return df[TwoNorms], df[(~TwoNorms)&(df.covEvent != 'Diploid')]

################################################################
# coverage segmentation
################################################################
def segment_genome_COV(dictCOV, gender, genoInfo, workDir, dataType, PON, doCovSeg):
    dictMAF, dictSV, dictDfSV  = {},{},{}
    # build
    if genoInfo.loc['1', 'chrSize'] == 249250621:
        build = 'hg19'
    elif genoInfo.loc['1', 'chrSize'] == 248956422:
        build = 'hg38'
        
    # no segmentation
    if (dataType != 'WGS'and PON == 'N') or (doCovSeg == 'False' or doCovSeg == False):
        print('********* No COV segmentation ... *********')
        return dictMAF, dictSV, dictDfSV
    
    print('********* COV segmentation ... *********')
    for CHR in dictCOV:
        if gender == 'M' and CHR == 'X':
            [PAR1_S,PAR1_E] = PAR_Xpos(build)[0]
            [PAR2_S,PAR2_E] = PAR_Xpos(build)[1]
            # keep nonPAR
            dfMAF, SV, dfSV = segmentation(dictCOV, 'X', genoInfo, workDir, START=PAR1_E, END=PAR2_S, dataCol='COV')
        # no need to segment female chrY
        elif gender == 'F' and CHR == 'Y':
            #print('***** Skipping female chrY *****')
            dfMAF, SV, dfSV = pd.DataFrame(), np.nan, pd.DataFrame()
        else:
            dfMAF, SV, dfSV = segmentation(dictCOV, CHR, genoInfo, workDir, dataCol='COV') 
        dictMAF[CHR], dictSV[CHR], dictDfSV[CHR] = dfMAF, SV, dfSV
            
    return dictMAF, dictSV, dictDfSV

def size2cutoff_COV(dataType, eventSize, covCount):
    # small segments are more likely to be merged
    if dataType == 'WGS':
        if eventSize > 1e7:     return 0.05
        elif eventSize > 1e6:   return 0.15
        elif eventSize > 1e5:   return 0.2
        else: return 0.3
    else:
        if covCount > 400: return 0.05
        if covCount > 200: return 0.1
        elif covCount > 100: return 0.15
        elif covCount > 50: return 0.2
        else: return 0.3
    
def segment_smooth_COV(dfMAF, SV, dataType, maxSEM, minSize=1e5, minBins=30):
    # To have enough heterozygotes for model fitting, min(minSize)=1e4, min(minBins)=10
    # mycn amp may have 14 bins
    minSize = np.max([1e4, minSize])
    minBins  = np.max([10, minBins])
        
    smoothed = False
    round_num = 0
    dfSmooth = pd.DataFrame()
    while smoothed == False:
        smoothed = True
        round_num += 1
        
        # retrieve segment positions
        df = pd.DataFrame([SV[:-1], SV[1:]], index=['BreakPt1','BreakPt2']).T       
        dfPOS = dfMAF[['POS']] 
        # discard break points, because it is ambiguous which event it belongs to
        df['BP1new'] = df.BreakPt1 + 1
        df['BP2new'] = df.BreakPt2 - 1
        df = pd.merge(df, dfPOS, left_on='BP1new', right_index=True, how='left')
        df = pd.merge(df, dfPOS, left_on='BP2new', right_index=True, how='left')
        
        df['Round'] = round_num
        df['covCount'] = df.BreakPt2 - df.BreakPt1 - 1 # exlude both ends
        df['minCovCount'] = df.covCount.rolling(window=2).min() # min maf/het count
        df['SIZE'] = df.POS_y - df.POS_x # event size
        df['minSIZE'] = df.SIZE.rolling(window=2).min()
        
        # compute metadata, such as mafMean
        for i in df.index:
            POS1 = df.loc[i,'POS_x']
            POS2 = df.loc[i,'POS_y']
            # median is more robust mean regarding outliers
            df.loc[i,'covMean'] = dfMAF[(dfMAF.POS>=POS1)&(dfMAF.POS<=POS2)].COV.mean()
            df.loc[i,'covSEM'] = dfMAF[(dfMAF.POS>=POS1)&(dfMAF.POS<=POS2)].COV.std()/df.loc[i,'covMean']
            df.loc[i,'covMean'] = dfMAF[(dfMAF.POS>=POS1)&(dfMAF.POS<=POS2)].COV.median()
        df['meanDiff'], df['DROP'] = df.covMean.diff(), np.nan
        df['meanDiff'] = np.abs(np.power(2, df['meanDiff']) - 1)
        
        # one segment left or begin with, such as germline samples
        if len(df)==1:
            df = df[['BreakPt1','BreakPt2','POS_x', 'POS_y','SIZE','covCount', 'covMean','covSEM','meanDiff']]
            return dfSmooth, SV, df
        
        for i in df.index:
            # vafCutoff for each row
            covCutoffRow = size2cutoff_COV(dataType, df.loc[i,'minSIZE'], df.loc[i,'minCovCount'])
            
            # if dropped a row, need to skip the next row
            if i != 0 and df.loc[i-1,'DROP'] == df.loc[i-1,'DROP']:
                continue 
        
            # drop small segments FIRST
            if df.loc[i,'covCount']<minBins or df.loc[i,'SIZE']<minSize:
                smoothed = False
                df.loc[i,'DROP'] = 'Small'
                df.loc[i,'dropBP'] = drop_BP(df, i, dataCol='covMean')
            # similar mafMean with the previous segment
            elif df.loc[i,'meanDiff'] < covCutoffRow:
                smoothed = False
                df.loc[i,'DROP'] = 'covMean'
                df.loc[i,'dropBP'] = df.loc[i,'BreakPt1']
            # too few SNP in WGS, probably a poor sequencing gap such as centromere
            elif dataType=='WGS' and df.loc[i,'covCount']/df.loc[i,'SIZE'] < 0.0002:
                smoothed = False
                df.loc[i,'DROP'] = 'Gap'
                df.loc[i,'dropBP'] = drop_BP(df, i, dataCol='covMean')
            # focal events such as mycn amplifcation may have high SEM
            #elif df.loc[i,'covSEM'] > maxSEM:
            #    smoothed = False
            #    df.loc[i,'DROP'] = 'maxSEM'
            #    df.loc[i,'dropBP'] = drop_BP(df, i, dataCol='covMean')
                   
        # drop breakpoints from SV after going over df
        if smoothed == False:
            SV = [x for x in SV if x not in df.dropBP.values]
        dfSmooth = pd.concat([dfSmooth, df], sort=True, ignore_index=True)

    # to avoid return SV below
    SVs = SV
    dfSVs = df[['BreakPt1','BreakPt2','POS_x', 'POS_y','SIZE','covCount', 'covMean','covSEM','meanDiff']]
    
    return dfSmooth, SVs, dfSVs

def smooth_genome_COV(dfCOV, dictMAF, dictSV, dictDfSV, gender, dataType, minSize, minBins):
    # No normalization
    if dictMAF == {}:
        return {},{},{}
        
    print('***** COV segments smoothing *****')
    # genome wide SEM      
    dfSVall = pd.DataFrame()
    for i in dictDfSV:
        dfSVall = pd.concat([dfSVall, dictDfSV[i]], axis = 0)
    maxSEM = dfSVall.SEM.median() + 2*dfSVall.SEM.std()
    #print('maxSEM:', maxSEM)
    
    dictSmooth, dictSVs, dictDfSVs = {},{},{}
    # no coverage segmentation
    if dictMAF == {}:
        return dictSmooth, dictSVs, dictDfSVs
    
    for CHR in dictMAF:
        #print('Working on chr' + CHR + ' ...')
        if gender == 'F' and CHR == 'Y':
            dfSmooth, SVs, dfSVs = pd.DataFrame(), [], pd.DataFrame()
        # no need to smooth chrY
        elif gender == 'M' and CHR == 'Y':
            dfSmooth, SVs = pd.DataFrame(), []
            rowY = dfCOV[dfCOV.CHR == 'Y']
            if len(rowY) == 0:
                dfSVs = pd.DataFrame([[0, 3e7, 3e7, 0, -1, np.nan]], 
                                     columns=['POS_x','POS_y', 'SIZE','covCount', 'covMean', 'covSEM'])
            else:
                rowY = dfCOV[dfCOV.CHR == 'Y'].index[0]
                log2Y = np.log2(dfCOV.loc[rowY, 'covMean'])
                dfSVs = pd.DataFrame([[0, 3e7, 3e7, len(dictMAF['Y']), log2Y, np.nan]], 
                                     columns=['POS_x','POS_y', 'SIZE','covCount', 'covMean', 'covSEM'])
        else:
            dfMAF = dictMAF[CHR]
            SV = dictSV[CHR]
            dfSmooth, SVs, dfSVs = segment_smooth_COV(dfMAF, SV, dataType, maxSEM, minSize, minBins)

        dictSmooth[CHR] = dfSmooth
        dictSVs[CHR] = SVs
        dictDfSVs[CHR] = dfSVs
        
    return dictSmooth, dictSVs, dictDfSVs

################################################################
# Compute coverage events using seqDepth
################################################################
def covEVENTS(dictDfSVs_COV, seqDepth, gender, minBins, cutoff):
    print('********** Post processing COV events ... **********')
    covSegments = []
    for CHR in dictDfSVs_COV:
        dfSVs_COV = dictDfSVs_COV[CHR]
        # no segments on a chr
        if dfSVs_COV.empty:
            continue
            
        for row in dfSVs_COV.index:
            # exclude Y
            if dfSVs_COV.loc[row, 'covCount'] < minBins:
                continue
            START, END = dfSVs_COV.loc[row, 'POS_x'], dfSVs_COV.loc[row, 'POS_y']
            coverage = np.power(2, dfSVs_COV.loc[row, 'covMean'])
            covEvent, covFrac, log2ratio = coverageEvent(seqDepth, coverage, gender, CHR, cutoff)
            if covEvent != 'Diploid':
                covSegments.append([CHR, START, int(END), int(END-START+1), \
                                    round(coverage, 3), covEvent, covFrac, log2ratio])
                #print(CHR, START, END, round(coverage, 3), covEvent, covFrac)
    df = pd.DataFrame(covSegments, columns=['CHR','START','END','SIZE','COV','covEvent','covFrac','log2ratio'])
    return df

################################################################
# cov segments - vaf segments
################################################################
def COVminusVAF(covSeg, vafSeg, dictCOV, minSize, minBins, seqDepth, gender, minTF, dataType, build):
    print('********** Merging coverage and VAF segments ... **********')
    if len(covSeg) == 0:
        print('No coverage segments ...')
        return pd.DataFrame()
    
    # remove nonPARX segment from vafSeg
    [PAR1_S,PAR1_E], [PAR2_S,PAR2_E] = PAR_Xpos(build)
    nonPARX = ((vafSeg.CHR == 'X') & (vafSeg.START == PAR1_E) & (vafSeg.END == PAR2_S))
    if nonPARX.sum() == 1:
        vafSeg = vafSeg[~nonPARX].copy()

    covSegBed = pybedtools.BedTool.from_dataframe(covSeg)
    vafSeg.END = vafSeg.END.astype(int)
    vafSegBed = pybedtools.BedTool.from_dataframe(vafSeg)
    # https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html
    # cov unique events, no overlaps at all
    df1 = covSegBed.subtract(vafSegBed, A=True)
    # cov - vaf
    df2 = covSegBed.subtract(vafSegBed)
    # (cov - vaf) - cov unique = hanging regions
    df = df2.subtract(df1, A=True)
    
    df1 = df1.to_dataframe()
    if len(df1) == 0:
        df1 = pd.DataFrame([],
               columns=['CHR', 'START', 'END', 'SIZE', 'COV', 'covEvent', 'covFrac', 'log2ratio'])
    else:
        df1.columns = ['CHR', 'START', 'END', 'SIZE', 'COV', 'covEvent', 'covFrac', 'log2ratio']
        df1.CHR = df1.CHR.astype(str)
    
    # filter small events
    df = df.to_dataframe()
    if df.empty:
        print('No overlap hanging regions ...')
        return df1
        
    df.columns = ['CHR', 'START', 'END', 'SIZE', 'COV', 'covEvent', 'covFrac', 'log2ratio']
    df.CHR = df.CHR.astype(str)
    df['SIZE'] = df['END'] - df['START']

    for i in df.index:
        CHR = str(df.loc[i, 'CHR'])
        dfCOV = dictCOV[CHR]
        START = df.loc[i, 'START']
        END = df.loc[i, 'END']
        df.loc[i, 'covCount'] = ((dfCOV.POS >= START) & (dfCOV.POS <= END)).sum()

    df = df[(df.covCount >= minBins) & (df['SIZE'] >= minSize)]
    if len(df) == 0:
        print('Found overlap, but small hanging regions ...')
        return df1

    # compute coverage and TF
    for i in df.index:
        CHR = str(df.loc[i, 'CHR'])
        dfCOV = dictCOV[CHR]
        START = df.loc[i, 'START']
        END = df.loc[i, 'END']
        covCountCOV, covMean, covSEM = computeCOV(dfCOV, START, END, dataType)
        covEvent, covFrac, log2ratio = coverageEvent(seqDepth, covMean, gender, CHR, minTF)
        df.loc[i, 'COV'] = covMean
        df.loc[i, 'covEvent'] = covEvent
        df.loc[i, 'covFrac'] = covFrac
        df.loc[i, 'log2ratio'] = log2ratio
    df = df[df.covFrac >= minTF]
    if len(df)>0:
        print('Found overlap hanging regions ...')
    else:
        print('Found overlap, but lowTF hanging regions ...')
    
    return pd.concat([df, df1], axis=0, sort=True, ignore_index=True)

################################################################
# Estimates of tumor fraction *** May need revision
################################################################
def tumorFrac(vafSeg, covSeg, seqDepth, minTF=0.1):
    print('********** Estimating tumor fraction ... **********')
    df = vafSeg.copy()
    # exclude amplification and samll events
    if ((df.vafFrac < 1.05) & (df.SIZE > 3e7)).sum() > 3:
        TF = df[(df.vafFrac < 1.05) & (df.SIZE > 2e7)].vafFrac.median()
        print('TF estimation senario 1: >3 TwoNorm')
    elif 0 < ((df.vafFrac < 1.05) & (df.SIZE > 3e7)).sum() <= 3:
        TF = df[(df.vafFrac < 1.05) & (df.SIZE > 2e7)].vafFrac.max()
        print('TF estimation senario 2: <=3 TwoNorm')
    # exclude amplification
    elif (df.vafFrac <= 1).sum()  > 0:
        TF = df[df.vafFrac < 1].vafFrac.max()
        print('TF estimation senario 3: small TwoNorm')
    elif ((df.CHR!= 'X') & (df.vafEvent != 'ROH')).sum() > 0: # exclude male chrX
        TF = df[(df.CHR!= 'X') & (df.vafEvent != 'ROH')].vafFrac.max()/2
        print('TF estimation senario 4: vafFrac > 1')
    elif len(covSeg) > 0:
        TF = np.abs(covSeg.covFrac).median()
        print('TF estimation senario 5: covEvents only')
    else:
        TF = 'Unknown'
        final = pd.DataFrame([],
               columns=['CHR', 'START', 'END', 'SIZE', 'covMean','covLog2Ratio', 'EVENT', 'Allelic_Imbalance', 
                        'chrA_freq', 'chrB_freq', 'chrA_CNV', 'chrB_CNV', 'chrA_TF', 'chrB_TF'])
        return 0, final
    print('TF estimation:', round(TF, 3))
    
    #print('Working on VAF segments ...')
    if len(df) == 0:
        df = pd.DataFrame([],
               columns=['CHR', 'START', 'END', 'SIZE', 'covMean','covLog2Ratio', 'EVENT', 'Allelic_Imbalance', 
                        'chrA_freq', 'chrB_freq', 'chrA_CNV', 'chrB_CNV', 'chrA_TF', 'chrB_TF'])
    else:
        df['chrA_CNVraw'] = df.chrA_freq/TF
        df['chrB_CNVraw'] = df.chrB_freq/TF
        df['chrA_CNV'],  df['chrB_CNV'] = np.nan, np.nan
        df['chrA_TF'],   df['chrB_TF'] = np.nan, np.nan
        df['Allelic_Imbalance'] = 'yes'
        for i in df.index:
            # Tumor fraction and CNV
            if df.loc[i, 'cnvEvent'] == 'cnLOH':
                df.loc[i, 'chrA_CNV'], df.loc[i, 'chrB_CNV'] = -1, 1
            elif df.loc[i, 'cnvEvent'] == 'Gain':
                if df.loc[i, 'chrB_CNVraw'] <= 1:
                    df.loc[i, 'chrA_CNV'], df.loc[i, 'chrB_CNV'] = np.nan, 1
                else:
                    df.loc[i, 'chrA_CNV'], df.loc[i, 'chrB_CNV'] = np.nan, round(df.loc[i, 'chrB_CNVraw'])
            elif df.loc[i, 'cnvEvent'] == 'Loss':
                if df.loc[i, 'chrB_CNVraw'] >= -1:
                    df.loc[i, 'chrA_CNV'], df.loc[i, 'chrB_CNV'] = np.nan, -1
                else:
                    df.loc[i, 'chrA_CNV'], df.loc[i, 'chrB_CNV'] = np.nan, round(df.loc[i, 'chrB_CNVraw'])

            elif df.loc[i, 'cnvEvent'] in ['complexGain', 'complexLoss']:
                # chrA
                if 0 < df.loc[i, 'chrA_CNVraw'] <= 1:
                    df.loc[i, 'chrA_CNV'] = 1
                elif -1 < df.loc[i, 'chrA_CNVraw'] < 0:
                    df.loc[i, 'chrA_CNV'] = -1
                else:
                    df.loc[i, 'chrA_CNV'] = round(df.loc[i, 'chrA_CNVraw'])

                # chrB
                if 0 < df.loc[i, 'chrB_CNVraw'] <= 1:
                    df.loc[i, 'chrB_CNV'] = 1
                elif -1 < df.loc[i, 'chrB_CNVraw'] < 0:
                    df.loc[i, 'chrB_CNV'] = -1
                else:
                    df.loc[i, 'chrB_CNV'] = round(df.loc[i, 'chrB_CNVraw'])

            # max loss = -2
            df.loc[i, 'chrA_CNV'] = np.max([df.loc[i, 'chrA_CNV'], -2])
            df.loc[i, 'chrB_CNV'] = np.max([df.loc[i, 'chrB_CNV'], -2])

            # tumor fraction
            if df.loc[i, 'cnvEvent'] == 'cnLOH':
                df.loc[i, 'chrB_TF'] = df.loc[i, 'chrB_freq']
            else:
                df.loc[i, 'chrA_TF'] = df.loc[i, 'chrA_freq']/df.loc[i, 'chrA_CNV']
                df.loc[i, 'chrB_TF'] = df.loc[i, 'chrB_freq']/df.loc[i, 'chrB_CNV']

        #print(df.columns)
        df = df[['CHR', 'START', 'END', 'SIZE', 'covMean', 'log2ratio', 'cnvEvent', 'Allelic_Imbalance', 
                 'chrA_freq', 'chrB_freq', 'chrA_CNV', 'chrB_CNV', 'chrA_TF', 'chrB_TF']]
        df.columns = ['CHR', 'START', 'END', 'SIZE', 'covMean','covLog2Ratio', 'EVENT', 'Allelic_Imbalance', 
                 'chrA_freq', 'chrB_freq', 'chrA_CNV', 'chrB_CNV', 'chrA_TF', 'chrB_TF']
    
    # coverage
    #print('Working on coverage segments ...')
    if len(covSeg) == 0:
        return round(TF, 3), df
    
    # filter OneNorm loss if less than 0.4
    if minTF >= 0.1:
        lowLoss = ((covSeg.covFrac > -0.4) & (covSeg.covEvent == 'Loss'))
        lowGain = ((covSeg.covFrac < 0.2) & (covSeg.covEvent == 'Gain'))
        df1 = covSeg[(~lowLoss) & (~lowGain)].copy()
    else:
        lowLoss = ((covSeg.covFrac > -0.4) & (covSeg.covEvent == 'Loss'))
        df1 = covSeg[~lowLoss].copy()
        
    if len(df1) == 0:
        return round(TF, 3), df
    df1['CNVraw'] = abs(df1.covFrac/TF)
    df1['Allelic_Imbalance'] = 'no'
    df1['chrA_freq'], df1['chrB_freq'] = np.nan, np.nan
    df1['chrA_CNV'], df1['chrB_CNV'] = np.nan, np.nan
    df1['chrA_TF'], df1['chrB_TF'] = 0, 0
    for i in df1.index:
        if df1.loc[i, 'CNVraw'] <= 1:
            df1.loc[i, 'CNV'] = 1
        else:
            df1.loc[i, 'CNV'] = round(df1.loc[i, 'CNVraw'])
        df1.loc[i, 'Tumor_Fraction'] = abs(df1.loc[i, 'covFrac']/df1.loc[i, 'CNV'])
    #print(df1.columns)
    if 'COV' in df1.columns:
        df1 = df1[['CHR', 'START', 'END', 'SIZE', 'COV', 'log2ratio', 'covEvent', 'Allelic_Imbalance', 
               'covFrac', 'CNV', 'Tumor_Fraction']]
    else:
        df1 = df1[['CHR', 'START', 'END', 'SIZE', 'covMean', 'log2ratio', 'covEvent', 'Allelic_Imbalance', 
               'covFrac', 'CNV', 'Tumor_Fraction']]
    df1.columns = ['CHR', 'START', 'END', 'SIZE', 'covMean', 'covLog2Ratio', 'EVENT', 'Allelic_Imbalance',
                'chrB_freq', 'chrB_CNV', 'chrB_TF']
    
    #print('Merging segments ...')
    final = pd.concat([df, df1], axis=0, sort=False, ignore_index=True)
    
    # sort
    final['CHRsort'] = np.nan
    for i in final.index:
        CHR = final.loc[i, 'CHR']
        if CHR == 'X':
            final.loc[i, 'CHRsort'] = 23
        elif CHR == 'Y':
            final.loc[i, 'CHRsort'] = 24
        else:
            final.loc[i, 'CHRsort'] = int(CHR)
            
    # formatting
    final.sort_values(by=['CHRsort', 'START'], inplace=True)
    final = final[['CHR', 'START', 'END', 'SIZE', 'covMean', 'covLog2Ratio', 'EVENT', 'Allelic_Imbalance', 
             'chrA_freq', 'chrA_CNV', 'chrA_TF', 'chrB_freq', 'chrB_CNV', 'chrB_TF']]

    for i in final.index:
        final.loc[i, 'covMean'] = ROUND(final.loc[i, 'covMean'])
        final.loc[i, 'covLog2Ratio'] = ROUND(final.loc[i, 'covLog2Ratio'])
        final.loc[i, 'chrA_freq'] = ROUND(final.loc[i, 'chrA_freq'])
        final.loc[i, 'chrB_freq'] = ROUND(final.loc[i, 'chrB_freq'])
        final.loc[i, 'chrA_TF'] = ROUND(final.loc[i, 'chrA_TF'])
        final.loc[i, 'chrB_TF'] = ROUND(final.loc[i, 'chrB_TF'])

        if final.loc[i, 'chrA_CNV'] > 0:
            final.loc[i, 'chrA_CNV'] = '+' + str(int(final.loc[i, 'chrA_CNV']))
        elif final.loc[i, 'chrA_CNV'] < 0:
            final.loc[i, 'chrA_CNV'] = str(int(final.loc[i, 'chrA_CNV']))
            
        if final.loc[i, 'chrB_CNV'] > 0:
            final.loc[i, 'chrB_CNV'] = '+' + str(int(final.loc[i, 'chrB_CNV']))
        elif final.loc[i, 'chrB_CNV'] < 0:
            final.loc[i, 'chrB_CNV'] = str(int(final.loc[i, 'chrB_CNV']))

    return round(TF, 3), final

def ROUND(value, decimal=3):
    if value == value:
        return round(value, decimal)
    else:
        return value

################################################################
# Segments visualization
################################################################
def plotSegments(dictSNV, dictDfSV, dictDfSVs, dictCOV, dictDfSV_COV, dictDfSVs_COV, CHR, \
                 XLIM=(0,1e9), markSize=3, ALPHA=1, \
                 plotDfSV=True, plotDfSVs=True, data='VAF'):
    dfSNV = dictSNV[CHR]
    dfSV  = dictDfSV[CHR]
    dfSVs = dictDfSVs[CHR]
    dfCOV = dictCOV[CHR]
    dfSV_COV = dictDfSV_COV[CHR]
    dfSVs_COV = dictDfSVs_COV[CHR]

    if XLIM==(0,1e9):
        if CHR == 'Y':
            XLIM = (0, 3e7)
        else:
            total_len = dfSNV.POS.max()-dfSNV.POS.min()
            minPOS = np.min([dfSNV.POS.min(), dfCOV.POS.min()])
            maxPOS = np.max([dfSNV.POS.max(), dfCOV.POS.max()])
            XLIM = (minPOS - 1e6 ,maxPOS + 1e6)
    
    # plot VAF & COV
    fig,ax = plt.subplots(figsize=(14,1.5))       
    ax.plot(dfCOV.POS, dfCOV.COV, '.', ms=markSize, alpha=ALPHA, c='darkorange')
    ax2 = ax.twinx()
    ax2.set(ylim=(-0.03,1.03), xlim=XLIM)
    ax2.plot(dfSNV.POS, dfSNV.VAF, '.', ms=markSize, alpha=ALPHA)
    ax.set_ylabel(CHR,fontweight ='bold',fontsize = 12,rotation=0,va="center",ha="right")
    
    # plot dfSV
    if data == 'VAF':
        dataSV = dfSV
        dataSVs = dfSVs
    else:
        dataSV = dfSV_COV
        dataSVs = dfSVs_COV
    
    if plotDfSV == True:
        for i in dataSV.index:
            ax2.axvline(x=dataSV.loc[i,'POS_x'], color='darkgreen', alpha=0.5)
            ax2.axvline(x=dataSV.loc[i,'POS_y'], color='darkgreen', alpha=0.5)
    if plotDfSVs == True:
        for i in dataSVs.index:
            ax2.axvline(x=dataSVs.loc[i,'POS_x'], color='red', alpha=1)
            ax2.axvline(x=dataSVs.loc[i,'POS_y'], color='red', alpha=1)


################################################################
# SubChrom visualization
################################################################
def plotCoverage(ax, df, genoInfo, seqDepth, smooth=True, xTickLabel=False):
    if smooth == True:
        covColumn = 'covSmooth'
    else:
        covColumn = 'covMean'
    for i in df.index:
        x1 = df.loc[i,'plotSTART']
        x2 = df.loc[i,'plotEND']
        y  = df.loc[i,covColumn]
        CHR = df.loc[i,'CHR']
        ax.fill_between([x1,x2], [y,y], color=genoInfo.loc[CHR, 'Color']) 
    
    # y-axis limit, ticks, labels
    covMax = np.max(df[covColumn])
    ax.set(ylim=(0,covMax*1.1))
    
    step = covMax/3.5 # divide to 3-4 segments
    LOG10 = np.floor(np.log10(step)) # number of digit minus 1
    INDEX = np.power(10, LOG10)
    STEP = np.round(step/INDEX)*INDEX
    ax.set_yticks(np.arange(0, covMax*1.1, step=STEP))
    ax.set_ylabel("Coverage", fontweight='bold')
    
    
    # x-axis limit, ticks, labels
    plotGenemoeSize = genoInfo.loc['Y', 'prevLen'] + genoInfo.loc['Y', 'chrSize']
    ax.set(xlim=(-3e6, plotGenemoeSize + 3e6))
    
    labels = genoInfo.index
    values = genoInfo.loc[labels].plotPos.values
    ax.set_xticks(values)
    if xTickLabel == True:        
        ax.set_xticklabels(labels)
        ax.set_xlabel('Chromosome', fontweight='bold')
    elif xTickLabel == False:
        ax.set_xticklabels([])

    # seqDepth (diploid depth)
    ax.axhline(y=seqDepth, color='black', linestyle = 'dotted',linewidth=1.5)
     
        
def plotVAF(ax, dictSNV, genoInfo, gender, seqDepth, minCOV=-1, xTickLabel=False):
    # build
    if genoInfo.loc['1', 'chrSize'] == 249250621:
        build = 'hg19'
    elif genoInfo.loc['1', 'chrSize'] == 248956422:
        build = 'hg38'
        
    # minCOV for plotting
    if minCOV == -1:
        minCOV = int(seqDepth/4)
    else:
        minCOV = int(minCOV)
    print('minCOV for plotting VAF:', minCOV)
    
    # number of het
    hetCount = genoInfo.hetCount_1Mb.max()
    if hetCount <= 10: # cfDNA 2-10
        markerSize=1.5
        plotAlpha=1
    elif hetCount <= 30: # WES 14
        markerSize=1
        plotAlpha=0.9
    elif hetCount <= 100:
        markerSize=1
        plotAlpha=0.7
    elif hetCount <= 300:
        markerSize=1
        plotAlpha=0.1
    else: # WGS
        markerSize=1
        plotAlpha=0.01
    #print(hetCount, markerSize, plotAlpha)
          
    for CHR in dictSNV:
        df = dictSNV[CHR]
        if gender == 'M' and CHR == 'X':
            [PAR1_S,PAR1_E] = PAR_Xpos(build)[0]
            [PAR2_S,PAR2_E] = PAR_Xpos(build)[1]
            df1 = df[(df.POS>PAR1_S) & (df.POS<PAR1_E)]
            # PAR1, more alpha due to few data
            ax.plot(df1.PlotPos, df1.VAF, 
                 'o', ms=markerSize,markeredgewidth=0,alpha=plotAlpha*3,
                 color=genoInfo.loc[CHR, 'Color'])
            # X
            df2 = df[(df.POS>PAR1_E) & (df.POS<PAR2_S)]
            ax.plot(df2.PlotPos, df2.VAF, 
                 'o', ms=markerSize,markeredgewidth=0,alpha=plotAlpha,
                 color=genoInfo.loc[CHR, 'Color'])
            # PAR2, more alpha due to few data
            df3 = df[(df.POS>PAR2_S) & (df.POS<PAR2_E)]
            ax.plot(df3.PlotPos, df3.VAF, 
                 'o', ms=markerSize,markeredgewidth=0,alpha=1,
                 color=genoInfo.loc[CHR, 'Color'])
        else:
            df = dictSNV[CHR].copy()
            df = df[df.COV >= minCOV]
            ax.plot(df.PlotPos, df.VAF, 
                     'o', ms=markerSize,markeredgewidth=0,alpha=plotAlpha,
                     color=genoInfo.loc[CHR, 'Color']) 
    
    # y-axis limit, ticks, labels
    ax.set(ylim=(-0.02, 1.02), yticks=([0, 0.5, 1]))
    if minCOV != int(seqDepth/4):
        ax.set_ylabel('VAF\n(minCOV=' + str(minCOV) + ')', fontweight='bold')
    else:
        ax.set_ylabel('VAF', fontweight='bold')
    
    # x-axis limit, ticks, labels
    plotGenemoeSize = genoInfo.loc['Y', 'prevLen'] + genoInfo.loc['Y', 'chrSize']
    ax.set(xlim=(-3e6, plotGenemoeSize + 3e6))
    
    labels = genoInfo.index
    values = genoInfo.loc[labels].plotPos.values
    ax.set_xticks(values)
    if xTickLabel == True:        
        ax.set_xticklabels(labels)
        ax.set_xlabel('Chromosome', fontweight='bold')
    elif xTickLabel == False:
        ax.set_xticklabels([])

def plotLog2Ratio(ax, dictCOV1, df, genoInfo, gender, seqDepth, xTickLabel=False):
    # build
    if genoInfo.loc['1', 'chrSize'] == 249250621:
        build = 'hg19'
    elif genoInfo.loc['1', 'chrSize'] == 248956422:
        build = 'hg38'
        
    # number of het
    hetCount = genoInfo.hetCount_1Mb.max()
    if hetCount <= 10: # cfDNA 2-10
        markerSize=2
        plotAlpha=0.5
    elif hetCount <= 30: # WES 14
        markerSize=1
        plotAlpha=0.5
    elif hetCount <= 100:
        markerSize=1
        plotAlpha=0.3
    elif hetCount <= 300:
        markerSize=1
        plotAlpha=0.1
    else: # WGS
        markerSize=1
        plotAlpha=0.01
    #print(hetCount, markerSize, plotAlpha)
        
    dictCOV = copy.deepcopy(dictCOV1)
    # auxilary lines
    ax.axhline(y=0, c='black', linewidth=1)
    for CHR in genoInfo.index:
        if CHR not in ['1']:
            ax.axvline(x=genoInfo.loc[CHR, 'prevLen']-3e6, c='black', linewidth=1, ls='dotted')
    
    # coverage points
    for CHR in dictCOV:
        dfCOV = dictCOV[CHR]
        dfCOV['log2COV'] = np.nan
        if gender == 'M' and CHR == 'X':
            [PAR1_S,PAR1_E] = PAR_Xpos(build)[0]
            [PAR2_S,PAR2_E] = PAR_Xpos(build)[1]
            PAR = dfCOV[(dfCOV.POS<=PAR1_E) | (dfCOV.POS>=PAR2_S)].index
            nonPAR = dfCOV[(dfCOV.POS>PAR1_E) & (dfCOV.POS<PAR2_S)].index
            dfCOV.loc[PAR, 'log2COV'] =  np.log2(dfCOV.loc[PAR, 'COV']/seqDepth)
            dfCOV.loc[nonPAR, 'log2COV'] =  np.log2(2 * dfCOV.loc[nonPAR, 'COV']/seqDepth)
        elif gender == 'M' and CHR == 'Y':
            dfCOV['log2COV'] = np.log2(2 * dfCOV.COV/seqDepth)
            markerSize=1
            plotAlpha=1
        elif gender == 'F' and CHR == 'Y':
            continue
        else:
            dfCOV['log2COV'] = np.log2(dfCOV.COV/seqDepth)
        ax.plot(dfCOV.PlotPos, dfCOV.log2COV, 'o', ms=markerSize, 
                markeredgewidth=0, alpha=plotAlpha, color='darkgrey')
  
    # lines
    for i in df.index:
        AI = df.loc[i, 'Allelic_Imbalance']
        if AI == 'no':
            if df.loc[i,'EVENT'] == 'Loss':
                COLOR = 'darkgreen'
            elif df.loc[i,'EVENT'] == 'Gain':
                COLOR = 'maroon'
            else:
                COLOR = 'darkorange'
        else:
            if df.loc[i,'EVENT'] in ['Loss', 'complexLoss']:
                COLOR = 'limegreen'
            elif df.loc[i,'EVENT'] in ['Gain', 'complexGain']:
                COLOR = 'red'
            elif df.loc[i,'EVENT'] == 'cnLOH':
                COLOR = 'blue'
            else:
                COLOR = 'orange'

        CHR = str(df.loc[i,'CHR'])
        x1 = df.loc[i,'START'] + genoInfo.loc[CHR, 'prevLen']
        x2 = df.loc[i,'END'] + genoInfo.loc[CHR, 'prevLen']
        y  = df.loc[i,'covLog2Ratio']
       
        # lines
        plotGap = 3e6
        if x2 - x1 < plotGap:
            ax.plot(x1+(x2-x1)/2, y, 'o', ms=2.5, color=COLOR)
        elif x2 - x1 < 2*plotGap:
            ax.plot([x1+plotGap/2,x2-plotGap/2], [y,y], c=COLOR, linewidth=2.5)
        else:
            ax.plot([x1+plotGap,x2-plotGap], [y,y], c=COLOR, linewidth=2.5)
    
    # y-axis limit, ticks, labels
    if len(df) == 0:
        ax.set(ylim=(-2, 2))
    else:
        lossMax = df.covLog2Ratio.min()
        lossMax = np.floor(lossMax)
        covMin = np.min([-2, lossMax-1])
        gainMax = df.covLog2Ratio.max()
        gainMax = np.round(gainMax)
        covMax = np.max([2, gainMax+1])
        #print(lossMax, gainMax)
        ax.set(ylim=(covMin, covMax))
        ax.set_yticks(np.arange(covMin, covMax+0.1, step=1)) # +0.1 to include covMax
    ax.set_ylabel("Copy ratio (log2)", fontweight='bold')
        
    # x-axis limit, ticks, labels
    plotGenemoeSize = genoInfo.loc['Y', 'prevLen'] + genoInfo.loc['Y', 'chrSize']
    ax.set(xlim=(-3e6, plotGenemoeSize + 3e6))
    
    labels = genoInfo.index
    values = genoInfo.loc[labels].plotPos.values
    ax.set_xticks(values)
    if xTickLabel == True:        
        ax.set_xticklabels(labels)
        ax.set_xlabel('Chromosome', fontweight='bold')
    elif xTickLabel == False:
        ax.set_xticklabels([])
        
def plotFreq(ax, df, genoInfo, xTickLabel=False):
    # auxilary lines
    for CHR in genoInfo.index:
        if CHR not in ['1']:
            ax.axvline(x=genoInfo.loc[CHR, 'prevLen']-3e6, c='black', linewidth=1, ls='dotted')
    for i in [0, 0.2, 0.4, 0.6,  0.8]:
        ax.axhline(y=i, c='grey', linewidth=1, alpha=0.2)
    
    # lines
    for i in df.index:
        AI = df.loc[i, 'Allelic_Imbalance']
        chrA_COLOR, chrB_COLOR = 'darkorange', 'darkorange'
        if AI == 'no':
            if df.loc[i,'EVENT'] == 'Loss':
                chrB_COLOR = 'darkgreen'
            elif df.loc[i,'EVENT'] == 'Gain':
                chrB_COLOR = 'maroon'
        else:
            if df.loc[i,'EVENT'] in ['Loss']:
                chrB_COLOR = 'limegreen'
            elif df.loc[i,'EVENT'] in ['Gain']:
                chrB_COLOR = 'red'
            elif df.loc[i,'EVENT'] == 'cnLOH':
                chrB_COLOR = 'blue'
            elif df.loc[i,'EVENT'] in ['complexLoss', 'complexGain']:
                if df.loc[i,'chrA_freq'] > 0:
                    chrA_COLOR = 'red'
                else:
                    chrA_COLOR = 'limegreen'
                if df.loc[i,'chrB_freq'] > 0:
                    chrB_COLOR = 'red'
                else:
                    chrB_COLOR = 'limegreen'

        CHR = str(df.loc[i,'CHR'])
        x1 = df.loc[i,'START'] + genoInfo.loc[CHR, 'prevLen']
        x2 = df.loc[i,'END'] + genoInfo.loc[CHR, 'prevLen']
        chrA_y = abs(df.loc[i,'chrA_TF'])
        chrB_y = abs(df.loc[i,'chrB_TF'])
       
        # plot lines
        plotGap = 3e6
        if x2 - x1 < plotGap:
            ax.plot(x1+(x2-x1)/2, chrA_y, 'o', ms=2.5, color=chrA_COLOR)
            ax.plot(x1+(x2-x1)/2, chrB_y, 'o', ms=2.5, color=chrB_COLOR)
        elif x2 - x1 < 2*plotGap:
            ax.plot([x1+plotGap/2,x2-plotGap/2], [chrA_y,chrA_y], color=chrA_COLOR, linewidth=2.5)
            ax.plot([x1+plotGap/2,x2-plotGap/2], [chrB_y,chrB_y], color=chrB_COLOR, linewidth=2.5)
        else:
            ax.plot([x1+plotGap,x2-plotGap], [chrA_y,chrA_y], color=chrA_COLOR, linewidth=2.5)
            ax.plot([x1+plotGap,x2-plotGap], [chrB_y,chrB_y], color=chrB_COLOR, linewidth=2.5)
    
    # y-axis limit, ticks, labels
    ax.set(ylim=(0, 1.03), yticks=([0, 0.2, 0.4, 0.6,  0.8, 1]))
    ax.set_ylabel("Tumor fraction", fontweight='bold')
    
    # x-axis limit, ticks, labels
    plotGenemoeSize = genoInfo.loc['Y', 'prevLen'] + genoInfo.loc['Y', 'chrSize']
    ax.set(xlim=(-3e6, plotGenemoeSize + 3e6))
    
    labels = genoInfo.index
    values = genoInfo.loc[labels].plotPos.values
    ax.set_xticks(values)
    if xTickLabel == True:        
        ax.set_xticklabels(labels)
        ax.set_xlabel('Chromosome', fontweight='bold')
    elif xTickLabel == False:
        ax.set_xticklabels([])
        
def plotFinal(dictSNV,dictCOV,dfCOV,df,genoInfo,gender,seqDepth,SAMPLE,dataType,TF,minCOV=-1,plotTF=True,SAVE=False):
    print('***** Working on SubChrom visualization ... *****')
    if plotTF:
        fig,ax = plt.subplots(4,1,figsize=(15,4*1.8))
        plt.subplots_adjust(hspace=0.12)
        ax[0].set_title(SAMPLE + '.' + dataType + '.SubChrom', fontweight='bold')
        plotCoverage(ax[0], dfCOV, genoInfo, seqDepth)
        plotVAF(ax[1], dictSNV, genoInfo, gender, seqDepth, minCOV)
        plotLog2Ratio(ax[2], dictCOV, df, genoInfo, gender, seqDepth)
        plotFreq(ax[3], df, genoInfo, xTickLabel=True)
        yPos = 0.04
        start = 0.135
        fig.add_artist(Text(0.78, yPos, text='Tumor fraction: '+ str(TF), size=11, va='center'))
    else:
        fig,ax = plt.subplots(3,1,figsize=(15,3*1.8))
        plt.subplots_adjust(hspace=0.12)
        ax[0].set_title(SAMPLE + '.' + dataType + '.SubChrom', fontweight='bold')
        plotCoverage(ax[0], dfCOV, genoInfo, seqDepth)
        plotVAF(ax[1], dictSNV, genoInfo, gender, seqDepth, minCOV)
        plotLog2Ratio(ax[2], dictCOV, df, genoInfo, gender, seqDepth, xTickLabel=True)
        #yPos = 0.015
        #start = 0.2
        yPos = 0.04
        start = 0.135
        fig.add_artist(Text(0.78, yPos, text='Tumor fraction: '+ str(TF), size=11, va='center'))
        
    # legend
    size  = 0.02
    size2 = 0.025
    space = 0.11
    fig.add_artist(lines.Line2D([start, start + size], [yPos, yPos], c='black', lw=1.5, ls='dotted'))
    fig.add_artist(lines.Line2D([start+space*1, start+space*1 + size], [yPos, yPos], c='red', lw=2.5))
    fig.add_artist(lines.Line2D([start+space*2, start+space*2 + size], [yPos, yPos], c='maroon', lw=2.5))
    fig.add_artist(lines.Line2D([start+space*3, start+space*3 + size], [yPos, yPos], c='limegreen', lw=2.5))
    fig.add_artist(lines.Line2D([start+space*4, start+space*4 + size], [yPos, yPos], c='darkgreen', lw=2.5))
    fig.add_artist(lines.Line2D([start+space*5, start+space*5 + size], [yPos, yPos], c='blue', lw=2.5))
    
    fig.add_artist(Text(start + size2, yPos, text='Diploid depth', size=11, va='center'))
    fig.add_artist(Text(start+space*1 + size2, yPos, text='Gain (with AI)', size=11, va='center', c='red'))
    fig.add_artist(Text(start+space*2 + size2, yPos, text='Gain (no AI)', size=11, va='center', c='maroon'))
    fig.add_artist(Text(start+space*3 + size2, yPos, text='Loss (with AI)', size=11, va='center', c='limegreen'))
    fig.add_artist(Text(start+space*4 + size2, yPos, text='Loss (no AI)', size=11, va='center', c='green'))
    fig.add_artist(Text(start+space*5 + size2, yPos, text='cnLOH', size=11, va='center', c='blue'))
    
    if SAVE == False:
        plt.show()
    else:
        if os.path.exists('results') == False:
            os.mkdir('results')
        plt.subplots_adjust(bottom=0.13, top=0.95, left=0.06, right=0.98)
        plt.savefig('results/' + SAMPLE + '.' + dataType + '.CNV.png', dpi=300)
        #plt.savefig('results/' + SAMPLE + '.' + dataType + '.CNV.pdf', format="pdf", dpi=300)
        plt.close()

################################################################
# plot allelic imbalance (AI) figure
################################################################
def plot_chr(df1,df2,yvalue,ax,color,markerSize, plotAlpha):
    ax.plot(df1.POS,df1[yvalue],'.',c=color,ms=markerSize, alpha=plotAlpha)
    ax.plot(df2.POS_new,df2[yvalue],'.',c=color,ms=markerSize, alpha=plotAlpha)
    ax.tick_params(axis='y', colors=color)
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.spines['bottom'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
def plot_chr_pair(dictSNV,dictCOV,final,chr_pair,ax1,seqDepth,genoInfo,dataType,dfGenes):
    # number of het
    hetCount = genoInfo.hetCount_1Mb.median()
    if hetCount <= 12: # cfDNA 2-10
        markerSize=1.5
        plotAlphaVAF=1
        plotAlphaCOV=1
    elif hetCount <= 30: # WES 14
        markerSize=1
        plotAlphaVAF=0.5
        plotAlphaCOV=0.1
    elif hetCount <= 100:
        markerSize=1
        plotAlphaVAF=0.3
        plotAlphaCOV=0.1
    elif hetCount <= 300:
        markerSize=1
        plotAlphaVAF=0.1
        plotAlphaCOV=0.1
    else: # WGS
        markerSize=1
        plotAlphaVAF=0.05
        plotAlphaCOV=0.05
    #print(hetCount, markerSize, plotAlpha)    
        
    df1=dictSNV[str(chr_pair[0])].copy()
    df1COV=dictCOV[str(chr_pair[0])].copy()
    df2=dictSNV[str(chr_pair[1])].copy()
    df2COV=dictCOV[str(chr_pair[1])].copy()
    # second chromosome size
    secondChrSize = genoInfo.loc[str(chr_pair[1]), 'chrSize']
    # new position for the second chromosome
    df2['POS_new']=df2.POS+(310000000-secondChrSize)
    df2COV['POS_new']=df2COV.POS+(310000000-secondChrSize)
    
    yTicks = np.arange(0, seqDepth*3, step=(round(seqDepth,-1) + 10)) # round up to the nearest ten and add 10
    ax1.set(ylim=(0,seqDepth*3),yticks=yTicks,xlim=(0,310000000))
    ax1.set_ylabel(str(chr_pair[0]),fontweight ='bold',fontsize = 12,rotation=0,va="center",ha="right")
    plot_chr(df1COV,df2COV,'COV',ax1,'darkorange',markerSize, plotAlphaCOV)
    
    ax2 = ax1.twinx()
    ax2.set(yticks=[0,0.5,1], ylim=(-0.03,1.03))
    ax2.set_ylabel(str(chr_pair[1]),fontweight ='bold',fontsize = 12,rotation=0,va="center",ha="left")
    plot_chr(df1,df2,'VAF',ax2,'C0',markerSize, plotAlphaVAF)
    
    # plot genes
    ax3 = ax1.twinx()
    dfGenes1 = dfGenes[dfGenes.CHR == str(chr_pair[0])].copy()
    dfGenes2 = dfGenes[dfGenes.CHR == str(chr_pair[1])].copy()
    dfGenes2['POS']=dfGenes2['POS']+(310000000-secondChrSize)
    plot_genes(dfGenes1, ax2)
    plot_genes(dfGenes2, ax2)
    
    # plot lines
    dfL1 = final[final.CHR == str(chr_pair[0])].copy()
    dfL2 = final[final.CHR == str(chr_pair[1])].copy()
    dfL2['START']=dfL2['START']+(310000000-secondChrSize)
    dfL2['END']=dfL2['END']+(310000000-secondChrSize)
    plot_AI_seg(dfL1, ax3)
    plot_AI_seg(dfL2, ax3)
    ax3.set(ylim=(0,seqDepth*3),yticks=yTicks,xlim=(0,310000000))
    ax3.set_yticks([])
    ax3.set_yticklabels([])
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
def plot_AI_seg(df, ax):
    for i in df.index:
        x1 = df.loc[i,'START']
        x2 = df.loc[i,'END']
        y = df.loc[i,'covMean']
        if df.loc[i, 'EVENT'] in ['Loss', 'complexLoss']:
            if df.loc[i, 'Allelic_Imbalance'] == 'yes':
                COLOR = 'limegreen'
            else:
                COLOR = 'green'     
        elif df.loc[i, 'EVENT'] in ['Gain', 'complexGain']:
            if df.loc[i, 'Allelic_Imbalance'] == 'yes':
                COLOR = 'red'
            else:
                COLOR = 'maroon'
        elif df.loc[i, 'EVENT'] == 'cnLOH':
            COLOR = 'blue'
        else:
            COLOR = 'black'
        if x2 - x1 < 1e6:
            ax.plot(x1+(x2-x1)/2, y, 'o', ms=3, color=COLOR)
        else:
            ax.plot([x1,x2], [y,y], c=COLOR, linewidth=2)

def plot_genes(dfGenes, ax):
    if not dfGenes.empty:
        for i in dfGenes.index:       
            ax.text(dfGenes.loc[i, 'POS'], -0.1, dfGenes.loc[i, 'GENE'], ha='center', va='top', fontsize=9)
            ax.plot(dfGenes.loc[i, 'POS'], 0.02, '.',color='black', ms=4)
            
def plot_AI(dictSNV,dictCOV,final,seqDepth,genoInfo,dataType,SAMPLE,geneList,SAVE=False,chr_pairs=[[1,22],[2,21],[3,20],[4,19],[5,18],[6,17],[7,16],[8,15],[9,14],[10,13],[11,12],['X','Y']]):  
    print('***** Working on AI visualization ... *****')
    
    if geneList != 'none':      
        dfGenes = pd.read_csv(geneList, sep='\t', names=['CHR','START','END','GENE'])
        dfGenes['POS'] = (dfGenes.START + dfGenes.END)/2
    else:
        dfGenes = pd.DataFrame([], columns=['CHR','START','END','GENE','POS'])
    
    pairNum = len(chr_pairs)
    fig,axs = plt.subplots(pairNum,1,figsize=(15,pairNum*0.85))
    plt.subplots_adjust(hspace=0.3)
    axs[0].set_title(SAMPLE + '.' + dataType + '.AI', fontweight='bold')
    # plot each chr pair
    row = 0
    for i in chr_pairs:
        plot_chr_pair(dictSNV,dictCOV,final, i,axs[row], seqDepth, genoInfo, dataType, dfGenes)
        row += 1
    
    # plot legend
    if SAVE == False:
        yPos = 0.09
    else:
        yPos = 0.02
    start = 0.32
    size  = 0.02
    size2 = 0.025
    space = 0.12
    fig.add_artist(patches.Circle((0.15, yPos), 0.002, color='darkorange'))
    fig.add_artist(patches.Circle((0.25, yPos), 0.002, color='C0'))
    fig.add_artist(lines.Line2D([start+space*0, start+space*0 + size], [yPos, yPos], c='red', lw=2.5))
    fig.add_artist(lines.Line2D([start+space*1, start+space*1 + size], [yPos, yPos], c='maroon', lw=2.5))
    fig.add_artist(lines.Line2D([start+space*2, start+space*2 + size], [yPos, yPos], c='limegreen', lw=2.5))
    fig.add_artist(lines.Line2D([start+space*3, start+space*3 + size], [yPos, yPos], c='darkgreen', lw=2.5))
    fig.add_artist(lines.Line2D([start+space*4, start+space*4 + size], [yPos, yPos], c='blue', lw=2.5))
    
    fig.add_artist(Text(0.16, yPos, text='Coverage', size=11, va='center', c='darkorange'))
    fig.add_artist(Text(0.26, yPos, text='VAF', size=11, va='center', c='C0'))
    fig.add_artist(Text(start+space*0 + size2, yPos, text='Gain (with AI)', size=11, va='center', c='red'))
    fig.add_artist(Text(start+space*1 + size2, yPos, text='Gain (no AI)', size=11, va='center', c='maroon'))
    fig.add_artist(Text(start+space*2 + size2, yPos, text='Loss (with AI)', size=11, va='center', c='limegreen'))
    fig.add_artist(Text(start+space*3 + size2, yPos, text='Loss (no AI)', size=11, va='center', c='darkgreen'))
    fig.add_artist(Text(start+space*4 + size2, yPos, text='cnLOH', size=11, va='center', c='blue'))

    if SAVE == False:
        plt.show()
    else:
        if os.path.exists('results') == False:
            os.mkdir('results')
        plt.subplots_adjust(bottom=0.06, top=0.96, left=0.06, right=0.95)
        plt.savefig('results/' + SAMPLE + '.' + dataType + '.focal.png', dpi=300)
        #plt.savefig('results/' + SAMPLE + '.' + dataType + '.focal.pdf', format="pdf", dpi=300)
        

################################################################
# argparse
################################################################
class NewlineHelpFormatter(argparse.RawTextHelpFormatter):
    def _split_lines(self, text, width):
        # Preserve newlines in the help message
        if text.startswith("R|"):
            return text[2:].splitlines()
        return argparse.RawTextHelpFormatter._split_lines(self, text, width)

# Define a custom function to handle both str and int
def str_or_int(value):
    try:
        return int(value)  # Try to convert to int
    except ValueError:
        return value  # If it's not an int, return as str

if __name__ == "__main__":
    # Create an ArgumentParser
    parser = argparse.ArgumentParser(
        description="SubChrom v0.1.0: Detection of Subclonal Chromosomal aberrations in NGS data",
        formatter_class=NewlineHelpFormatter
    )

    # Define required arguments and options
    parser.add_argument("-s", "--sample", dest='samplename', type=str, required=True)
    parser.add_argument("-d", "--data",   dest='datatype', type=str, required=True, 
                       help="Sequencing data type. Options: WGS, WES, cfDNA, panel, etc")
    
    # Define optional arguments and options
    parser.add_argument("-g", dest='genome', type=str, default='hg38',
                       help="Genome build. Options: hg38 (default), hg19")
    parser.add_argument("-n", dest='PON', type=str, default='../data/cfDNA_PoN.txt',
                       help="Normal or panel of normal")
    parser.add_argument("-cs", dest='doCOVseg', type=str, default='True',
                       help="Perform coverage segmentation")
    parser.add_argument("-rs", dest='doROHseg', type=str, default='True',
                       help="Perform ROH segmentation")
    parser.add_argument("-md", dest='marker', type=str, default='/SubChrom/data/SNPmarker',
                       help="/path/to/SNPmarker Default: /SubChrom/data/SNPmarker")
    
    parser.add_argument("-minTF", dest='minTF', type=float, default=0.1,
                       help="Minimal tumor fraction to report a CNV event. Default: 0.1")
    parser.add_argument("-minSize", dest='minSize', type=int, default=1e5,
                       help="Minimal size to report a CNV event. Default: 100000")
    parser.add_argument("-minBins", dest='minBins', type=int, default=30,
                       help="Minimal data bins to report a CNV event. Default: 30")
    parser.add_argument("-gender", dest='gender', type=str, default='atuo',
                       help="""R|Gender of the sample. Options: Male/M, Female/F, auto
    Default: auto for automatic detection""") 
    parser.add_argument("-dd", dest='diploid_depth', type=str_or_int, default='auto',
                       help="""R|How to estimate the diploid depth.
    Option: chr1...chr22, chrX, a specific value (such as 500). 
    Default: auto (auto optimization)""")
    
    parser.add_argument("-covWin", dest='covWindow', type=int, default=2000000,
                       help="Coverage window size (bp) for visualization. Default: 2000000. Minimum: 500000")  
    parser.add_argument("-plotTF", dest='plotTF', type=str, default=True,
                       help="Plot tumor fraction or not. Default: True")
    parser.add_argument("-genes", dest='geneList', type=str, default='/SubChrom/data/geneList.bed',
                       help="/path/to/geneList.bed Default: /SubChrom/data/geneList.bed")
    
    # Parse command-line arguments
    args = parser.parse_args()
    
    # Access the value
    workDir = '.'
    SAMPLE  = args.samplename
    dataType = args.datatype
    
    build = args.genome
    PON = args.PON
    doCOVseg = args.doCOVseg
    doROHseg = args.doROHseg
    SNPmarker = args.marker
    
    minTF = args.minTF
    minSize = args.minSize
    minBins = args.minBins
    GENDER = args.gender
    dipDep = args.diploid_depth
    
    covWindow = args.covWindow    
    plotTF = args.plotTF
    geneList = args.geneList

    variant_path = workDir + '/' + SAMPLE + '.' + dataType + '.snp.txt'
    bedGraph = workDir + '/' + SAMPLE + '.' + dataType + '.bedGraph'
    
    # creat folder
    if os.path.exists('results') == False:
        os.mkdir('results')
    else:
        shutil.rmtree('results')
        os.mkdir('results')
    
    # load SNP and coverage raw data
    dictSNVall = SNPmarkers(variant_path, build, workDir, SNPmarker)
    dictSNV, genoInfo, rawDepth = filter_markers(dictSNVall, build, dataType, bedGraph)
    gender, chrY_status = get_gender(GENDER, genoInfo, dataType)
    dictCOV = load_coverage(dictSNV, bedGraph, genoInfo, dataType, PON)
    dfCOV = coverage_genome(dictCOV, dictSNV, genoInfo, gender, dataType, chrY_status, covWindow)
    
    # VAF segmentation
    dictMAF, dictSV, dictDfSV = segment_genome(dictSNV, gender, chrY_status, genoInfo, workDir)
    dictSmooth, dictSVs, dictDfSVs = \
        smooth_genome(dictSNV, dictMAF, dictSV, gender, dataType, minSize, minBins, build)
    
    if str(doROHseg) == 'True':
        dictMAF_ROH, dictSV_ROH, dictDfSV_ROH = \
            segment_genome(dictSNV, gender, chrY_status, genoInfo, workDir, dataCol='ROH', dataType=dataType)
        dictSmooth_ROH, dictSVs_ROH, dictDfSVs_ROH = \
            smooth_genome(dictSNV, dictMAF_ROH, dictSV_ROH, gender, dataType, 3e6, minBins, build, dataCol='ROH')
        dictDfSVs_Merged = merge_genome(dictSNV, dictMAF, dictDfSVs, dictDfSVs_ROH, minSize, minBins)
    else:
        dictDfSVs_Merged = dictDfSVs
        
    # VAF events
    vafSegments = subchrom_genome(dictSNV, dictCOV, dfCOV, dictDfSVs_Merged, genoInfo, chrY_status, dataType)
    seqDepth = diploidDepth(dictCOV, vafSegments, dataType, dipDep)
    vafSeg, rawCovSeg = vafEVENTS(vafSegments, dictSNV, genoInfo, seqDepth, gender, dataType, minTF)
    
    # coverage events and merge with VAF
    if str(doCOVseg) == 'True':
        dictMAF_COV, dictSV_COV, dictDfSV_COV = segment_genome_COV(dictCOV, gender, genoInfo, workDir, dataType, PON, doCOVseg)
        dictSmooth_COV, dictSVs_COV, dictDfSVs_COV = smooth_genome_COV(dfCOV, dictMAF_COV, dictSV_COV, dictDfSV_COV, gender, dataType, minSize, minBins)
        covSegments = covEVENTS(dictDfSVs_COV, seqDepth, gender, minBins, minTF)
        covSeg = COVminusVAF(covSegments, vafSeg, dictCOV, minSize, minBins, seqDepth, gender, minTF, dataType, build)
    else:
        covSeg = rawCovSeg
    TF, final = tumorFrac(vafSeg, covSeg, seqDepth, minTF)
    
    # plot and save
    plotFinal(dictSNV,dictCOV, dfCOV,final,genoInfo,gender,seqDepth,SAMPLE,dataType,TF,plotTF=plotTF, SAVE=True)
    plot_AI(dictSNV, dictCOV, final, seqDepth, genoInfo, dataType, SAMPLE, geneList, SAVE=True)
    
    vafSegments.to_csv('results/' + SAMPLE + '.' + dataType + '.vafSegments.txt', sep='\t', index=False)
    vafSeg.to_csv('results/' + SAMPLE + '.' + dataType + '.vafEvents.txt', sep='\t', index=False)
    final.to_csv('results/' + SAMPLE + '.' + dataType + '.SubChrom.txt', sep='\t', index=False)
    covSeg.to_csv('results/' + SAMPLE + '.' + dataType + '.covSegments.txt', sep='\t', index=False)
    
    # output sample info
    with open('results/' + SAMPLE + '.' + dataType + '.info.txt', 'w') as info:
        info.write("seqDepth:\t" + str(seqDepth) + "\n")
        info.write("gender:\t" + gender + "\n")
        info.write("chrY_status:\t" + chrY_status + "\n")
        info.write("TumorFrac:\t" + str(TF) + "\n")
