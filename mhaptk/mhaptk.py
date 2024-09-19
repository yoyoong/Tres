import math
import tabix
import pandas as pd
import matplotlib
import os

import matplotlib.pyplot as plt
import argparse, sys, os, time
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import scipy.stats  as ss
import seaborn as sns
import matplotlib.gridspec as gs
import matplotlib as mpl
import matplotlib.ticker as ticker
from pathlib import Path
import glob,re
from tqdm import tqdm
from scipy.stats import binom,chi2
import gzip
import matplotlib.gridspec as gridspec

import matplotlib.colors as mcolors
import random

matplotlib.use('Agg')
plt.switch_backend('agg')
version = '1.0'

class cpgAnnotation:
    def __init__(self, iChr, posArray):
        self.iChr = iChr
        self.posArray = posArray

class HT:
    def __init__(self, hChr, hStart, hEnd, HapMet, count, strand):
        self.hChr = hChr
        self.hStart = int(hStart)
        self.hEnd = int(hEnd)
        self.HapMet = HapMet
        self.count = int(count)
        self.WC = strand

class Stat:
    def __init__(self,
                 mhap_path: str,
                 cpg_path: str,
                 maxK: int,
                 minK: int,
                 stranded: str,
                 K: int):
        self.mhap_path = mhap_path
        self.cpg_path = cpg_path
        self.stranded = stranded
        self.K = K
        self.startK = minK
        self.endK = maxK
        self.data = tabix.open(mhap_path)
        self.statslist = ["MM", 'CHALM', 'PDR', 'MHL', 'MBS', 'MCR', 'Entropy', 'R2']
        self.statsbool = {"MM": False, 'CHALM': False, 'PDR': False, 'MHL': False, 'MBS': False, 'MCR': False,
                          'Entropy': False}
    def getBed(self,bed_path, line=0):
        self.bed = pd.read_csv(bed_path, sep='\t', header=None).iloc[line, :]
        self.start = int(self.bed[1])
        self.end = int(self.bed[2])
        self.Chr = self.bed[0]
        self.len = self.end - self.start
        self.name = str(self.bed[0]) + ':' + str(self.bed[1]) + '-' + str(self.bed[2])
    def Region(self,string):
        self.Chr = string.split(':')[0]
        self.start = int(string.split(':')[1].split('-')[0])-1
        self.end = int(string.split(':')[1].split('-')[1])
        self.len= self.end - self.start
        self.name = self.Chr + ":" + str(self.start) + "-" + str(self.end)

    def getrecord(self, strand: str = 'both'):

        self.record = []
        self.strand = strand

        try:
            records = self.data.query(self.Chr, self.start, self.end)
        except:
        #     print('record query error')
            print(self.Chr, self.start, self.end, 'the region can not be loaded')
            return self.record
        for rd in records:

            if strand == 'plus' and rd[5] != '+':
                continue
            if strand == 'minus' and rd[5] != '-':
                continue

            self.record.append(HT(rd[0], rd[1], rd[2], rd[3], rd[4], rd[5]))

        return self.record

    def aimming(self, stats_):

        for k in self.statslist:
            if k in stats_:
                self.statsbool[k] = True
            else:
                self.statsbool[k] = False
        #print('config reset as: ', self.statsbool)
        return self.statsbool
    def Kmer_gen(self, s, k=4):  # r:str
        '''
        To cut the kmers from a long sequense.
        - input:
         - s: input sequence, type : str
         - k: kmer length , type:int , default=k=4
        - output
         - ret_ : return as a bag of kmer(string type) in a list.
        '''
        ret_ = []
        for i in range(len(s) - k + 1):
            ret_.append(s[i: i + k])
        return ret_
    def calculating(self):
        self.kmer = self.K  # kmer片段长度
        self.result = {}  # 结果dictionay

        for a in self.statsbool.keys():
            if self.statsbool[a] == True:
                self.result[a] = ''

        if self.record ==[]:
            return self.result
        tBase = 0
        mBase = 0
        K4plus = 0
        nDR = 0
        nMR = 0
        cBase = 0
        nReads = 0
        self.Total_dict = {}
        self.epiallell_dic = {}

        if self.record != []:
            for val in self.record:
                nReads += val.count
                tBase += val.count * len(val.HapMet)
                mBase += val.HapMet.count('1') * val.count
                if '1' in val.HapMet:
                    cBase += val.HapMet.count('0') * val.count
                if len(val.HapMet) >= self.kmer:
                    K4plus += val.count
                    if '1' in val.HapMet:
                        nMR += val.count

                    if len(set(val.HapMet)) > 1:
                        nDR += val.count
                if self.statsbool['Entropy']:
                    if len(val.HapMet) >= self.kmer:
                        rk_ = self.Kmer_gen(val.HapMet, self.kmer)
                        for rk in rk_:
                            try:
                                self.epiallell_dic[rk] += val.count
                            except:
                                self.epiallell_dic[rk] = val.count
            self.nReads, self.mBase, self.tBase, self.K4plus, self.nDR, self.nMR, self.cBase = nReads, mBase, tBase, K4plus, nDR, nMR, cBase

        else:
            self.nReads, self.mBase, self.tBase, self.K4plus, self.nDR, self.nMR, self.cBase=0, 0, 0, 0, 0, 0,0
        '''
        ——————————按照statsbool的情况计算统计量，要分情况计算————————————
        '''
        if self.statsbool['PDR']:
            if self.K4plus == 0:
                self.result['PDR'] = ''
            else:
                self.result['PDR'] = self.nDR /self.K4plus #if self.K4plus != 0 else 0
        if self.statsbool['MHL']:
            cpgAnn = self.tabixCPG(self.cpg_path)
            self.buildBinaryMatrix()
            self.result['MHL'] = self.getMHL()

        if self.statsbool['MM']:  # 简单
            if self.tBase ==0:
                self.result['MM'] = ''
            else:
                self.result['MM'] = self.mBase / self.tBase #if self.tBase != 0 else 0 # 甲基化位点数/所有位点数
        if self.statsbool['CHALM']:  # 简单
            if self.K4plus == 0:
                self.result['CHALM'] = ''
            else:
                self.result['CHALM'] = self.nMR / self.K4plus #if self.K4plus != 0 else 0 # 含甲基化的片段比例
        if self.statsbool['MBS']:
            self.result['MBS'] = self.get_mark_MBS(self.record)
        if self.statsbool['MCR']:
            if self.tBase ==0:
                self.result['MCR'] = ''
            else:
                self.result['MCR'] = self.cBase/self.tBase #if self.tBase != 0 else 0
        if self.statsbool['Entropy']:
            temp = 0
            all = sum(self.epiallell_dic.values())
            for i in self.epiallell_dic.values():
                temp += i/all*np.log2(i/all)

            self.result['Entropy'] = abs(-1/4 * temp)


        return self.result

    def Kmer(self, s):
        ret_ = []
        for i in range(len(s) - self.K + 1):
            ret_.append(s[i: i + self.K])
        return ret_

    def tabixCPG(self, htGZ, shift=500):
        self.cpg_path = htGZ
        tb = tabix.open(htGZ)
        records = tb.query(self.Chr, self.start - shift, self.end + shift)
        posArray = []
        for rd in records:
            if len(rd) < 3:
                continue
            posArray.append(int(rd[1]))
        self.cpgAnn = cpgAnnotation(self.Chr, posArray)

        return self.cpgAnn

    def get_r_MBS(self, r):
        Li = len(r)
        r_split = r.split("0")
        res = 0
        for lij in r_split:
            res = res + len(lij) ** 2
        res = res / Li ** 2
        return res

    def get_mark_MBS(self, r, cut=4):
        res = 0
        n_r = 0
        for i in r:
            if len(i.HapMet) < cut:
                continue
            res = res + self.get_r_MBS(i.HapMet) * int(i.count)
            n_r = n_r + int(i.count)
        try:
            res = res / n_r
            return res
        except:
            return  -1

    def getMHL(self):
        Nc = np.shape(self.MC)[1]
        obsK = max(np.sum(self.MC, 1))
        maxK = min([obsK, self.endK])
        if (self.startK > maxK):
            print("Error: startK is too large.\n")
            sys.exit(-1)

        uCount = np.zeros(maxK, dtype='int')
        mCount = np.zeros(maxK, dtype='int')
        tCount = np.zeros(maxK, dtype='int')

        for k in range(maxK):
            for j in range(Nc - k):

                xM0 = self.M0[..., j:(j + k + 1)]
                xM1 = self.M1[..., j:(j + k + 1)]
                xMC = self.MC[..., j:(j + k + 1)]
                uCount[k] += sum(((np.sum(xM0, 1) == (k + 1)) + 0) * self.count)
                mCount[k] += sum(((np.sum(xM1, 1) == (k + 1)) + 0) * self.count)
                tCount[k] += sum(((np.sum(xMC, 1) == (k + 1)) + 0) * self.count)

        mHL = sumFactor = 0.0
        for k in range(self.startK - 1, maxK):
            mHL += (k + 1.0) * mCount[k] / tCount[k]
            sumFactor += k + 1.0
        mHL = round(mHL / sumFactor * 100000) / 100000.0

        return mHL

    def info_to_file(self, k=4):
        tBase = 0
        mBase = 0
        K4plus = 0
        nDR = 0
        nMR = 0
        cBase = 0
        nReads = 0

        if self.record != []:
            for val in self.record:
                nReads += val.count
                tBase += val.count * len(val.HapMet)
                mBase += val.HapMet.count('1') * val.count
                if '1' in val.HapMet:
                    cBase += val.HapMet.count('0') * val.count
                if len(val.HapMet) >= k:
                    K4plus += val.count
                    if '1' in val.HapMet:
                        nMR += val.count

                    if len(set(val.HapMet)) > 1:
                        nDR += val.count

            self.nReads, self.mBase, self.tBase, self.K4plus, self.nDR, self.nMR, self.cBase = nReads, mBase, tBase, K4plus, nDR, nMR, cBase
            df = pd.DataFrame([self.Chr, self.start, self.end, nReads, mBase, tBase, K4plus, nDR, nMR]).T
            df.columns = ['chr', 'start', 'end', 'nReads', 'mBase', 'tBase', 'K4plus', 'nDR', 'nMR']
        else:
            self.nReads, self.mBase, self.tBase, self.K4plus, self.nDR, self.nMR = 0, 0, 0, 0, 0, 0
            df = pd.DataFrame([self.Chr, self.start, self.end, 0, 0, 0, 0, 0, 0]).T
            df.columns = ['chr', 'start', 'end', 'nReads', 'mBase', 'tBase', 'K4plus', 'nDR', 'nMR']

        return df

    def buildBinaryMatrix(self):
        posArray = self.cpgAnn.posArray
        posDict = dict()
        Nc = len(posArray)
        Nr = len(self.record)

        for i in range(Nc):
            posDict[self.cpgAnn.iChr + ":" + str(posArray[i])] = i

        self.MC = np.zeros((Nr, Nc), dtype='int')
        self.M0 = np.zeros((Nr, Nc), dtype='int')
        self.M1 = np.zeros((Nr, Nc), dtype='int')

        self.count = np.zeros(Nr, dtype='int')
        strand = np.zeros(Nr, dtype='int')

        for i in range(Nr):
            ht = self.record[i]
            pos = ht.hChr + ":" + str(ht.hStart)

            self.count[i] = ht.count
            if ht.WC == "+":
                strand[i] = 1
            # print(pos, posDict)
            if pos not in posDict:
                print("Haplotype positions are out of range.")
                sys.exit(-1)
            Idx = posDict[pos]
            HapMet = ht.HapMet
            for j in range(len(HapMet)):
                self.MC[i, Idx + j] = 1
                if HapMet[j] == "0":
                    self.M0[i, Idx + j] = 1
                elif HapMet[j] == "1":
                    self.M1[i, Idx + j] = 1
                else:
                    print("Haplotypes must be 0/1 only.")
                    sys.exit(-1)
        return [self.MC, self.M0, self.M1, self.count, strand]

    def count2pos(self, MC, M0, M1, count, pi, pj):
        if pi >= np.shape(MC)[1] or pj >= np.shape(MC)[1]:
            print("pi or pj are out of range.")
            sys.exit(-1)

        iAll = MC[..., pi] + MC[..., pj] == 2
        iN00 = M0[..., pi] + M0[..., pj] == 2
        iN01 = M0[..., pi] + M1[..., pj] == 2
        iN10 = M1[..., pi] + M0[..., pj] == 2
        iN11 = M1[..., pi] + M1[..., pj] == 2

        N00 = sum(count[np.logical_and(iAll, iN00)])
        N01 = sum(count[np.logical_and(iAll, iN01)])
        N10 = sum(count[np.logical_and(iAll, iN10)])
        N11 = sum(count[np.logical_and(iAll, iN11)])

        iAlli = MC[..., pi] == 1
        iAllj = MC[..., pj] == 1
        iNi0 = M0[..., pi] == 1
        iNi1 = M1[..., pi] == 1
        iNj0 = M0[..., pj] == 1
        iNj1 = M1[..., pj] == 1

        Ni0 = sum(count[np.logical_and(iAlli, iNi0)])
        Ni1 = sum(count[np.logical_and(iAlli, iNi1)])
        Nj0 = sum(count[np.logical_and(iAllj, iNj0)])
        Nj1 = sum(count[np.logical_and(iAllj, iNj1)])
        return [N00, N01, N10, N11], [Ni0, Ni1, Nj0, Nj1]

    def get_r2(self, N00, N01, N10, N11, Ni0, Ni1, Nj0, Nj1, cov=10):

        n = N00 + N01 + N10 + N11

        if n < cov:
            try:
                Pi = Ni1 / (Ni1 + Ni0)
                Pj = Nj1 / (Nj1 + Nj0)
                return (np.nan, np.nan)
            except:
                return (np.nan, np.nan)

        Pa = (N10 + N11) / n
        Pb = (N01 + N11) / n
        Pi = Ni1 / (Ni1 + Ni0)
        Pj = Nj1 / (Nj1 + Nj0)

        if (Pa * (1 - Pa) * Pb * (1 - Pb) == 0):
            return (0, 0.5)

        r = (N11 / n - Pa * Pb) / np.sqrt(Pa * (1 - Pa) * Pb * (1 - Pb))
        if r >= 1:
            return (r ** 2, 1)
        if r <= -1:
            return (r ** 2, 0)
        t = r * np.sqrt((n - 1) / (1 - r ** 2))

        pval = ss.t.cdf(t, n - 1)
        return (r ** 2, pval)

def main():

    parser = argparse.ArgumentParser()

    subparsers1 =  parser.add_subparsers()
    gw = subparsers1.add_parser('genomeWide', help='calculate methylation metrics for mHaps that cover each CpG site across the genome')
    gw.add_argument('--tag', type=str, required=True, help='prefix of the output file(s)')
    gw.add_argument('--mhapPath', type=str,required=True, help='input file, mhap.gz format, generated by mHapTools and indexed')
    gw.add_argument('--cpgPath', type=str, required=True, help='genomic CpG file, gz format and indexed')
    gw.add_argument('--metrics', required=True, nargs='*',help='mHap-level metrics, including MM, PDR, CHALM, MHL, MCR, MBS, Entropy, and R2')
    gw.add_argument('--outputDir', type=str, required=True, help='output directory, created in advance')
    gw.add_argument('--minK',  default=np.inf,type=int, help='minimum k-mer length for MHL [1]')
    gw.add_argument('--maxK',default=np.inf,type=int, help='maximum k-mer length for MHL [10]')
    gw.add_argument('--K',  default=np.inf,type=int, help='k-mer length for entropy, PDR, and CHALM, can be 3, 4, or 5 [4]')
    gw.add_argument('--strand', type=str, default='both', help='strand information, one of plus, minus and both [both]')
    gw.set_defaults(func='genomeWide')



    stat = subparsers1.add_parser('stat', help='calculate methylation metrics for mHaps that cover predefined regions')
    stat.add_argument('--metrics', nargs='*',default=None, help='mHap-level metrics, including MM, PDR, CHALM, MHL, MCR, MBS, and Entropy [None]')
    stat.add_argument('--mhapPath', type=str, required=True, help='input file, mhap.gz format, generated by mHapTools and indexed')
    stat.add_argument('--cpgPath', type=str, required=True, help='genomic CpG file, gz format and Indexed')
    stat.add_argument('--region', type=str, default=None, help='one region, in the format of chr:start-end')
    stat.add_argument('--bedPath', type=str, default=None, help='input BED file')
    stat.add_argument('--outputFile', type=str, required=True, help='output file name')
    stat.add_argument('--minK',  default=np.inf,type=int, help='minimum k-mer length for MHL [1]')
    stat.add_argument('--maxK', default=np.inf,type=int, help='maximum k-mer length for MHL [10]')
    stat.add_argument('--K',  default=np.inf,type=int,
                    help='k-mer length for entropy, PDR, and CHALM, can be 3, 4, or 5 [4]')
    stat.add_argument('--strand', type=str, default='both', help='strand information, one of plus, minus and both [both]')
    stat.set_defaults(func='stat')


    tanghulu = subparsers1.add_parser('tanghulu', help='plot the DNA methylation status for mHaps in a region')
    tanghulu.add_argument('--mhapPath', type=str,required=True ,help='input file, mhap.gz format, generated by mHapTools and indexed')
    tanghulu.add_argument('--simulation',action='store_true',help='indicates whether mHaps should be simulated')
    tanghulu.add_argument('--cpgPath', type=str,  required=True, help='genomic CpG file, gz format and Indexed')
    tanghulu.add_argument('--region', type=str, required=True, help='one region, in the format of chr:start-end')
    tanghulu.add_argument('--merge', action='store_true', help='indicates whether identical mHaps should be merged')
    tanghulu.add_argument( "--outcut", type=int, default=2000, help="the max length of region to plot, default is 2000")
    tanghulu.add_argument('--outputFile', type=str,required=True, help='output file name')
    tanghulu.set_defaults(func='tanghulu')

    R2 = subparsers1.add_parser('R2', help='calculate linkage disequilibrium between CpG sites within predefined regions')
    R2.add_argument('--tag', type=str, required=True, help='prefix of the output file(s)')
    R2.add_argument('--mhapPath', type=str, required=True, help='input file, mhap.gz format, sorted by samtools')
    R2.add_argument('--cpgPath', type=str, required=True, help='genomic CpG file, gz format and Indexed')
    R2.add_argument('--region', type=str, required=True, help='one region, in the format of chr:start-end')
    R2.add_argument('--outputDir', type=str, required=True, help='output directory name')
    R2.add_argument('--mHapView',action='store_true',  help='plot linkage disequilibrium patterns of pair-wise CpGs')
    R2.add_argument('--longrange', action='store_true', help='indicates whether generate a file in longrange format')
    R2.add_argument('--strand', type=str, default='both', help='strand information, one of plus, minus and both [both]')
    R2.set_defaults(func='R2')


    MHBDiscovery = subparsers1.add_parser('MHBDiscovery', help='identification of methylation haplotype blocks within a region or genome-wide')
    MHBDiscovery.add_argument('--mhapPath', type=str,required=True, help='input file, mhap.gz format, generated by mHapTools and indexed')
    MHBDiscovery.add_argument('--cpgPath', type=str, required=True, help='genomic CpG file, gz format and Indexed')
    MHBDiscovery.add_argument('--region', type=str, help='one region, in the format of chr:start-end')
    MHBDiscovery.add_argument('--bedPath', type=str,default=None, help='input BED file')
    MHBDiscovery.add_argument('--outputFile', type=str, required=True, help='output file name')
    MHBDiscovery.add_argument("--window", type=int, default=5, required=False, help="size of core window [5]")
    MHBDiscovery.add_argument("--r_square", type=float, default=0.5, required=False, help="R-square cutoff [0.5]")
    MHBDiscovery.add_argument("--p_value", type=float, default=0.05, required=False, help="P-value cutoff [0.05]")
    MHBDiscovery.set_defaults(func='MHBDiscovery')

    args = parser.parse_args()
    try:
        args.func
    except:
        print('mHapTk:A comprehensive tool kit for analysis of DNA methylation haplotypes')
        print('version:', version)
        sys.exit()

    if args.func == 'genomeWide':
        assert args.strand == 'both' or args.strand == 'plus' or args.strand == 'minus', '--stranded should be both plus minus'

        if 'MHL' not in args.metrics and (args.maxK != np.inf or args.minK != np.inf):
            print('Warning: --maxK and --minK is only for mhl')
        if not('PDR' in args.metrics or 'CHALM' in args.metrics or 'Entropy' in args.metrics) and args.K != np.inf:
            print('Warning: --K is only for PDR CHALM Entropy')
        if args.maxK == np.inf:
            args.maxK = 10
        if args.minK == np.inf:
            args.minK = 1
        if args.K == np.inf:
            args.K = 4
        assert isinstance(args.maxK, int), 'maxK should be int'
        assert isinstance(args.K, int), 'K should be int'
        assert isinstance(args.minK, int), 'minK should be int'
        assert Path(args.outputDir).exists(), 'outDir does not exist'
        assert 3<= args.K <=5, 'K：the default is 4 and values must be between 3 and 5'
        assert  args.maxK > args.minK, 'maxK should be larger than minK'
        assert 1 <= args.maxK <= 10, 'maxK should be in 1 to 10'
        assert 1 <= args.minK <= 10, 'minK should be in 1 to 10'
        resultPath = args.outputDir + '/' + args.tag + '_'
        print('the stat u chose is ', args.metrics)
        gw = GenomeWide(args.mhapPath,
                        args.cpgPath,
                        args.maxK,
                        args.minK,
                        args.strand,
                        args.K)
        for stat in args.metrics:
            print(stat)
            if stat not in ['MM', 'CHALM', 'PDR', 'MHL', 'MBS', 'MCR', 'Entropy', 'R2']:
                print('you input a wrong stat')
                print('the right input like', gw.statslist)
            else:
                if stat == 'MM':
                    Time = time.time()
                    MM = gw.MM()
                    print('MM was done')
                    # MM.to_csv(resultPath + 'MM GW.csv', sep='\t', index=False, header=None)
                    MM_new = pd.concat([MM.iloc[:, 0], MM.iloc[:, 1] - 1, MM.iloc[:, 1], MM.iloc[:, -1]], axis=1).round(8)
                    MM_new[(1 - np.isnan(MM_new.iloc[:, 3])).astype(np.bool_)].to_csv(resultPath + 'MM.bedGraph',
                                                                                      index=False, header=None,
                                                                                      sep='\t')
                    print('MM time span:', time.time() - Time)
                if stat == 'CHALM':
                    Time = time.time()
                    CHALM = gw.CHALM()
                    print('CHALM was done')
                    CHALM_new = pd.concat(
                        [CHALM.iloc[:, 0], CHALM.iloc[:, 1] - 1, CHALM.iloc[:, 1], CHALM.iloc[:, -1]], axis=1).round(8)
                    CHALM_new[(1 - np.isnan(CHALM_new.iloc[:, 3])).astype(np.bool)].to_csv(
                        resultPath + 'CHALM.bedGraph', index=False, header=None, sep='\t')
                    print('CHALM time span:', time.time() - Time)
                if stat == 'PDR':
                    Time = time.time()
                    PDR = gw.PDR()
                    print('PDR was done ')
                    # PDR.to_csv(resultPath + 'PDR GW.csv', sep='\t', index=False, header=None)
                    PDR_new = pd.concat([PDR.iloc[:, 0], PDR.iloc[:, 1] - 1, PDR.iloc[:, 1], PDR.iloc[:, -1]],
                                        axis=1).round(8)
                    PDR_new[(1 - np.isnan(PDR_new.iloc[:, 3])).astype(np.bool_)].to_csv(
                        resultPath + 'PDR.bedGraph', index=False, header=None, sep='\t')
                    print('PDR time span:', time.time() - Time)
                if stat == 'MHL':
                    Time = time.time()
                    MHL = gw.MHL()
                    print('MHL was done')
                    # MHL.to_csv(resultPath + 'MHL GW.csv', sep='\t', index=False, header=None)
                    MHL_new = pd.concat([MHL.iloc[:, 0], MHL.iloc[:, 1] - 1, MHL.iloc[:, 1], MHL.iloc[:, -1]],
                                        axis=1).round(8)
                    MHL_new[(1 - np.isnan(MHL_new.iloc[:, 3])).astype(np.bool_)].to_csv(
                        resultPath + 'MHL.bedGraph', index=False, header=None, sep='\t')
                    print('MHL time span:', time.time() - Time)
                if stat == 'MBS':
                    Time = time.time()
                    MBS = gw.MBS()
                    print('MBS was done ')
                    # MBS.to_csv(resultPath + 'MBS GW.csv', sep='\t', index=False, header=None)
                    MBS_new = pd.concat([MBS.iloc[:, 0], MBS.iloc[:, 1] - 1, MBS.iloc[:, 1], MBS.iloc[:, -1]],
                                        axis=1).round(8)
                    MBS_new[(1 - np.isnan(MBS_new.iloc[:, 3])).astype(np.bool_)].to_csv(
                        resultPath + 'MBS.bedGraph', index=False, header=None, sep='\t')
                    print('MBS time span:', time.time() - Time)
                if stat == 'MCR':
                    Time = time.time()
                    MCR = gw.MCR()
                    print('MCR was done')
                    # MCR.to_csv(resultPath + 'MCR GW.csv', sep='\t', index=False, header=None)
                    MCR_new = pd.concat([MCR.iloc[:, 0], MCR.iloc[:, 1] - 1, MCR.iloc[:, 1], MCR.iloc[:, -1]],
                                        axis=1).round(8)
                    MCR_new[(1 - np.isnan(MCR_new.iloc[:, 3])).astype(np.bool_)].to_csv(
                        resultPath + 'MCR.bedGraph', index=False, header=None, sep='\t')
                    print('MCR time span:', time.time() - Time)
                if stat == 'Entropy':
                    Time = time.time()
                    Entropy = gw.Entropy()
                    print('Entropy was done')
                    # Entropy.to_csv(resultPath + 'Entropy GW.csv', sep='\t', index=False, header=None)
                    Entropy_new = pd.concat(
                        [Entropy.iloc[:, 0], Entropy.iloc[:, 1] - 1, Entropy.iloc[:, 1], Entropy.iloc[:, -1]],
                        axis=1).round(8)
                    Entropy_new[(1 - np.isnan(Entropy_new.iloc[:, 3])).astype(np.bool_)].to_csv(
                        resultPath + 'Entropy.bedGraph', index=False, header=None, sep='\t')
                    print('Entropy time span:', time.time() - Time)
                if stat == 'R2':
                    Time = time.time()
                    R2 = gw.R2()
                    print('R2 was done')
                    # R2.to_csv(resultPath + 'R2 GW.csv', sep='\t', index=False, header=None)
                    # R2_new = pd.concat([R2.iloc[:, 0], R2.iloc[:, 1] - 1, R2.iloc[:, 1], R2.iloc[:, -1]], axis=1)
                    # print(R2_new)
                    # print(np.isnan(R2.iloc[:, 3]))
                    R2_new = pd.concat([R2.iloc[:, 0], R2.iloc[:, 1] - 1, R2.iloc[:, 1], R2.iloc[:, -1]], axis=1)
                    R2_new = R2_new.dropna().round(8)
                    R2_new.to_csv(resultPath + 'R2.bedGraph', index=False, header=None, sep='\t')
                    # R2_new[(1 - np.isnan(R2.iloc[:, 3])).astype(np.bool_)].to_csv(
                    #     resultPath + 'R2.bedGraph', index=False, header=None, sep='\t')
                    print('R2 time span:', time.time() - Time)

    if args.func == 'R2':
        assert Path(args.outputDir).exists(), 'outDir does not exist'
        assert args.strand == 'both' or args.strand == 'plus' or args.strand == 'minus', '--stranded should be both plus minus'

        resultPath = args.outputDir + '/' + args.tag + '_'
        M = R2_c(args.mhapPath,
                     args.cpgPath,
                     args.strand)
        M.Region(args.region)
        M.tabixCPG(args.cpgPath, shift=500)

        M.getrecord(strand=args.strand)

        position = M.Chr + "_" + str(M.start + 1) + "_" + str(M.end)
        cpgAnn = M.tabixCPG(args.cpgPath)
        [MC, M0, M1, count, strand_] = M.buildBinaryMatrix()

        samplename = os.path.basename(args.mhapPath)
        outtxt = resultPath + position + ".cpg_sites_rsquare.txt"
        outpdf1 = resultPath + position + ".cpg_sites_rsquare.pdf"
        outpdf2 = resultPath + position + ".cpg_sites_rsquare_hp.pdf"
        if args.longrange:
            outlongrange = resultPath  + M.Chr + '_' + str( M.start) + '_' + str(M.end) + '_' +'longrange'


        sites_pos = np.where(sum(MC) != 0)[0]

        xx_matrix = np.zeros([len(sites_pos), len(sites_pos)], dtype=int)
        pval_matrix = np.zeros([len(sites_pos), len(sites_pos)], dtype=float)
        ref_posArray = [cpgAnn.posArray[int(i)] for i in np.where(sum(MC) != 0)[0]]
        R_matrix = np.zeros([len(sites_pos), len(sites_pos)], dtype=float)
        with open(outtxt, 'w+') as f:
            if args.longrange:
                f2 = open(outlongrange, 'w+')
            f.write("\t".join(['chr', 'posi', 'posj','N00', 'N01', 'N10', 'N11', 'Chr', 'r2', 'pvalue']))
            f.write("\n")
            for i in range(0, len(sites_pos)):
                pi = sites_pos[i]

                if not M.start < ref_posArray[i] <= M.end:
                    continue
                for j in [j for j in range(0, len(sites_pos)) if sites_pos[j] >= sites_pos[i]]:
                    if not M.start < ref_posArray[j] <= M.end:
                        continue
                    pj = sites_pos[j]


                    (N00, N01, N10, N11), (Ni0, Ni1, Nj0, Nj1) = M.count2pos(MC, M0, M1, count, pi, pj)
                    r2, pval = M.D_r2(N11, N10, N01, N00)
                    R_matrix[i][j] = r2
                    pval_matrix[i][j] = pval
                    xx_matrix[i][j] = 1
                    if i != j:
                        f.write(
                            "\t".join([cpgAnn.iChr, str(ref_posArray[i]), str(ref_posArray[j]), str(N00), str(N01), str(N10), str(N11),
                                        str(round(r2,8)), str(pval)]))
                        if args.longrange:
                            f2.write(f'{cpgAnn.iChr}\t{ref_posArray[i]}\t{ref_posArray[i]+1}\t{cpgAnn.iChr}:{ref_posArray[j]}-{ref_posArray[j]+1},{round(r2,8)}')
                            f2.write('\n')
                        f.write("\n")

        rd1 = []
        rd2 = []
        for i in range(R_matrix.shape[0]):
            if sum(R_matrix[i]) != 0:
                rd1.append(i)
            if sum(R_matrix[:, i]) != 0:
                rd2.append(i)

        R_matrix = R_matrix[rd1][:, rd2]
        xx_matrix = xx_matrix[rd1][:, rd2]

        M.paint_rsquare_plot(samplename, R_matrix, xx_matrix, MC, outpdf1)
        M.paint_rsquare_heatmap( R_matrix, MC, outpdf2)

        if args.mHapView:
            if args.cpgPath is None:
                print('lack of cpg path')
            else:
                M.tabixCPG(args.cpgPath)
                hp_cpg = M.hp_CPG(resultPath)
            outpdf_join = resultPath + position + "_join.pdf"
            M.join_pic( R_matrix, xx_matrix, MC, hp_cpg, outpdf_join)

    if args.func == 'stat':
        if (args.metrics is not None) and ('MHL' not in args.metrics)  and (args.maxK != np.inf or args.minK != np.inf):
            print('Warning: --maxK and --minK is only for mhl')
        if (args.metrics is not None) and not('PDR' in args.metrics or 'CHALM' in args.metrics or 'Entropy' in args.metrics) and args.K != np.inf:
            print('Warning: --K is only for PDR CHALM Entropy')
        if args.maxK == np.inf:
            args.maxK = 10
        if args.minK == np.inf:
            args.minK = 1
        if args.K == np.inf:
            args.K = 4
        assert isinstance(args.maxK, int), 'maxK should be int'
        assert isinstance(args.K, int), 'K should be int'
        assert isinstance(args.minK, int), 'minK should be int'
        assert Path(os.path.dirname(args.outputFile)).exists(), 'outDir does not exist'
        assert  args.maxK > args.minK, 'maxK should be larger than minK'
        assert 3 <= args.K <= 5, 'K：the default is 4 and values must be between 3 and 5'
        assert  args.region is  None or args.bedPath is  None, 'U should only inuput bedPath or region'
        assert (args.region is not None) or (args.bedPath is not None), 'U should input bedPath or region'
        assert 1<= args.maxK <=10, 'maxK should be in 1 to 10'
        assert 1 <= args.minK <= 10, 'minK should be in 1 to 10'
        assert args.strand == 'both' or args.strand == 'plus' or args.strand == 'minus', '--strand should be both plus minus'

        resultPath = args.outputFile
        M = Stat(args.mhapPath,
                        args.cpgPath,
                        args.maxK,
                        args.minK,
                        args.strand,
                        args.K)
        if args.bedPath is not None:

            lines = pd.read_csv(args.bedPath, sep='\t', header=None).shape[0]

            f_info = open(resultPath, 'w+')
            f_info.write('chr\tstart\tend\tnReads\tmBase\tcBase\ttBase\tK4plus\tnDR\tnMR')

            if args.metrics:
                if 'MHl' in args.metrics:
                    cpgAnn = M.tabixCPG(args.cpgPath)
                    M.buildBinaryMatrix()
                for i in args.metrics:
                    f_info.write('\t')
                    f_info.write(i)
            f_info.write('\n')

            for i in tqdm(range(lines), desc='Read bed lines'):
                M.getBed(args.bedPath, i)
                M.getrecord(strand=args.strand)
                # --------------------把mhap的数据计算储存---------------------------------------------
                # columns = ['chr', 'start', 'end', 'nReads', 'mBase', 'tBase', 'K4plus', 'nDR', 'nMR']
                M.info_to_file()
                f_info.write(f'{M.Chr}\t{M.start}\t{M.end}\t{M.nReads}\t{M.mBase}\t{M.cBase}\t{M.tBase}\t{M.K4plus}\t{M.nDR}\t{M.nMR}')


                # ---------------------选择计算的统计量--------------------------------------------

                if args.metrics is not None:
                    stats_list = []
                    for stats in args.metrics:
                        stats_list.append(stats)
                    M.aimming(stats_list)
                    if 'MHL' in args.metrics:
                        dic_stat = M.calculating()
                    else:
                        dic_stat = M.calculating()

                    for i, val in enumerate(args.metrics):

                            f_info.write(f'\t{round(dic_stat[val],8)}')

                f_info.write('\n')
            f_info.close()

        elif args.region is not None:
            M.Region(args.region)
            M.getrecord(args.strand)
            M.info_to_file()

            f = open(resultPath, 'w+')
            f.write('chr\tstart\tend\tnReads\tmBase\tcBase\ttBase\tK4plus\tnDR\tnMR')
            if args.metrics:
                stats_list = []
                for stats in args.metrics:
                    stats_list.append(stats)
                M.aimming(stats_list)
                dic_stat = M.calculating()
                print(dic_stat)
                for key in dic_stat:
                    f.write('\t' + key)

            f.write('\n')
            f.write(f'{M.Chr}\t{M.start}\t{M.end}\t{M.nReads}\t{M.mBase}\t{M.cBase}\t{M.tBase}\t{M.K4plus}\t{M.nDR}\t{M.nMR}')
            if args.metrics:
                for i,key in enumerate(dic_stat):
                        f.write(f'\t{round(dic_stat[key],8)}')
            f.write('\n')

            f.close()

        else:
            print('you should input the region you need if you dont know what to input.Please use help.')

    if args.func == 'tanghulu':
        assert Path(os.path.dirname(args.outputFile)).exists(), 'outDir does not exist'
        resultPath = args.outputFile
        M = Tanghulu(args.mhapPath,args.cpgPath)
        M.Region(args.region)
        M.getrecord()
        M.tabixCPG(args.cpgPath, shift=500)

        position = M.Chr + "_" + str(M.start + 1) + "_" + str(M.end)
        if not args.simulation:
            if M.len > args.outcut:
                print("The region is larger than " + str(
                    args.outcut) + ", it's not recommanded to do tanghulu plotting and system will exit right now...")
                time.sleep(0.1)
                sys.exit()
            [MC, M0, M1, count, strand_] = M.buildBinaryMatrix()
            samplename = os.path.basename(args.mhapPath)

            outpdf = resultPath


            M.paint_tanghulu_plot(args, MC, M0, M1, count, strand_, samplename, outpdf)
        else:
            a, _, _, _,ref_posArray = M.simulate()
            mm = M.MM()
            plot_x = []
            plot_y = []
            meth = []

            for i in range(a.shape[0]):
                for j in range(a.shape[1]):
                    plot_x.append(j)
                    plot_y.append(i)
                    meth.append(a[i, j])

            plt.clf()
            fig = plt.gcf()
            plt.figure(dpi=200, figsize=(8, 15))
            xinches = 10
            yinches = 6 / 24 * 40
            ax = plt.axes([0.1, 0.1, .7, .8])
            plt.title(f'Average methylation:{mm}', fontdict={'size': 22})
            for x, y, m in zip(plot_x, plot_y, meth):
                inner = 'w' if m == 0 else 'k'
                ax.scatter(x, y, c=inner, linewidths=1, edgecolors='k', s=800, zorder=2)
            for i in range(a.shape[0]):
                for j in range(a.shape[1] - 1):
                    ax.plot([j, j + 1], [i, i], zorder=1, c='k')
            for i in range(a.shape[1]):
                ax.text(i , -2, ref_posArray[i], rotation=60, fontdict={'size': 13},horizontalalignment='center')
            ax.set_xticks([])
            ax.set_yticks([])
            plt.savefig(resultPath, dpi=200)

    if args.func == "MHBDiscovery":
        assert Path(os.path.dirname(args.outputFile)).exists(), 'outDir does not exist'
        Time = time.time()
        resultPath = args.outputFile
        M = MHB()
        irList = []
        if args.region is not None:
            chunks = re.split(':|-', args.region)
            queryIR = IR(chunks[0], chunks[1], chunks[2], '')
            irList.append(queryIR)
        else:
            irList = M.loadIR(args.bedPath, BED=True)

        OUT = open(resultPath , "w")
        for queryIR in irList:
            htList = M.tabixHT(args.mhapPath, queryIR)
            if len(htList) == 0:
                continue
            cpgAnn = M.tabixCPG(args.cpgPath, queryIR, shift=500)
            [MC, M0, M1, count, strand] = M.buildBinaryMatrix(htList, cpgAnn)
            irHB = M.getMHB(MC, M0, M1, count, cpgAnn, window=args.window, r_square=args.r_square, p_value=args.p_value)
            for ir in irHB:
                outString = '\t'.join([ir.iChr, str(ir.iStart), str(ir.iEnd)])
                # print(outString)
                OUT.write(outString + "\n")
        OUT.close()
        print('MHB time span:', time.time() - Time)







if __name__ == '__main__':

    main()