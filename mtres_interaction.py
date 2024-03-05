import argparse
import numpy
import os
import pandas
import sys
import warnings

import CytoSig
from statsmodels.stats.multitest import multipletests

from Util import read_expression

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_path", type=str, required=False, help="Gene expression file.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.Tres_data/sc_cohorts/Breast.GSE156728.10x.pickle.gz')
parser.add_argument('-R', "--response_data", type=str, required=False, help="Precomputed response data frame.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.qc_data/Prolifertion/Breast.GSE156728.10x.Prolifertion.csv')
parser.add_argument('-S', "--signaling_data", type=str, required=False, help="Precomputed signaling data frame.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.qc_data/Signaling/Breast.GSE156728.10x.Signaling.csv')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.qc_data/qc_result')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='Breast.GSE156728.10x')

parser.add_argument('-C', "--count_threshold", type=int, default=100, required=False, help="Minimal number of cells needed for regression [100].")
parser.add_argument('-RK', "--response_key", type=str, default='Proliferation', required=False, help="Name of response in the data table [Proliferation].")
parser.add_argument('-SK', "--signaling_key", type=str, default='TGFB1', required=False,
                    help="Name of signaling in the data table [TGFB1].") # if null, calculate all cytokine
parser.add_argument("-QC", "--run_qc", help="whether run QC workflow", action="store_true")
args = parser.parse_args()

err_tol = 1e-8

def interaction_test(expression, X, y):
    signal = X.loc[:, 'pivot']
    
    failed = []
    merge = []

    # 对每个基因分别求Tres分数
    for gid, arr in expression.iterrows():        
        X.loc[:,'partner'] = arr # b * G
        X.loc[:,'interaction'] = arr * signal # c * G * suppression
    
        # other covariates have no sufficient variation
        if arr.std() < err_tol or X.loc[:,'interaction'].std() < err_tol: continue
        
        try:
            y = pandas.DataFrame(y)
            result = CytoSig.ridge_significance_test(X, y, alpha=0, alternative="two-sided", nrand=0, flag_normalize=False, verbose=False)
        
        except ArithmeticError:
            failed.append(gid)
            continue
        
        tvalue = result[2].loc['interaction'].iloc[0]
        pvalue = result[3].loc['interaction'].iloc[0]
        
        merge.append(pandas.Series([tvalue, pvalue], index=['t', 'p'], name=gid))    
    
    result = pandas.concat(merge, axis=1, join='inner').transpose()
    result['q'] = multipletests(result['p'], method='fdr_bh')[1]
    
    return result


# read response data frame
response_key = args.response_key
signaling_key = args.signaling_key

result_response = pandas.read_csv(args.response_data, sep='\t', index_col=0)
if response_key not in result_response.index:
	sys.stderr.write('Fail to find %s row in file %s.\n' % (response_key, args.response_data))

# read signaling data frame
result_signaling = pandas.read_csv(args.signaling_data, sep='\t', index_col=0)
if signaling_key not in result_signaling.index:
	sys.stderr.write('Fail to find %s row in file %s.\n' % (signaling_key, args.signaling_data))

expression = read_expression(args.expression_path)

# check column names
if not result_response.columns.equals(result_signaling.columns):
    print('%s and %s have different column names.\n' % (args.response_data, args.signaling_data))
    sys.exit(1)

if not result_response.columns.equals(expression.columns):
    print('%s and %s have different column names.\n' % (args.response_data, args.expression_path))
    sys.exit(1)


flag_group = ['.'.join(v.split('.')[:2]) for v in expression.columns]
expression_group = expression.groupby(flag_group, axis=1)

report = []
if args.run_qc:
    for title, expression_sub in expression_group:
        N = expression_sub.shape[1]
        print('process', title, N)
        y = (result_response.loc[response_key]).loc[expression_sub.columns]
        y = y.to_frame(name=None)
        X = pandas.DataFrame(numpy.zeros((N,1)), columns = ['pivot'], index = expression_sub.columns)
        X.loc[:,'pivot'] = result_signaling.loc[signaling_key, expression_sub.columns]
        result = CytoSig.ridge_significance_test(X, y, alpha=10000, alternative="two-sided", nrand=1000, verbose=1)
        tvalue = result[2].loc['pivot'].iloc[0]
        pvalue = result[3].loc['pivot'].iloc[0]
        report.append(pandas.Series([tvalue, pvalue], index=['t', 'p'], name=title))

    report = pandas.concat(report, axis=1, join='inner').transpose()
    report_filename = os.path.join(args.output_file_directory, f'{args.output_tag}.qc.csv')
    report.to_csv(report_filename, sep='\t', index_label='ID')
    sys.exit(1)

merge = []
for title, expression_sub in expression_group:
    N = expression_sub.shape[1]
    if N < args.count_threshold: continue
    
    print('process', title, N)
    
    # remove rows all zeros
    flag_nonzero = (expression_sub == 0).mean(axis=1) < 1
    if sum(~flag_nonzero) > 0:
        expression_sub = expression_sub.loc[flag_nonzero]
    
    y = (result_response.loc[response_key]).loc[expression_sub.columns]

    # regression scaffold
    X = pandas.DataFrame(numpy.zeros((N,4)), columns = ['const', 'pivot', 'partner', 'interaction'], index = expression_sub.columns)
    X.loc[:, 'const'] = 1 # d * 1
    X.loc[:,'pivot'] = result_signaling.loc[signaling_key, expression_sub.columns] # a * suppression
    result = interaction_test(expression_sub, X, y)
    result.columns += ('.%s.%s' % (title, signaling_key))
    merge.append(result)

result = pandas.concat(merge, axis=1, join='outer')
result_filename = os.path.join(args.output_file_directory, f'{args.output_tag}.tres.csv')
result.to_csv(result_filename, sep='\t', index_label='ID')
print("Process end!")