import argparse
import time
import resource
import os, sys, pandas, numpy, pathlib
import CytoSig
from scipy import stats
from statsmodels.stats.multitest import multipletests

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_file", type=str, default=None, required=False, help="Gene expression file.")
parser.add_argument('-R', "--response_data", type=str, default=None, required=True, help="Precomputed response data frame.")
parser.add_argument('-S', "--signaling_data", type=str, default=None, required=True, help="Precomputed signaling data frame.")
parser.add_argument('-O', "--output_tag", type=str, default=None, required=False, help="Prefix for output files.")

parser.add_argument('-C', "--count_threshold", type=int, default=100, required=False, help="Minimal number of cells needed for regression [100].")
parser.add_argument('-RK', "--response_key", type=str, default='Proliferation', required=False, help="Name of response in the data table [Proliferation].")
parser.add_argument('-SK', "--signaling_key", type=str, default='TGFB1', required=False, help="Name of signaling in the data table [TGFB1].")
parser.add_argument("-QC", "--run_qc", help="whether run QC workflow", action="store_true")
args = parser.parse_args()

err_tol = 1e-8

def read_expression(input_file):
    # read input
    try:
        f = os.path.basename(input_file)
        if f.find('.pickle') >= 0:
            print(f)
            expression = pandas.read_pickle(input_file)
        else:
            expression = pandas.read_csv(input_file, sep='\t', index_col=0)
    except:
        sys.stderr.write('Fail to open input file %s\n' % input_file)
        sys.exit(1)
    
    # gene and sample names must be unique
    assert expression.index.value_counts().max() == 1
    assert expression.columns.value_counts().max() == 1
    
    print('input matrix dimension', expression.shape)
    return expression

def interaction_test(expression, X, y):
    signal = X.loc[:, 'pivot']
    
    failed = []
    merge = []
    
    for gid, arr in expression.iterrows():        
        X.loc[:,'partner'] = arr
        X.loc[:,'interaction'] = arr * signal
    
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

expression = read_expression(args.expression_file)

# check column names
if not result_response.columns.equals(result_signaling.columns):
    print('%s and %s have different column names.\n' % (args.response_data, args.signaling_data))
    sys.exit(1)

if not result_response.columns.equals(expression.columns):
    print('%s and %s have different column names.\n' % (args.response_data, args.expression_file))
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
    report.to_csv(args.output_tag + '_qc.tsv', sep='\t', index_label='ID')
    sys.exit(1)

merge = []
for title, expression_sub in expression_group:
    N = expression_sub.shape[1]
    if N < args.count_threshold: continue
    
    print('process', title, N)
    
    # remove rows all zeros
    flag_nonzero = (expression_sub == 0).mean(axis=1) < 1
    if sum(~flag_nonzero) > 0: expression_sub = expression_sub.loc[flag_nonzero]
    
    y = (result_response.loc[response_key]).loc[expression_sub.columns]

    # regression scaffold
    X = pandas.DataFrame(numpy.zeros((N,4)), columns = ['const', 'pivot', 'partner', 'interaction'], index = expression_sub.columns)
    X.loc[:, 'const'] = 1
    
    X.loc[:,'pivot'] = result_signaling.loc[signaling_key, expression_sub.columns]
    result = interaction_test(expression_sub, X, y)
    result.columns += ('.%s.%s' % (title, signaling_key))
    merge.append(result)

result = pandas.concat(merge, axis=1, join='outer')
result.to_csv(args.output_tag + '_res.tsv', sep='\t', index_label='ID')
