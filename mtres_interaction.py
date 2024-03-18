import argparse
import numpy
import os
import pandas as pd
import sys
import warnings
from tqdm.autonotebook import tqdm

import CytoSig
from statsmodels.stats.multitest import multipletests

from Util import read_expression

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_path", type=str, required=False, help="Gene expression file.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/1.gem_data/UVM_GSE139829')
parser.add_argument('-R', "--response_data", type=str, required=False, help="Precomputed response data frame.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-1.Prolifertion/UVM_GSE139829.csv')
parser.add_argument('-S', "--signaling_data", type=str, required=False, help="Precomputed signaling data frame.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-2.Signaling/UVM_GSE139829.csv')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/4.Interaction')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='HNSC_GSE139324')
parser.add_argument('-C', "--count_threshold", type=int, default=100, required=False, help="Minimal number of cells needed for regression [100].")
parser.add_argument('-RK', "--response_key", type=str, default='Proliferation', required=False, help="Name of response in the data table [Proliferation].")
parser.add_argument('-SK', "--signaling_key", type=str, default='', required=False, help="Name of signaling in the data") # if null, calculate all cytokine
parser.add_argument('-CT', "--celltype", type=str, default='CD8T', required=False, help="cell type")
args = parser.parse_args()

expression_path = args.expression_path
response_data = args.response_data
signaling_data = args.signaling_data
output_file_directory = args.output_file_directory
output_tag = args.output_tag
count_threshold = args.count_threshold
response_key = args.response_key
signaling_key = args.signaling_key
celltype = args.celltype

err_tol = 1e-8
def interaction_test(expression, X, y):
    signal = X.loc[:, 'pivot']
    failed = []
    merge = []

    # 对每个基因分别求Tres分数
    for gid, arr in expression.iterrows():
        X.loc[:, 'partner'] = arr  # b * G
        X.loc[:, 'interaction'] = arr * signal  # c * G * suppression

        # other covariates have no sufficient variation
        if arr.std() < err_tol or X.loc[:, 'interaction'].std() < err_tol: continue

        try:
            y = pd.DataFrame(y)
            result = CytoSig.ridge_significance_test(X, y, alpha=0, alternative="two-sided", nrand=0,
                                                     flag_normalize=False, verbose=False)

        except ArithmeticError as e:
            # print(f"{gid} ArithmeticError:", str(e))
            failed.append(gid)
            continue

        tvalue = result[2].loc['interaction'].iloc[0]
        pvalue = result[3].loc['interaction'].iloc[0]
        merge.append(pd.Series([tvalue, pvalue], index=['t', 'p'], name=gid))

    if len(merge) > 0:
        result = pd.concat(merge, axis=1, join='inner').transpose()
        result['q'] = multipletests(result['p'], method='fdr_bh')[1]
        return result
    else:
        return None

result_response = pd.read_csv(response_data, sep='\t', index_col=0)
if response_key not in result_response.index:
    sys.stderr.write('Fail to find %s row in file %s.\n' % (response_key, response_data))

result_signaling = pd.read_csv(signaling_data, sep='\t', index_col=0)
signaling_list = []
if not signaling_key:
    signaling_list = list(result_signaling.index)
else:
    if signaling_key not in result_signaling.index:
	    sys.stderr.write('Fail to find %s row in file %s.\n' % (signaling_key, signaling_data))
    signaling_list.append(signaling_key)
if not result_response.columns.equals(result_signaling.columns):
    print('%s and %s have different column names.\n' % (response_data, signaling_data))
    sys.exit(1)

if not os.path.isdir(expression_path):
    expression = read_expression(expression_path)
    if not result_response.columns.equals(expression.columns):
        print('%s and %s have different column names.\n' % (response_data, expression_path))
        sys.exit(1)

    celltype_list = [v.split('.')[0] for v in expression.columns]
    if celltype not in celltype_list:
        print(f"This dataset has not {celltype} celltype.")
        sys.exit(1)

    flag_group_filtered = [v for v in expression.columns if v.split('.')[0] == celltype]
    expression_filtered = expression[flag_group_filtered]
else:
    expression_list = os.listdir(expression_path)
    celltype_list = [os.path.basename(v).split('.')[1] for v in expression_list]
    if celltype not in celltype_list:
        print(f"This dataset has not {celltype} celltype.")
        sys.exit(1)

    tag = expression_path.split('/')[-1]
    expression_filtered = pd.read_csv(os.path.join(expression_path, f'{tag}.{celltype}.csv'), index_col=0)

flag_group = ['.'.join(v.split('.')[:2]) for v in expression_filtered.columns]
expression_group = expression_filtered.groupby(flag_group, axis=1)
result_all = pd.DataFrame(columns=['GeneID', 'Cytokine', 'SampleID', 't', 'p', 'q'])
for sample, expression_sub in tqdm(expression_group, desc="Sample processing"):
    # filter the cell type
    sample_celltype = sample.split(".")[0]
    if sample_celltype != celltype: continue

    # filter the cell count
    N = expression_sub.shape[1]
    if N < count_threshold:
        print(f"Sample:{sample} cell number is {N}, less than {count_threshold}.")
        continue

    # remove rows all zeros
    flag_nonzero = (expression_sub == 0).mean(axis=1) < 1
    if sum(~flag_nonzero) > 0:
        expression_sub = expression_sub.loc[flag_nonzero]
    y = (result_response.loc[response_key]).loc[expression_sub.columns]

    for signaling_name in signaling_list:
        # regression scaffold
        X = pd.DataFrame(numpy.zeros((N, 4)), columns = ['const', 'pivot', 'partner', 'interaction'], index = expression_sub.columns)
        X.loc[:, 'const'] = 1 # d * 1
        X.loc[:, 'pivot'] = result_signaling.loc[signaling_name, expression_sub.columns] # a * suppression
        result_sub = interaction_test(expression_sub, X, y)
        if result_sub is not None:
            result_sub['GeneID'] = result_sub.index
            result_sub['Cytokine'] = signaling_name
            result_sub['SampleID'] = sample
            result_all = pd.concat([result_all, result_sub], ignore_index=True)
            print(f"Signal: {signaling_name}, Sample: {sample} calculate end.")
        else:
            print(f"Sample: {sample}, Signal: {signaling_name} no result.")

if result_all.shape[0] > 0:
    result_all.set_index(['Cytokine', 'GeneID', 'SampleID'], inplace=True)
    result_all.sort_values(['Cytokine', 'GeneID', 'SampleID'], inplace=True)
    result_filename = os.path.join(output_file_directory, f'{output_tag}.{celltype}.csv')
    result_all.to_csv(result_filename)
print("Process end!")