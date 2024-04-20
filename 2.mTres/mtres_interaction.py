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
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.new_gem_data/ALL_GSE153697')
parser.add_argument('-R', "--response_path", type=str, required=False, help="Precomputed response data frame.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-2.Polarization/ALL_GSE153697.csv')
parser.add_argument('-S', "--signaling_path", type=str, required=False, help="Precomputed signaling data frame.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2.Signaling/ALL_GSE153697.csv')
parser.add_argument('-CTR', "--celltype_mapping_rules_file", type=str, required=False, help="for paper data use",
                    default='')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction/dataset_interaction')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='ALL_GSE153697')
parser.add_argument('-C', "--count_threshold", type=int, default=100, required=False, help="Minimal number of cells needed for regression [100].")
parser.add_argument('-RK', "--response_key", type=str, default='Polarization', required=False, help="Name of response in the data table [Proliferation].")
parser.add_argument('-SK', "--signaling_key", type=str, default='', required=False, help="Name of signaling in the data") # if null, calculate all cytokine
parser.add_argument('-CT', "--celltype", type=str, default='Macrophage', required=False, help="cell type")
args = parser.parse_args()

expression_path = args.expression_path
response_path = args.response_path
signaling_path = args.signaling_path
celltype_mapping_rules_file = args.celltype_mapping_rules_file
output_file_directory = args.output_file_directory
output_tag = args.output_tag
count_threshold = args.count_threshold
response_key = args.response_key
signaling_key = args.signaling_key
celltype = args.celltype

if celltype == 'Macrophage':
    celltype_in_column = 'Mono/Macro'
    celltype_in_file = 'Mono_Macro'
else:
    celltype_in_column = celltype
    celltype_in_file = celltype

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

response_data = pd.read_csv(response_path, sep='\t', index_col=0)

signaling_data = pd.read_csv(signaling_path, sep='\t', index_col=0)
signaling_list = []
if not signaling_key:
    signaling_list = list(signaling_data.index)
else:
    if signaling_key not in signaling_data.index:
	    sys.stderr.write('Fail to find %s row in file %s.\n' % (signaling_key, signaling_data))
    signaling_list.append(signaling_key)
# if not response_data.columns.equals(signaling_data.columns):
#     print('%s and %s have different column names.\n' % (response_data, signaling_data))
#     sys.exit(1)

if not os.path.isdir(expression_path):
    expression = read_expression(expression_path)
    if not response_data.columns.equals(expression.columns):
        print('%s and %s have different column names.\n' % (response_data, expression_path))
        sys.exit(1)

    if celltype_mapping_rules_file:
        celltype_mapping_rules_df = pd.read_csv(celltype_mapping_rules_file, sep='\t', index_col=1)
        celltype_real = celltype_mapping_rules_df.loc[os.path.basename(expression_path)]['FILE_CELLTYPE']
        if isinstance(celltype_real, str):
            flag_group_filtered = [v for v in expression.columns if v.split('.')[0] == celltype_real]
        else:
            flag_group_filtered = [v for v in expression.columns if v.split('.')[0] in list(celltype_real)]
    else:
        celltype_list = [v.split('.')[0] for v in expression.columns]
        if celltype_in_column not in celltype_list:
            print(f"This dataset has not {celltype} celltype.")
            sys.exit(1)
        flag_group_filtered = [v for v in expression.columns if v.split('.')[0] == celltype_in_column]
    expression_filtered = expression[flag_group_filtered]
else:
    expression_list = os.listdir(expression_path)
    celltype_list = [os.path.basename(v).split('.')[1] for v in expression_list]
    if celltype_in_file not in celltype_list:
        print(f"This dataset has not {celltype} celltype.")
        sys.exit(1)

    tag = expression_path.split('/')[-1]
    expression_filtered = pd.read_csv(os.path.join(expression_path, f'{tag}.{celltype_in_file}.csv'), index_col=0)

flag_group = ['.'.join(v.split('.')[:2]) for v in expression_filtered.columns]
expression_group = expression_filtered.groupby(flag_group, axis=1)
result_all = []
for sample, expression_sub in tqdm(expression_group, desc="Processing sample"):
    # filter the cell type
    sample_celltype = sample.split(".")[0]
    if celltype_mapping_rules_file:
        if isinstance(celltype_real, str):
            if sample_celltype != celltype_real: continue
        else:
            if sample_celltype not in list(celltype_real): continue
    else:
        if sample_celltype != celltype_in_column: continue

    # filter the cell count
    N = expression_sub.shape[1]
    if N < count_threshold:
        print(f"Sample:{sample} cell number is {N}, less than {count_threshold}.")
        continue

    # remove rows all zeros
    flag_nonzero = (expression_sub == 0).mean(axis=1) < 1
    if sum(~flag_nonzero) > 0:
        expression_sub = expression_sub.loc[flag_nonzero]
    y = (response_data.loc[response_key]).loc[expression_sub.columns]

    for signaling_name in tqdm(signaling_list, desc="Processing signaling"):
        # regression scaffold
        X = pd.DataFrame(numpy.zeros((N, 4)), columns = ['const', 'pivot', 'partner', 'interaction'], index = expression_sub.columns)
        X.loc[:, 'const'] = 1 # d * 1
        X.loc[:, 'pivot'] = signaling_data.loc[signaling_name, expression_sub.columns] # a * suppression
        result_sub = interaction_test(expression_sub, X, y)
        result_sub.columns += ('.%s.%s' % (sample, signaling_name))
        result_all.append(result_sub)

if len(result_all) > 0:
    result = pd.concat(result_all, axis=1, join='outer')
    result_filename = os.path.join(output_file_directory, f'{output_tag}.{celltype}.csv')
    result.to_csv(result_filename)
print("Process end!")