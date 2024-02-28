import argparse
import numpy as np
import pandas as pd
import os
import pandas
import sys
import warnings

import CytoSig
import tqdm

from Util import read_expression

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_path", type=str, required=False, help="Gene expression file.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.Tres_data/sc_cohorts')
parser.add_argument('-R', "--response_path", type=str, required=False, help="Precomputed response data frame.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.qc_data/Prolifertion')
parser.add_argument('-S', "--signaling_path", type=str, required=False, help="Precomputed signaling data frame.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.qc_data/Signaling')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.qc_data/qc_result')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='Breast.GSE156728')
parser.add_argument('-RK', "--response_key", type=str, default='Proliferation', required=False, help="Name of response in the data table [Proliferation].")
parser.add_argument('-CT', "--celltype", type=str, default='CD8', required=False, help="cell type")
args = parser.parse_args()

expression_path = args.expression_path
response_path = args.response_path
signaling_path = args.signaling_path
output_tag = args.output_tag
response_key = args.response_key
celltype = args.celltype
output_file_directory = args.output_file_directory

qc_result = pd.DataFrame(columns=["Dataset", "SampleID", "Cytokine", "t", "p"])
expression_list = sorted(os.listdir(expression_path))
response_list = sorted(os.listdir(response_path))
signaling_list = sorted(os.listdir(signaling_path))

for (expression_file, signaling_file, response_file) in list(zip(expression_list, signaling_list, response_list)):
    # if not expression_file == "Nasopharyngeal.GSE162025.10x.pickle.gz": continue
    expression_filename = os.path.join(expression_path, expression_file)
    response_filename = os.path.join(response_path, response_file)
    signaling_filename = os.path.join(signaling_path, signaling_file)
    # expression_data = read_expression(expression_filename)
    response_data = pandas.read_csv(response_filename, sep='\t', index_col=0)
    signaling_data = pandas.read_csv(signaling_filename, sep='\t', index_col=0)
    # check data
    if response_key not in response_data.index:
        sys.stderr.write('Fail to find %s row in file %s.\n' % (response_key, response_file))
    if not response_data.columns.equals(signaling_data.columns):
        print('%s and %s have different column names.\n' % (response_file, signaling_file))
        continue
    # if not response_data.columns.equals(expression_data.columns):
    #     print('%s and %s have different column names.\n' % (response_file, args.expression_file))
    #     continue

    for signaling_key in signaling_data.index:
    # for signaling_key in ['TGFB1']:
        flag_group = ['.'.join(v.split('.')[:2]) for v in response_data.columns]
        expression_group = response_data.groupby(flag_group, axis=1)
        for title, expression_sub in expression_group:
            if celltype not in title: continue # filter the cell type
            N = expression_sub.shape[1]
            y = (response_data.loc[response_key]).loc[expression_sub.columns]
            y = y.to_frame(name=None)
            X = pandas.DataFrame(np.zeros((N, 1)), columns=['pivot'], index=expression_sub.columns)
            X.loc[:, 'pivot'] = signaling_data.loc[signaling_key, expression_sub.columns]
            try:
                result = CytoSig.ridge_significance_test(X, y, alpha=10000, alternative="two-sided", nrand=1000, verbose=0)
            except ArithmeticError:
                sys.stderr.write('CytoSig regression failed.\n')
                continue
            # result = CytoSig.ridge_significance_test(X, y, alpha=0, alternative="two-sided", nrand=0, verbose=0)
            tvalue = result[2].loc['pivot'].iloc[0]
            pvalue = result[3].loc['pivot'].iloc[0]
            new_row = {'Dataset': output_tag, 'SampleID': title, 'Cytokine': signaling_key, 't': tvalue, 'p': pvalue}
            qc_result.loc[len(qc_result)] = new_row

    print(f"{expression_file} process end.")

qc_result_filename = os.path.join(output_file_directory, f'{celltype}.qc_result.csv')
qc_result.to_csv(qc_result_filename, sep='\t')
print("Process end!")
