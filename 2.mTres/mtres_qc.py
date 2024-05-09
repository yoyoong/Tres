import argparse
import numpy as np
import pandas as pd
import os
import pandas
import sys
import warnings

import CytoSig
from tqdm.autonotebook import tqdm

from Util import read_expression

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_path", type=str, required=False, help="Gene expression file.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.new_gem_data')
parser.add_argument('-R', "--response_path", type=str, required=False, help="Precomputed response data frame.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-3.Neutrophils_response')
parser.add_argument('-S', "--signaling_path", type=str, required=False, help="Precomputed signaling data frame.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2.Signaling')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/4.qc_result')
parser.add_argument('-RK', "--response_key", type=str, default='Neutrophils_signature', required=False,
                    help="Name of response in the data table [Proliferation].")
parser.add_argument('-CT', "--celltype", type=str, required=False, help="cell type", default='Neutrophils')
parser.add_argument('-CTR', "--cohort_celltype_mapping_file", type=str, required=False, help="Celltype mapping rules file, .txt format",
                    default='')
args = parser.parse_args()

expression_path = args.expression_path
response_path = args.response_path
signaling_path = args.signaling_path
response_key = args.response_key
celltype = args.celltype
output_file_directory = args.output_file_directory
cohort_celltype_mapping_file = args.cohort_celltype_mapping_file

if celltype == 'Macrophage':
    celltype_in_column = 'Mono/Macro'
    celltype_in_file = 'Mono_Macro'
else:
    celltype_in_column = celltype
    celltype_in_file = celltype

if cohort_celltype_mapping_file:
    all_mapping_rules_df = pd.read_csv(cohort_celltype_mapping_file, sep='\t')
    celltype_mapping_rules_df = all_mapping_rules_df[all_mapping_rules_df['CELLTYPE'] == celltype]
    celltype_mapping_rules_dict = {key: value for key, value in
                                   zip(all_mapping_rules_df.iloc[:, 1], all_mapping_rules_df.iloc[:, 2])}

qc_result = pd.DataFrame(columns=["Dataset", "SampleID", "Cytokine", "t", "p"])
if os.path.isfile(expression_path):
    dataset_tag = os.path.basename(expression_path).split('.')[0]
    if not os.path.exists(response_path) or not os.path.exists(signaling_path):
        sys.exit()

    # expression_data = read_expression(expression_filename)
    response_data = pandas.read_csv(response_path, sep='\t', index_col=0)
    signaling_data = pandas.read_csv(signaling_path, sep='\t', index_col=0)

    # check data
    if response_key not in response_data.index:
            sys.stderr.write('Fail to find %s row in file %s.\n' % (response_key, response_filename))
    # if not response_data.columns.equals(signaling_data.columns):
    #     print('%s and %s have different column names.\n' % (response_filename, signaling_filename))
    #     continue
    # if not response_data.columns.equals(expression_data.columns):
    #     print('%s and %s have different column names.\n' % (response_file, args.expression_file))
    #     continue

    # get the real cell type
    if cohort_celltype_mapping_file and expression_basename in celltype_mapping_rules_dict.keys():
        dataset_celltype = celltype_mapping_rules_dict[expression_basename]
    else:
        dataset_celltype = celltype_in_column

    for signaling_key in tqdm(signaling_data.index, desc="Processing Cytokine"):
        flag_group = ['.'.join(v.split('.')[:2]) for v in response_data.columns]
        response_group = response_data.groupby(flag_group, axis=1)
        for title, response_sub in response_group:
            sample_celltype = title.split(".")[0].replace(":", "")
            if dataset_celltype != sample_celltype: continue  # filter the cell type

            N = response_sub.shape[1]
            y = (response_data.loc[response_key]).loc[response_sub.columns]
            y = y.to_frame(name=None)
            X = pandas.DataFrame(np.zeros((N, 1)), columns=['pivot'], index=response_sub.columns)
            X.loc[:, 'pivot'] = signaling_data.loc[signaling_key, response_sub.columns]
            try:
                result = CytoSig.ridge_significance_test(X, y, alpha=10000, alternative="two-sided", nrand=1000,
                                                         verbose=0)
            except ArithmeticError:
                sys.stderr.write('CytoSig regression failed.\n')
                continue

            tvalue = result[2].loc['pivot'].iloc[0]
            pvalue = result[3].loc['pivot'].iloc[0]
            new_row = {'Dataset': dataset_tag, 'SampleID': title, 'Cytokine': signaling_key, 't': tvalue,
                       'p': pvalue}
            qc_result.loc[len(qc_result)] = new_row
else:
    expression_list = sorted(os.listdir(expression_path))
    for expression_file in tqdm(expression_list, desc="Processing expression files"):
        expression_basename = os.path.basename(expression_file)
        file_suffix = ".pickle.gz"
        if "tisch" in expression_path:
            file_suffix = ".csv"
        dataset_tag = expression_basename.split(file_suffix)[0]

        expression_filename = os.path.join(expression_path, f'{dataset_tag}{file_suffix}')
        response_filename = os.path.join(response_path, f'{dataset_tag}.csv')
        signaling_filename = os.path.join(signaling_path, f'{dataset_tag}.csv')
        if not os.path.exists(response_filename) or not os.path.exists(signaling_filename):
            continue

        # expression_data = read_expression(expression_filename)
        response_data = pandas.read_csv(response_filename, sep='\t', index_col=0)
        signaling_data = pandas.read_csv(signaling_filename, sep='\t', index_col=0)

        # check data
        if response_key not in response_data.index:
            sys.stderr.write('Fail to find %s row in file %s.\n' % (response_key, response_filename))
        # if not response_data.columns.equals(signaling_data.columns):
        #     print('%s and %s have different column names.\n' % (response_filename, signaling_filename))
        #     continue
        # if not response_data.columns.equals(expression_data.columns):
        #     print('%s and %s have different column names.\n' % (response_file, args.expression_file))
        #     continue

        # get the real cell type
        if cohort_celltype_mapping_file and expression_basename in celltype_mapping_rules_dict.keys():
            dataset_celltype = celltype_mapping_rules_dict[expression_basename]
        else:
            dataset_celltype = celltype_in_column

        for signaling_key in signaling_data.index:
            flag_group = ['.'.join(v.split('.')[:2]) for v in response_data.columns]
            response_group = response_data.groupby(flag_group, axis=1)
            for title, response_sub in response_group:
                sample_celltype = title.split(".")[0].replace(":", "")
                if dataset_celltype != sample_celltype: continue # filter the cell type

                N = response_sub.shape[1]
                y = (response_data.loc[response_key]).loc[response_sub.columns]
                y = y.to_frame(name=None)
                X = pandas.DataFrame(np.zeros((N, 1)), columns=['pivot'], index=response_sub.columns)
                X.loc[:, 'pivot'] = signaling_data.loc[signaling_key, response_sub.columns]
                try:
                    result = CytoSig.ridge_significance_test(X, y, alpha=10000, alternative="two-sided", nrand=1000, verbose=0)
                except ArithmeticError:
                    sys.stderr.write('CytoSig regression failed.\n')
                    continue

                tvalue = result[2].loc['pivot'].iloc[0]
                pvalue = result[3].loc['pivot'].iloc[0]
                new_row = {'Dataset': dataset_tag, 'SampleID': title, 'Cytokine': signaling_key, 't': tvalue, 'p': pvalue}
                qc_result.loc[len(qc_result)] = new_row

qc_result_filename = os.path.join(output_file_directory, f'qc_result.{celltype}.csv')
qc_result.to_csv(qc_result_filename, sep='\t')
print("Process end!")
