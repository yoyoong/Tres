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
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")

celltype = 'NK'
expression_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.new_gem_data'
signaling_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/2.Signaling'
output_file_directory = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/4.qc_result'
count_threshold = 100
cohort_celltype_mapping_file = ''

if celltype == 'CD8T':
    response_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-1.Proliferation'
    response_key = 'Proliferation'
elif celltype == 'Macrophage':
    response_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-2.Polarization'
    response_key = 'Polarization'
elif celltype == 'Neutrophils':
    response_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-3.Neutrophils_response'
    response_key = 'Neut_IFN-15'
elif celltype == 'NK':
    response_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-4.NK_response'
    response_key = 'NK_signature'
elif celltype == 'NFkB' or celltype == 'Hif1a':
    signaling_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-4.NK_response'
    response_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-4-2.{celltype}_response'
    response_key = f'{celltype}_signature'

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

qc_result = pd.DataFrame(columns=["Dataset", "SampleID", "Cytokine", "correlation", "p"])
if os.path.isfile(expression_path):
    dataset_tag = os.path.basename(expression_path).split('.')[0]
    if not os.path.exists(response_path) or not os.path.exists(signaling_path):
        sys.exit()

    # expression_data = read_expression(expression_filename)
    response_data = pandas.read_csv(response_path, sep='\t', index_col=0)
    signaling_data = pandas.read_csv(signaling_path, sep='\t', index_col=0)
    if len(response_data) < 1 or len(signaling_data) < 1:
        print(f'{dataset_tag} data is null')
        sys.exit()

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
            response = np.array((response_data.loc[response_key]).loc[response_sub.columns])
            signaling = np.array(signaling_data.loc[signaling_key, response_sub.columns])
            if len(response) < count_threshold or len(signaling) < count_threshold:
                # print("Length of response or signaling < 2.")
                continue

            correlation, pvalue = pearsonr(response, signaling)
            new_row = {'Dataset': dataset_tag, 'SampleID': title, 'Cytokine': signaling_key, 'correlation': correlation, 'p': pvalue}
            qc_result.loc[len(qc_result)] = new_row
else:
    expression_list = sorted(os.listdir(expression_path))
    if celltype == 'Neutrophils':
        expression_list.append('Gao2024')
    for expression_file in tqdm(expression_list, desc="Processing expression files"):
        expression_basename = os.path.basename(expression_file)
        file_suffix = ".pickle.gz"
        if "tisch" in expression_path:
            file_suffix = ".csv"
        dataset_tag = expression_basename.split(file_suffix)[0]

        # expression_filename = os.path.join(expression_path, f'{dataset_tag}{file_suffix}')
        response_filename = os.path.join(response_path, f'{dataset_tag}.csv')
        signaling_filename = os.path.join(signaling_path, f'{dataset_tag}.csv')
        if not os.path.exists(response_filename) or not os.path.exists(signaling_filename):
            continue

        # expression_data = read_expression(expression_filename)
        response_data = pandas.read_csv(response_filename, sep='\t', index_col=0)
        signaling_data = pandas.read_csv(signaling_filename, sep='\t', index_col=0)
        if len(response_data) < 1 or len(signaling_data) < 1:
            print(f'{dataset_tag} data is null')
            continue

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
            if signaling_key == 'study bias': continue
            flag_group = ['.'.join(v.split('.')[:2]) for v in response_data.columns]
            response_group = response_data.groupby(flag_group, axis=1)
            for title, response_sub in response_group:
                sample_celltype = title.split(".")[0].replace(":", "")
                # if dataset_celltype != sample_celltype: continue # filter the cell type
                response = np.array((response_data.loc[response_key]).loc[response_sub.columns])
                signaling = np.array(signaling_data.loc[signaling_key, response_sub.columns])
                if len(response) < count_threshold or len(signaling) < count_threshold:
                    # print("Length of response or signaling < 2.")
                    continue

                response = np.nan_to_num(response, nan=0.0, copy=True)
                signaling = np.nan_to_num(signaling, nan=0.0, copy=True)
                correlation, pvalue = pearsonr(response, signaling)
                new_row = {'Dataset': dataset_tag, 'SampleID': title, 'Cytokine': signaling_key, 'correlation': correlation, 'p': pvalue}
                qc_result.loc[len(qc_result)] = new_row

qc_result_filename = os.path.join(output_file_directory, f'qc_result_correlation.{celltype}.csv')
qc_result.to_csv(qc_result_filename, sep='\t')
print("Process end!")
