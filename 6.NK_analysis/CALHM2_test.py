import argparse
import numpy as np
import pandas as pd
import os
from scipy.stats import pearsonr
import warnings

import CytoSig
from tqdm.autonotebook import tqdm
warnings.filterwarnings("ignore")

signature_list = ["NK", "NK_act"]
for signature in signature_list:
    if signature == "NK":
        response_dir = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-4.NK_response'
        signature_name = 'NK_signature'
    elif signature == "NK_act":
        response_dir = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-4-0.NK_act_response'
        signature_name = 'NK_act_signature'

    expression_dir = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.new_gem_data'
    count_threshold = 100

    test_result = pd.DataFrame(columns=["Dataset", "SampleID", "correlation", "p"])
    response_list = sorted(os.listdir(response_dir))
    for response_file in tqdm(response_list, desc="Processing response files"):
        dataset_tag = response_file.split('.')[0]
        response_filename = os.path.join(response_dir, response_file)
        expression_filename = os.path.join(expression_dir, dataset_tag, f'{dataset_tag}.NK.csv')

        response_data = pd.read_csv(response_filename, sep='\t', index_col=0)
        expression_data = pd.read_csv(expression_filename, index_col=0)
        if 'CALHM2' not in expression_data.index.values.tolist(): continue

        flag_group = ['.'.join(v.split('.')[:2]) for v in response_data.columns]
        response_group = response_data.groupby(flag_group, axis=1)
        expression_group = expression_data.groupby(flag_group, axis=1)
        for sample, response_sub in response_group:
            CALHM2_expression_value = expression_group.get_group(sample).loc['CALHM2', :]
            signature_value = response_sub.loc[signature_name, :]
            if len(CALHM2_expression_value) < count_threshold or len(signature_value) < count_threshold:
                continue
            correlation, pvalue = pearsonr(CALHM2_expression_value, signature_value)
            new_row = {'Dataset': dataset_tag, 'SampleID': sample, 'correlation': correlation, 'p': pvalue}
            test_result.loc[len(test_result)] = new_row

    output_dir = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/6.NK_analysis'
    test_result.to_csv(os.path.join(output_dir, f'CALHM2_test.{signature}.csv'))
    print(f'{signature} process end')