import os
import numpy as np
import pandas as pd
import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

dataset_list = ["GSE28221"]
for dataset in dataset_list:
    output_dir = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/4.Macrophage_analysis/{dataset}'
    gem_df = pd.read_csv(os.path.join(output_dir, f'{dataset}.gem.csv'), index_col=0, header=0)

    survival_df = pd.DataFrame(columns=['sample_name', 'time', 'death'])
    survival_df['sample_name'] = gem_df.columns.values

    if dataset == 'GSE15654' or dataset == 'GSE21374' or dataset == 'GSE65682':
        if dataset == 'GSE15654':
            series_matrix_path = os.path.join(output_dir, 'GSE15654_series_matrix.txt')
        elif dataset == 'GSE21374':
            series_matrix_path = os.path.join(output_dir, 'GSE21374_series_matrix.txt')
        elif dataset == 'GSE65682':
            series_matrix_path = os.path.join(output_dir, 'GSE65682_series_matrix.txt')

        with open(series_matrix_path, mode='r') as file:
            lines = file.readlines()
            for line in lines:
                if dataset == 'GSE15654':
                    time_key = 'days to death'
                    death_key = 'death'
                elif dataset == 'GSE21374':
                    time_key = 'time from biopsy to failure/censoring (days)'
                    death_key = 'failed=1/non failed=0'
                elif dataset == 'GSE65682':
                    time_key = 'time_to_event_28days'
                    death_key = 'mortality_event_28days'

                if line.find(f'!Sample_characteristics_ch1	\"{time_key}') >= 0:
                    item_list = line.replace("\"", "").split('\t')
                    time_list = [np.nan if item.split(": ")[-1].strip() == 'NA' else int(item.split(": ")[-1])
                                 for item in item_list if item.startswith(time_key)]
                if line.find(f'!Sample_characteristics_ch1	\"{death_key}') >= 0:
                    item_list = line.replace("\"", "").split('\t')
                    death_list = [np.nan if item.split(": ")[-1].strip() == 'NA' else int(item.split(": ")[-1])
                                  for item in item_list if item.startswith(death_key)]
                if line.find('!series_matrix_table_begin') >= 0:
                    break
        survival_df['time'] = np.array(time_list)
        survival_df['death'] = np.array(death_list)
    else:
        annotation_path = os.path.join(output_dir, 'GSE28221_readme.txt')
        annotation_df = pd.read_csv(annotation_path, sep='\t', index_col=2, header=0)
        ffff = annotation_df.loc[survival_df['sample_name'], 'characteristics: time to survival']

        survival_df['time'] = np.array(annotation_df.loc[survival_df['sample_name'], 'characteristics: time to survival'])
        survival_df['death'] = np.array(annotation_df.loc[survival_df['sample_name'], 'characteristics: survival status, 0 = censored, 1 = death'])

    bulk_response_path = os.path.join(output_dir, f'{dataset}.bulk_response.csv')
    bulk_response_data = pd.read_csv(bulk_response_path, delimiter='\t', index_col=0)
    for signature_name in ['SMART_C13', 'SMART_C3']:
        signature_key = signature_name.split('_')[-1] + '_response'
        survival_df[signature_key] = survival_df['sample_name'].apply(
            lambda x: bulk_response_data.loc[signature_name][x])

        def get_response_group(x):
            if x[signature_key] >= 0:
                return "t>0"
            else:
                return "t<0"
        survival_df[f'{signature_key}_group'] = survival_df.apply(lambda x: get_response_group(x), axis=1)

    survival_df.set_index('sample_name', inplace=True)
    for flag in [1, 2, 3]:
        tres_signature_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction/Tres_signature_{flag}.positive.csv'
        tres_signature_df = pd.read_csv(tres_signature_path, index_col=0, header=0)

        # calculate the correlation of expression and tres signature
        survival_df[f'correlation_{flag}'] = None
        common_gene_list = gem_df.index.intersection(tres_signature_df.index)
        tres_signature = tres_signature_df.loc[common_gene_list]['Tres']
        gem_filtered = gem_df.loc[common_gene_list]
        for sample in gem_df.columns.values:
            gem_data = gem_filtered[sample]
            correlation, _ = pearsonr(np.array(gem_data), np.array(tres_signature))
            survival_df.loc[sample, f'correlation_{flag}'] = correlation

        def get_corr_group(x):
            if x[f'correlation_{flag}'] >= 0:
                return "Corr>0"
            else:
                return "Corr<0"

        survival_df[f'corr_group_{flag}'] = survival_df.apply(lambda x: get_corr_group(x), axis=1)

    survival_df.to_csv(os.path.join(output_dir, f'{dataset}.survival.csv'))
    print(f"{dataset} process end")
