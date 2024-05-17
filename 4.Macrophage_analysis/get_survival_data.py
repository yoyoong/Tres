import os
import numpy as np
import pandas as pd
import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.preprocessing import MaxAbsScaler
from sklearn.metrics import mean_absolute_error

celltype = 'Neutrophils'
dataset_list = ["GSE14333", "GSE15654", "GSE17538", "GSE21374", "GSE28221", "GSE65218", "GSE65682", "GSE112927", "GSE33113", "GSE31595"]
# dataset_list = ["GSE31595"]
for dataset in dataset_list:
    output_dir = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/4.Macrophage_analysis/{dataset}'
    gem_df = pd.read_csv(os.path.join(output_dir, f'{dataset}.gem.csv'), index_col=0, header=0)

    survival_df = pd.DataFrame(columns=['sample_name', 'time', 'death'])
    if dataset == 'GSE33113':
        gem_df = gem_df.iloc[:, :90]

    survival_df['sample_name'] = gem_df.columns.values

    if (dataset == 'GSE15654' or dataset == 'GSE21374' or dataset == 'GSE65682' or dataset == 'GSE112927'
            or dataset == 'GSE65218' or dataset == 'GSE17538' or dataset == 'GSE33113' or dataset == 'GSE31595'):
        series_matrix_path = os.path.join(output_dir, f'{dataset}_series_matrix.txt')
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
                elif dataset == 'GSE112927':
                    time_key = 'follow up days'
                    death_key = 'death censored graft loss'
                elif dataset == 'GSE65218':
                    time_key = 'follow_up_time_in_years_until_31.5.2014'
                    death_key = 'status 31.5.2014 (1 died, 0 alive)'
                elif dataset == 'GSE17538':
                    time_key = 'dfs_time'
                    death_key = 'dfs_event (disease free survival; cancer recurrence)'
                elif dataset == 'GSE33113':
                    time_key = 'time to meta or recurrence'
                    death_key = 'meta or recurrence within 3 years'
                elif dataset == 'GSE31595':
                    time_key = 'relapse free survival time'
                    death_key = 'recurrence'

                if line.find(f'!Sample_characteristics_ch1	\"{time_key}') >= 0:
                    item_list = line.replace("\"", "").split('\t')
                    time_list = [np.nan if item.split(": ")[-1].strip() == 'NA' else float(item.split(": ")[-1])
                                 for item in item_list if item.startswith(time_key)]
                if line.find(f'!Sample_characteristics_ch1	\"{death_key}') >= 0:
                    item_list = line.replace("\"", "").split('\t')
                    if dataset == 'GSE17538':
                        death_list = [0 if item.split(": ")[-1].strip() == 'no recurrence' else 1
                                      for item in item_list if item.startswith(death_key)]
                    elif dataset == 'GSE33113' or dataset == 'GSE31595':
                        death_list = [0 if item.split(": ")[-1].strip() == 'no' else 1
                                      for item in item_list if item.startswith(death_key)]
                    else:
                        death_list = [np.nan if item.split(": ")[-1].strip() == 'NA' else int(item.split(": ")[-1])
                                      for item in item_list if item.startswith(death_key)]
                if line.find('!series_matrix_table_begin') >= 0:
                    break
        survival_df['time'] = np.array(time_list)
        survival_df['death'] = np.array(death_list)

    elif dataset == 'GSE14333':
        series_matrix_path = os.path.join(output_dir, f'{dataset}_series_matrix.txt')
        time_list = []
        death_list = []
        with open(series_matrix_path, mode='r') as file:
            lines = file.readlines()
            for line in lines:
                if line.find(f'!Sample_characteristics_ch1') >= 0:
                    item1_list = line.replace("\"", "").split('\t')
                    for item1 in item1_list:
                        item2_list = item1.split('; ')
                        for item2 in item2_list:
                            if item2.startswith("DFS_Time"):
                                if item2.split(": ")[-1].strip() == 'NA':
                                    time_list.append(np.nan)
                                else:
                                    time_list.append(float(item2.split(": ")[-1]))
                            if item2.startswith("DFS_Cens"):
                                if item2.split(": ")[-1].strip() == 'NA':
                                    death_list.append(np.nan)
                                else:
                                    death_list.append(int(item2.split(": ")[-1]))
            if line.find('!series_matrix_table_begin') >= 0:
                break
        survival_df['time'] = np.array(time_list)
        survival_df['death'] = np.array(death_list)

    elif dataset == 'GSE28221':
        annotation_path = os.path.join(output_dir, 'GSE28221_readme.txt')
        annotation_df = pd.read_csv(annotation_path, sep='\t', index_col=2, header=0)
        ffff = annotation_df.loc[survival_df['sample_name'], 'characteristics: time to survival']

        survival_df['time'] = np.array(annotation_df.loc[survival_df['sample_name'], 'characteristics: time to survival'])
        survival_df['death'] = np.array(annotation_df.loc[survival_df['sample_name'], 'characteristics: survival status, 0 = censored, 1 = death'])


    bulk_response_path = os.path.join(output_dir, f'{dataset}.bulk_response.csv')
    bulk_response_data = pd.read_csv(bulk_response_path, delimiter='\t', index_col=0)
    for signature_name in ['SMART_C13', 'SMART_C3', 'SMART_C14']:
    # for signature_name in ['Polarization']:
        signature_key = signature_name.split('_')[-1] + '_response'
        survival_df[signature_key] = survival_df['sample_name'].apply(
            lambda x: bulk_response_data.loc[signature_name][x])

        # cutoff = np.median(survival_df[signature_key])
        # def get_response_group(x):
        #     if x[signature_key] >= cutoff:
        #         return "high"
        #     else:
        #         return "low"
        def get_response_group(x):
            if x[signature_key] >= 0:
                return "t>0"
            else:
                return "t<0"
        survival_df[f'{signature_key}_group'] = survival_df.apply(lambda x: get_response_group(x), axis=1)

    survival_df.set_index('sample_name', inplace=True)
    for flag in [1, 2, 3]:
        tres_signature_tag = f'Tres_signature_{flag}.positive'
        if celltype == 'Macrophage':
            tres_signature_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-2.Macrophage_Interaction/{tres_signature_tag}.csv'
        elif celltype == 'Neutrophils':
            # tres_signature_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/5-3.Neutrophils_Interaction/{tres_signature_tag}.csv'
            tres_signature_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.neutrophil_data/{tres_signature_tag}.csv'
        tres_signature_df = pd.read_csv(tres_signature_path, index_col=0, header=0)

        # calculate the correlation of expression and tres signature
        survival_df[f'correlation_{flag}'] = None
        # survival_df[f'mse'] = None
        common_gene_list = gem_df.index.intersection(tres_signature_df.index)
        tres_signature = tres_signature_df.loc[common_gene_list]['Tres']
        gem_filtered = gem_df.loc[common_gene_list]
        for sample in gem_df.columns.values:
            gem_data = gem_filtered[sample]
            correlation, _ = pearsonr(np.array(gem_data), np.array(tres_signature))
            survival_df.loc[sample, f'correlation_{flag}'] = correlation

            # scaler = MaxAbsScaler()
            # gem_data_normalized = scaler.fit_transform(gem_data.values.reshape(-1, 1)).ravel()
            # tres_signature_normalized = scaler.fit_transform(tres_signature.values.reshape(-1, 1)).ravel()
            # mse = mean_absolute_error(gem_data, tres_signature)
            # survival_df.loc[sample, f'mse'] = mse

        def get_corr_group1(x):
            if x[f'correlation_{flag}'] >= 0:
                return "Corr>0"
            else:
                return "Corr<0"
        survival_df[f'corr_group_{flag}'] = survival_df.apply(lambda x: get_corr_group1(x), axis=1)

        # cutoff = np.median(survival_df[f'correlation_{flag}'])
        # def get_corr_group2(x):
        #     if x[f'correlation_{flag}'] >= cutoff:
        #         return "Corr high"
        #     else:
        #         return "Corr low"
        # survival_df[f'corr_group2_{flag}'] = survival_df.apply(lambda x: get_corr_group2(x), axis=1)

        survival_df.to_csv(os.path.join(output_dir, 'survival', f'{dataset}.{celltype}.{tres_signature_tag}.csv'))
    print(f"{dataset} process end")
