import numpy as np
import pandas as pd
import os
from tqdm.autonotebook import tqdm

dataset_list = ['TCGA-COAD', 'TCGA-LUAD', 'TCGA-PAAD', 'TCGA-STAD']
for dataset in dataset_list:
    rawdata_dir = f'/sibcb1/bioinformatics/hongyuyang/dataset/Tres/3-5.B_analysis/{dataset}/rawdata'
    metadata_path = f'/sibcb1/bioinformatics/hongyuyang/dataset/Tres/3-5.B_analysis/{dataset}/metadata.json'
    clinical_path = f'/sibcb1/bioinformatics/hongyuyang/dataset/Tres/3-5.B_analysis/{dataset}/clinical.tsv'
    outputdir = f'/sibcb1/bioinformatics/hongyuyang/dataset/Tres/3-5.B_analysis/{dataset}'
    
    metadata = pd.read_json(metadata_path)
    rawdata_files = os.listdir(rawdata_dir)
    
    # get gene expression matrix
    gem_df = pd.DataFrame()
    for rawdata_filename in tqdm(rawdata_files, desc=f'Process {dataset}'):
        try:
            case_id = metadata[metadata['file_name'] == rawdata_filename].values[0, 2][0]['case_id']
        except:
            print(f"Get {case_id} fail.")
            continue
        rawdata_path = os.path.join(rawdata_dir, rawdata_filename)
        rawdata = pd.read_csv(rawdata_path, index_col=0, header=0, delimiter='\t', comment='#', skiprows=[2, 3, 4, 5])
        rawdata.set_index('gene_name', inplace=True)
        rawdata_filter = rawdata.groupby(rawdata.index).max()
        count_value = pd.DataFrame(rawdata_filter['unstranded'])
        count_value.columns = [case_id]
        gem_df = pd.concat([gem_df, count_value], axis=1, join='outer')
    
    # normalize
    col_sums = gem_df.sum() # 计算每列的总和
    non_zero_cols = col_sums[col_sums != 0].index # 找到非零列
    for col in non_zero_cols: # 对非零列进行归一化
        gem_df[col] = gem_df[col] / gem_df[col].sum() * 10000
    # log
    gem_df = np.log2(gem_df + 1)
    gem_df.to_csv(os.path.join(outputdir, 'gem.csv'))
    print(f"{dataset} process end.")