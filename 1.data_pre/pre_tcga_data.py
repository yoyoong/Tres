import numpy as np
import pandas as pd
import os
from tqdm.autonotebook import tqdm

rawdata_dir = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/6.NK_analysis/TCGA-SKCM/rawdata'
metadata_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/6.NK_analysis/TCGA-SKCM/metadata.json'
clinical_path = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/6.NK_analysis/TCGA-SKCM/clinical.tsv'
outputdir = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/6.NK_analysis/TCGA-SKCM'

metadata = pd.read_json(metadata_path)
rawdata_files = os.listdir(rawdata_dir)

# get gene expression matrix
gem_df = pd.DataFrame()
for rawdata_filename in tqdm(rawdata_files, desc='Process'):
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