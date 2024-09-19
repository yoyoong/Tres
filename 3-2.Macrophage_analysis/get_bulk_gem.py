import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr



dataset_list = ["GSE14333", "GSE15654", "GSE17538", "GSE21374", "GSE28221", "GSE65218", "GSE65682", "GSE112927", "GSE33113", "GSE31595"]
for dataset in dataset_list:
    output_dir = f'/sibcb1/bioinformatics/hongyuyang/dataset/Tres/4.Macrophage_analysis/{dataset}'

    if dataset == 'GSE15654':
        series_matrix_path = os.path.join(output_dir, 'GSE15654_series_matrix.txt')
        gene_map_path = os.path.join(output_dir, 'GPL8432-11703.txt')
        gem_df = pd.read_csv(series_matrix_path, sep='\t', index_col=0, header=0, comment="!")
        gene_map_df = pd.read_csv(gene_map_path, sep='\t', index_col=0, header=0, comment="#")

        # convert ILMN id to gene name
        gem_df = gem_df.assign(gene_name=lambda x: gene_map_df.loc[x.index]['ILMN_Gene'])
        gem_df = gem_df.groupby('gene_name').max()  # if a gene has mutiple value, use the max

    elif (dataset == 'GSE21374' or dataset == 'GSE14333' or dataset == 'GSE17538'
          or dataset == 'GSE33113' or dataset == 'GSE31595'):
        series_matrix_path = os.path.join(output_dir, f'{dataset}_series_matrix.txt')
        gene_map_path = os.path.join(output_dir, 'GPL570-55999.txt')
        gem_df = pd.read_csv(series_matrix_path, sep='\t', index_col=0, header=0, comment="!")
        gene_map_df = pd.read_csv(gene_map_path, sep='\t', index_col=0, header=0, comment="#")

        def get_genename(x):
            symbol_dict = dict(gene_map_df.loc[x.index]['Gene Symbol'])
            for key, symbol in symbol_dict.items():
                if str(symbol).find('///') >= 0:  # if a gene has mutiple value, use the first
                    symbol_dict[key] = symbol.split('///')[0].strip()
            return pd.Series(symbol_dict)

        gem_df = gem_df.assign(gene_name=lambda x: get_genename(x))
        gem_df = gem_df.groupby('gene_name').max()  # if a gene has mutiple value, use the max

    elif dataset == 'GSE65682':
        series_matrix_path = os.path.join(output_dir, 'GSE65682_series_matrix.txt')
        gene_map_path = os.path.join(output_dir, 'GPL13667-15572.txt')

        gem_df = pd.read_csv(series_matrix_path, sep='\t', index_col=0, header=0, comment="!")
        gene_map_df = pd.read_csv(gene_map_path, sep='\t', index_col=0, header=0, comment="#")

        def get_genename(x):
            symbol_dict = dict(gene_map_df.loc[x.index]['Gene Symbol'])
            for key, symbol in symbol_dict.items():
                if str(symbol).find('///') >= 0:  # if a gene has mutiple value, use the first
                    symbol_dict[key] = symbol.split('///')[0].strip()
            return pd.Series(symbol_dict)

        gem_df = gem_df.assign(gene_name=lambda x: get_genename(x))
        gem_df = gem_df.groupby('gene_name').max()  # if a gene has mutiple value, use the max

    elif dataset == 'GSE28221':
        gem_path = os.path.join(output_dir, 'GSE28221_GeneID_matched_gene_expression_dataset_of_IPF_cohorts.txt')
        gem_df = pd.read_csv(gem_path, sep='\t', index_col=1, header=0)
        gem_df = gem_df.drop('Gene ID', axis=1)
        # flag = [item for item in gem_df.columns if item.startswith('IPF')]
        # gem_df = gem_df[flag]

    elif dataset == 'GSE112927':
        gem_path = os.path.join(output_dir, 'GSE112927_RAW')
        gem_filelist = sorted(os.listdir(gem_path))
        gem_df = pd.DataFrame()
        for gem_file in gem_filelist:
            sub_gem_df = pd.read_csv(os.path.join(gem_path, gem_file), sep='\t', index_col=0, header=0)
            sub_gem_df = sub_gem_df.groupby(sub_gem_df.index).agg('max')
            if len(gem_df) == 0:
                gem_df = sub_gem_df
            else:
                gem_df = gem_df.join(sub_gem_df, how='inner')

    elif dataset == 'GSE65218':
        series_matrix_path = os.path.join(output_dir, 'GSE65218_series_matrix.txt')
        gene_map_path = os.path.join(output_dir, 'GPL159888.txt')
        gem_df = pd.read_csv(series_matrix_path, sep='\t', index_col=0, header=0, comment="!")
        gene_map_df = pd.read_csv(gene_map_path, sep='\t', index_col=0, header=0, comment="#")

        # convert ILMN id to gene name
        gem_df = gem_df.assign(gene_name=lambda x: gene_map_df.loc[x.index]['GENE SYMBOL'])
        gem_df = gem_df.groupby('gene_name').max()  # if a gene has mutiple value, use the max

    if dataset == 'GSE112927' or dataset == 'GSE14333':
        gem_df = gem_df.subtract(gem_df.mean(axis=1), axis=0)
    else:
        if len(gem_df) > len(set(gem_df.index)):
            gem_df = gem_df.groupby(gem_df.index).agg('max')
        # normalized
        gem_df *= 1E5 / gem_df.sum()
        gem_df = np.log2(gem_df + 1)
        gem_df = gem_df.subtract(gem_df.mean(axis=1), axis=0)
    gem_df.to_csv(os.path.join(output_dir, f'{dataset}.gem.csv'))
    print(f"{dataset} process end")
