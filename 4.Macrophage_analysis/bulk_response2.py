import sys
import argparse
import numpy as np
import os
import pandas as pd
import warnings
from tqdm.autonotebook import tqdm

import CytoSig

warnings.filterwarnings("ignore")

# dataset_list = ["GSE14333", "GSE15654", "GSE17538", "GSE21374", "GSE28221", "GSE65218", "GSE65682", "GSE112927", "GSE33113", "GSE31595"]
dataset_list = ['Braun2020', 'E-MTAB-6270', 'GSE106128', 'GSE135222', 'GSE67501', 'GSE91061', 'GSE93157_LUSC',
                'GSE93157_nonsqNSCLC', 'Nathanson2017', 'PRJEB23709', 'PRJNA482620', 'GSE115821', 'GSE100797',
                'GSE126044', 'GSE145996', 'GSE78220', 'GSE93157_HNSC', 'GSE93157_Melanoma', 'GSE96619', 'phs000452', 'PRJEB25780']
for dataset in dataset_list:
    expression_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/5.Analysis_data/{dataset}/{dataset}.Express.tsv'
    genesets_GMT_file = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/SMaRT_geneset.txt'
    output_file_directory = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/5.Analysis_data/{dataset}'
    output_tag = f'{dataset}.bulk_response'

    def profile_macrophage_geneset_signature(expression, geneset_file, signature_name):
        # all gene sets
        signature = []
        fin = open(geneset_file)
        for l in fin:
            fields = l.strip().split('\t')
            s = fields[2:]
            signature.append(pd.Series(np.ones(len(s)), index=s, name=fields[0]))
        fin.close()

        # 此时signature是每个kegg通路的基因名，行为通路，列为基因
        signature = pd.concat(signature, axis=1, join='outer', sort=False)  # 将通路的所有基因都作为一列，形成一个基因数*通路数的矩阵
        signature.fillna(0, inplace=True)  # 将通路中不存在的基因的值由nan变成0

        common = expression.index.intersection(signature.index)  # 求表达矩阵和通路基因共有的基因
        signature, expression = signature.loc[common], expression.loc[common]  # 只取共有的基因

        background = signature.mean(axis=1)  # 求每个基因对所有通路的平均值
        background.name = f'{signature_name} study bias'

        X = signature.loc[:, signature_name].to_frame().mean(axis=1)  # 求每个基因对name_file中的通路的平均值
        X.name = signature_name

        X = pd.concat([background, X], axis=1, join='inner')

        # regression
        result = CytoSig.ridge_significance_test(X, expression, alpha=0, alternative="two-sided", nrand=0, verbose=1)

        return result[2]

    expression = pd.read_csv(expression_path, index_col=0, header=0, delimiter='\t')
    expression.set_index('GENE_SYMBOL', inplace=True)
    expression = expression.dropna()

    # normalized
    expression *= 1E5 / expression.sum()
    expression = np.log2(expression + 1)
    expression = expression.subtract(expression.mean(axis=1), axis=0)

    result = []
    for signature_name in ['SMART_C13', 'SMART_C3']:
        sub_result = profile_macrophage_geneset_signature(expression, genesets_GMT_file, signature_name)
        result.append(sub_result)
    result = pd.concat(result)

    response_filename = os.path.join(output_file_directory, f'{output_tag}.csv')
    result.to_csv(response_filename, sep='\t')
    print(f"{dataset} Process end!")
