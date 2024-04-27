import sys
import argparse
import numpy
import os
import pandas as pd
import warnings
from tqdm.autonotebook import tqdm

import CytoSig

warnings.filterwarnings("ignore")

dataset_list = ["GSE14333", "GSE15654", "GSE17538", "GSE21374", "GSE28221", "GSE65218", "GSE65682", "GSE112927", "GSE33113", "GSE31595"]
for dataset in dataset_list:
    expression_path = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/4.Macrophage_analysis/{dataset}/{dataset}.gem.csv'
    genesets_GMT_file = '/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/SMaRT_geneset.txt'
    output_file_directory = f'/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/4.Macrophage_analysis/{dataset}'
    output_tag = f'{dataset}.bulk_response'

    def profile_macrophage_geneset_signature(expression, geneset_file, signature_name):
        # all gene sets
        signature = []
        fin = open(geneset_file)
        for l in fin:
            fields = l.strip().split('\t')
            s = fields[2:]
            signature.append(pd.Series(numpy.ones(len(s)), index=s, name=fields[0]))
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

    expression = pd.read_csv(expression_path, index_col=0, header=0)
    result = []
    for signature_name in ['SMART_C13', 'SMART_C3', 'SMART_C14']:
        sub_result = profile_macrophage_geneset_signature(expression, genesets_GMT_file, signature_name)
        result.append(sub_result)
    result = pd.concat(result)
    # c13_result = profile_macrophage_geneset_signature(expression, genesets_GMT_file, 'SMART_C13')
    # c3_result = profile_macrophage_geneset_signature(expression, genesets_GMT_file, 'SMART_C3')
    # result = pd.concat([c13_result, c3_result])
    # new_row = result.loc["SMART_C13"] - result.loc["SMART_C3"]
    # new_row.name = 'Polarization'
    # result = pd.concat([result, new_row.to_frame().T])

    response_filename = os.path.join(output_file_directory, f'{output_tag}.csv')
    result.to_csv(response_filename, sep='\t')
    print(f"{dataset} Process end!")
