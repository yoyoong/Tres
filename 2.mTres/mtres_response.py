import sys
import argparse
import numpy
import os

import numpy as np
import pandas as pd
import warnings
from tqdm.autonotebook import tqdm
import random

import CytoSig

from Util import read_expression

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_path", type=str, required=False, help="Gene expression file or directory.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/1.new_gem_data/AEL_GSE142213')
parser.add_argument('-G', "--genesets_GMT_file", type=str, required=False, help="Background gene sets in GMT format.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/Tres_kegg.NK.txt')
parser.add_argument('-S', "--signature_name_file", type=str, required=False, help="Names of the signatures, one name in one line.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/signature_name_file.txt')
parser.add_argument('-CT', "--celltype", type=str, default='NK', required=False, help="cell type")
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch_data/3-4.NK_response')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='AEL_GSE142213')
args = parser.parse_args()

expression_path = args.expression_path
genesets_GMT_file = args.genesets_GMT_file
signature_name_file = args.signature_name_file
celltype = args.celltype
output_file_directory = args.output_file_directory
output_tag = args.output_tag

all_signature_name = pd.read_csv(signature_name_file, index_col=0, header=0)
geneset_name = list(pd.Series(all_signature_name.loc[celltype]['signature']))
real_celltype = celltype
if "_" in real_celltype:
    real_celltype = real_celltype.split("_")[0]

if real_celltype == 'Macrophage':
    celltype_in_column = 'Mono/Macro'
    celltype_in_file = 'Mono_Macro'
else:
    celltype_in_column = real_celltype
    celltype_in_file = real_celltype

print("Process start!")
def profile_geneset_signature(expression, geneset_file, geneset_name, signature_name):
    # all gene sets
    signature = []
    fin = open(geneset_file)
    for l in fin:
        fields = l.strip().split('\t')
        s = fields[2:]
        signature.append(pd.Series(numpy.ones(len(s)), index=s, name=fields[0]))
    fin.close()

    # 此时signature是每个kegg通路的基因名，行为通路，列为基因
    signature = pd.concat(signature, axis=1, join='outer', sort=False) # 将通路的所有基因都作为一列，形成一个基因数*通路数的矩阵
    signature.fillna(0, inplace=True) # 将通路中不存在的基因的值由nan变成0
    
    common = expression.index.intersection(signature.index) # 求表达矩阵和通路基因共有的基因
    signature, expression = signature.loc[common], expression.loc[common] # 只取共有的基因
    
    background = signature.mean(axis=1) # 求每个基因对所有通路的平均值
    background.name = 'study bias'
    
    X = signature.loc[:, geneset_name].mean(axis=1) # 求每个基因对name_file中的通路的平均值
    X.name = signature_name
    
    X = pd.concat([background, X], axis=1, join='inner')
    
    # regression
    result = CytoSig.ridge_significance_test(X, expression, alpha=0, alternative="two-sided", nrand=0, verbose=1)
    
    return result[2]


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


result = pd.DataFrame()
# get the expression by celltype
expression_list = sorted(os.listdir(expression_path))
celltype_list = [v.split('.')[1] for v in expression_list]
if celltype_in_file not in celltype_list:
    print(f"This dataset has not {celltype_in_file} celltype.")
    sys.exit(1)
tag = expression_path.split('/')[-1]
expression = pd.read_csv(os.path.join(expression_path, f'{tag}.{celltype_in_file}.csv'), index_col=0)

if celltype == 'CD8T':
    result = profile_geneset_signature(expression, genesets_GMT_file, geneset_name, 'Proliferation')
elif celltype == 'Macrophage':
    c13_result = profile_macrophage_geneset_signature(expression, genesets_GMT_file, 'SMART_C13')
    c3_result = profile_macrophage_geneset_signature(expression, genesets_GMT_file, 'SMART_C3')
    result = pd.concat([c13_result, c3_result])
    new_row = result.loc["SMART_C13"] - result.loc["SMART_C3"]
    new_row.name = 'Polarization'
    result = pd.concat([result, new_row.to_frame().T])
elif celltype == 'Neutrophils':
    result = profile_geneset_signature(expression, genesets_GMT_file, geneset_name, 'Neut_IFN-15')
elif celltype == 'NK':
    result = profile_geneset_signature(expression, genesets_GMT_file, geneset_name, 'NK_signature')
elif celltype == 'NK_act':
    result = profile_geneset_signature(expression, genesets_GMT_file, geneset_name, 'NK_act_signature')


response_filename = os.path.join(output_file_directory, f'{output_tag}.csv')
result.to_csv(response_filename, sep='\t')
print("Process end!")
