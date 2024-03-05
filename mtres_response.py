import argparse
import numpy
import os
import pandas as pd
import warnings
from tqdm.autonotebook import tqdm

import CytoSig

from Util import read_expression

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_path", type=str, required=False, help="Gene expression file or directory.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/1.gem_data/CRC_GSE120909_mouse_aPD1.csv')
parser.add_argument('-G', "--genesets_GMT_file", type=str, required=False, help="Background gene sets in GMT format.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/Tres.kegg')
parser.add_argument('-S', "--signature_name_file", type=str, required=False, help="Names of the signatures, one name in one line.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.model_file/signature_name_file.txt')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/2.tisch2_data/2-1.Prolifertion')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='CRC_GSE120909_mouse_aPD1')
args = parser.parse_args()

expression_path = args.expression_path
genesets_GMT_file = args.genesets_GMT_file
signature_name_file = args.signature_name_file
output_file_directory = args.output_file_directory
output_tag = args.output_tag

def profile_geneset_signature(expression, geneset_file, name_file):
    # all gene sets
    signature = []
    fin = open(geneset_file)
    for l in fin:
        fields = l.strip().split('\t')
        s = fields[2:]
        signature.append(pd.Series(numpy.ones(len(s)), index=s, name=fields[0]))
    fin.close()

    # gene set names
    names = []
    fin = open(name_file)
    for l in fin:
        fields = l.strip().split('\t')
        names.append(fields[0])
    fin.close()

    # 此时signature是每个kegg通路的基因名，行为通路，列为基因
    signature = pd.concat(signature, axis=1, join='outer', sort=False) # 将通路的所有基因都作为一列，形成一个基因数*通路数的矩阵
    signature.fillna(0, inplace=True) # 将通路中不存在的基因的值由nan变成0
    
    common = expression.index.intersection(signature.index) # 求表达矩阵和通路基因共有的基因
    signature, expression = signature.loc[common], expression.loc[common] # 只取共有的基因
    
    background = signature.mean(axis=1) # 求每个基因对所有通路的平均值
    background.name = 'study bias'
    
    X = signature.loc[:, names].mean(axis=1) # 求每个基因对name_file中的通路的平均值
    X.name = 'Proliferation'
    
    X = pd.concat([background, X], axis=1, join='inner')
    
    # regression
    result = CytoSig.ridge_significance_test(X, expression, alpha=0, alternative="two-sided", nrand=0, verbose=1)
    
    return result[2]

result = pd.DataFrame()
if not os.path.isdir(expression_path):
    expression = read_expression(expression_path)
    result = profile_geneset_signature(expression, genesets_GMT_file, signature_name_file)
else:
    expression_list = sorted(os.listdir(expression_path))
    for expression_file in tqdm(expression_list, desc="Calculate"):
        expression = read_expression(os.path.join(expression_path, expression_file))
        sub_result = profile_geneset_signature(expression, genesets_GMT_file, signature_name_file)
        result = pd.concat([result, sub_result], axis=1)
        print(f"{expression_file} calculate end.")

prolifertion_filename = os.path.join(output_file_directory, f'{output_tag}.csv')
result.to_csv(prolifertion_filename, sep='\t')
print("Process end!")
