import argparse
import numpy
import os
import pandas
import warnings

import CytoSig

from Util import read_expression

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument('-E', "--expression_file", type=str, required=False, help="Gene expression file.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/0.Tres_data/sc_cohorts/Breast.GSE156728.10x.pickle.gz')
parser.add_argument('-G', "--genesets_GMT_file", type=str, required=False, help="Background gene sets in GMT format.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/1.model_data/Tres.kegg')
parser.add_argument('-S', "--signature_name_file", type=str, required=False, help="Names of the signatures, one name in one line.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/1.model_data/signature_name_file.txt')
parser.add_argument('-D', "--output_file_directory", type=str, required=False, help="Directory for output files.",
                    default='/sibcb2/bioinformatics2/hongyuyang/dataset/Tres/Breast/GSE156728')
parser.add_argument('-O', "--output_tag", type=str, required=False, help="Prefix for output files.", default='Breast.GSE156728')
args = parser.parse_args()

def profile_geneset_signature(expression, geneset_file, name_file):
    # all gene sets
    signature = []
    fin = open(geneset_file)
    for l in fin:
        fields = l.strip().split('\t')
        s = fields[2:]
        signature.append(pandas.Series(numpy.ones(len(s)), index=s, name=fields[0]))
    fin.close()

    # gene set names
    names = []
    fin = open(name_file)
    for l in fin:
        fields = l.strip().split('\t')
        names.append(fields[0])
    fin.close()

    # pandas data frame
    # 此时signature是每个kegg通路的基因名，行为通路，列为基因
    signature = pandas.concat(signature, axis=1, join='outer', sort=False) # 将通路的所有基因都作为一列，形成一个基因数*通路数的矩阵
    signature.fillna(0, inplace=True) # 将通路中不存在的基因的值由nan变成0
    
    common = expression.index.intersection(signature.index) # 求表达矩阵和通路基因共有的基因
    signature, expression = signature.loc[common], expression.loc[common] # 只取共有的基因
    
    background = signature.mean(axis=1) # 求每个基因对所有通路的平均值
    background.name = 'study bias'
    
    X = signature.loc[:, names].mean(axis=1) # 求每个基因对name_file中的通路的平均值
    X.name = 'Proliferation'
    
    X = pandas.concat([background, X], axis=1, join='inner')
    
    # regression
    result = CytoSig.ridge_significance_test(X, expression, alpha=0, alternative="two-sided", nrand=0, verbose=1)
    
    return result[2]

expression = read_expression(args.expression_file)
result = profile_geneset_signature(expression, args.genesets_GMT_file, args.signature_name_file)
prolifertion_filename = os.path.join(args.output_file_directory, f'{args.output_tag}.Prolifertion.csv')
result.to_csv(prolifertion_filename, sep='\t')
print("Process end!")
